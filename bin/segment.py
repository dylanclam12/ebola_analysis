import argparse
import tempfile
from typing import Literal

import cv2
import dask
import numpy as np
import spatialdata as sd
import xarray as xr
from dask import delayed
from skimage import rank
from skimage.exposure import adjust_gamma
from skimage.filters import gaussian, threshold_multiotsu
from skimage.morphology import disk
import geopandas as gpd

import sopa


def segment_protein_channel(
    sdata,
    protein_channel: str,
    num_classes: int = 3,
    gamma: float = 1,
    sigma: float = 0,
    chunks: int = 5000,
    depth: int = 100,
):
    """
    Segment protein channel using multi-Otsu thresholding.

    Parameters:
    -----------
    sdata : SpatialData
        The SpatialData object containing the image data.
    protein_channel : str
        The name of the protein channel to segment.
    num_classes : int, optional
        The number of classes for multi-Otsu thresholding. Default is 3.
    chunks : int, optional
        The chunk size for processing. Default is 5000.
    depth : int, optional
        The depth parameter for sd.map_raster. Default is 100.

    Returns:
    --------
    None. The segmented image is added to the sdata object.
    """
    protein_input = sdata["image"]["scale0"].chunk(chunks).sel(c=protein_channel).image

    def apply_multiotsu(img, num_classes=3, gamma=1, sigma=0):
        if gamma != 1:
            new_img = adjust_gamma(img, gamma)

        if sigma > 0:
            new_img = gaussian(new_img, sigma=sigma)

        try:
            thresholds = threshold_multiotsu(new_img, classes=num_classes)
            regions = np.digitize(new_img, bins=thresholds)
            return regions
        except ValueError:
            return np.zeros_like(img)

    segmented_da = sd.map_raster(
        protein_input,
        apply_multiotsu,
        func_kwargs={"num_classes": num_classes, "gamma": gamma, "sigma": sigma},
        dims=("y", "x"),
        depth=depth,
        trim=True,
        meta=np.array(()),
    )

    label_name = f"{protein_channel.lower()}_multiotsu"
    shape_name = f"{protein_channel.lower()}_shapes"
    sdata[label_name] = sd.models.Labels2DModel.parse(
        (segmented_da.compute() == (num_classes - 1)).astype(int)[0]
    )
    protein_shapes = gpd.GeoDataFrame(
        sd.to_polygons(sdata[label_name]).explode(ignore_index=True).geometry
    )
    sdata[shape_name] = sd.models.Shapes2DModel.parse(protein_shapes)

    if f"shapes/{shape_name}" in sdata.element_paths_on_disk():
        sdata.delete_element_from_disk(shape_name)
    sdata.write_element(shape_name)


def preprocess_image(
    img,
    kernel_size=(5, 5),
    iterations=3,
    normalize=False,
    footprint_radius=100,
):
    kernel = np.ones(kernel_size, np.uint8)
    new_img = cv2.erode(img, kernel, iterations=iterations)

    if normalize:
        footprint = disk(footprint_radius)
        new_img = rank.equalize(new_img, footprint=footprint)

    return new_img


@delayed
def process_patch_delayed(segmentation_type, patch_dir, patch_index):
    segmentation_type.write_patch_cells(patch_dir, patch_index)


def preprocess_image_channel(
    sdata,
    chunks=5000,
    depth=100,
    channel="DAPI",
    normalize=False,
    footprint_radius=100,
    num_workers=8,
):
    dapi_input = sdata["image"]["scale0"].chunk(chunks).sel(c=channel).image

    # DAPI preprocessing
    img_da = sd.map_raster(
        dapi_input,
        preprocess_image,
        func_kwargs={
            "normalize": normalize,
            "footprint_radius": footprint_radius,
        },
        dims=("y", "x"),
        depth=depth,
        trim=True,
        meta=np.array(()),
    )

    out_name = f"{str(channel).lower()}_processed"
    sdata[out_name] = sd.models.Image2DModel.parse(
        img_da.compute(num_workers=num_workers).expand_dims(dim={"c": 1})
    )

    if f"images/{out_name}" in sdata.element_paths_on_disk():
        sdata.delete_element_from_disk(out_name)

    sdata.write_element(out_name)


def segment_image(
    sdata,
    mode: Literal["cell", "nucleus"],
    patch_width: int,
    patch_overlap: int,
    diameter: int,
    flow_threshold: float,
    num_workers: int,
    nucleus_img_key: str,
    cell_img_key: str = None,
):
    # Define keys for each mode
    if mode == "cell":
        img_keys = [cell_img_key, nucleus_img_key]  # Order matterss
        model = "cyto3"
        shape_key = "cell_boundaries"
    elif mode == "nucleus":
        img_keys = [nucleus_img_key]
        model = "nuclei"
        shape_key = "nucleus_boundaries"

    # Preparation of input for Cellpose and SOPA patch creation
    cellpose_img = xr.concat(
        [sdata[img_key] for img_key in img_keys],
        dim="c",
    )
    sdata["cellpose_img"] = sd.models.Image2DModel.parse(
        cellpose_img, dims=["c", "y", "x"], c_coords=img_keys
    )
    patches = sopa.segmentation.Patches2D(
        sdata, "cellpose_img", patch_width=patch_width, patch_overlap=patch_overlap
    )

    seg_method = sopa.segmentation.methods.cellpose_patch(
        model_type=model,
        diameter=diameter,
        channels=img_keys,
        flow_threshold=flow_threshold,
        cellpose_model_kwargs={"gpu": True},
    )
    segmentation_type = sopa.segmentation.StainingSegmentation(
        sdata,
        seg_method,
        img_keys,
        image_key="cellpose_img",
        clip_limit=0,
        gaussian_sigma=0,
    )

    # Parallel processing
    with tempfile.TemporaryDirectory() as patch_dir:
        tasks = [
            process_patch_delayed(segmentation_type, patch_dir, i)
            for i in range(len(patches))
        ]
        dask.compute(*tasks, num_workers=num_workers)

        # Post-processing
        shapes = sopa.segmentation.StainingSegmentation.read_patches_cells(patch_dir)
        shapes = sopa.segmentation.shapes.solve_conflicts(shapes)

        sopa.segmentation.StainingSegmentation.add_shapes(
            sdata, shapes, "cellpose_img", shape_key
        )

    if f"shapes/{shape_key}" in sdata.element_paths_on_disk():
        sdata.delete_element_from_disk(shape_key)
    sdata.write_element(shape_key)


def segment_dataset(
    zarr_path: str,
    nucleus_diameter: int = 100,
    nucleus_flow_threshold: float = 0.9,
    cell_diameter: int = 130,
    cell_flow_threshold: float = 0.9,
    num_workers: int = 1,
):
    sdata = sd.read_zarr(zarr_path)

    # Number of workers for segmentation, heuristic 1/2 of total workers
    seg_num_workers = int(num_workers / 2) if num_workers > 1 else 1

    # Segment nuclei
    preprocess_image_channel(
        sdata,
        channel="DAPI",
        normalize=False,
        chunks=5000,
        depth=50,
        num_workers=num_workers,
    )
    segment_image(
        sdata,
        mode="nucleus",
        nucleus_img_key="dapi_processed",
        cell_img_key=None,
        patch_width=3000,
        patch_overlap=200,
        diameter=nucleus_diameter,
        flow_threshold=nucleus_flow_threshold,
        num_workers=seg_num_workers,
    )

    # Segment cells
    preprocess_image_channel(
        sdata,
        channel="PolyT",
        normalize=True,
        footprint_radius=100,
        chunks=5000,
        depth=200,
        num_workers=num_workers,
    )
    segment_image(
        sdata,
        mode="cell",
        nucleus_img_key="dapi_processed",
        cell_img_key="polyt_processed",
        patch_width=3000,
        patch_overlap=200,
        diameter=cell_diameter,
        flow_threshold=cell_flow_threshold,
        num_workers=seg_num_workers,
    )
    print(f"Segmentation complete. Updated data saved to: {zarr_path}")


def main():
    parser = argparse.ArgumentParser(description="Segment and analyze spatial data")
    parser.add_argument("--input", required=True, help="Input Zarr file path")
    args = parser.parse_args()

    segment_dataset(args.input)


if __name__ == "__main__":
    main()
