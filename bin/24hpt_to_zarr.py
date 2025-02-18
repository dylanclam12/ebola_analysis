import shutil
from pathlib import Path

from spatialdata_io import merscope


def main():
    projdir = Path("/mnt/d/ebola")
    sample_path = Path(f"{projdir}/data/raw/202406121440_HEK293T-VS219-24HPT-CompA-S1_VMSC12502/region_0")
    zarr_path = Path(f"{projdir}/data/24hpt.zarr")

    sdata = merscope(sample_path, z_layers=0)

    if zarr_path.exists():
        shutil.rmtree(zarr_path)
    sdata.write(zarr_path)

    print(f"Saved merscope data:\n{sdata}")


if __name__ == "__main__":
    main()
