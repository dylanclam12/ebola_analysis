import shutil
from pathlib import Path

from spatialdata_io import merscope


def main():
    projdir = Path("/mnt/d/ebola")
    sample_path = Path(
        f"{projdir}/data/raw/202406121421_HEK293T-VS219-48HPT-CompA-S1-JM_VMSC07201/region_0"
    )
    zarr_path = Path(f"{projdir}/data/48hpt.zarr")

    sdata = merscope(sample_path, z_layers=0)

    if zarr_path.exists():
        shutil.rmtree(zarr_path)
    sdata.write(zarr_path)

    print(f"Saved merscope data:\n{sdata}")


if __name__ == "__main__":
    main()
