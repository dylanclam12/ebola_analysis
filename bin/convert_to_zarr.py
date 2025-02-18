import argparse
import shutil
from pathlib import Path

from spatialdata_io import merscope


def convert_to_zarr(input_dir, output_file):
    input_path = Path(input_dir)
    output_path = Path(output_file)

    sdata = merscope(input_path, z_layers=0)

    if output_path.exists():
        shutil.rmtree(output_path)
    sdata.write(output_path)

    print(f"Saved merscope data:\n{sdata}")


def main():
    parser = argparse.ArgumentParser(description="Convert MERSCOPE data to Zarr format")
    parser.add_argument(
        "--input", required=True, help="Input directory containing MERSCOPE data"
    )
    parser.add_argument("--output", required=True, help="Output Zarr file path")
    args = parser.parse_args()

    convert_to_zarr(args.input, args.output)


if __name__ == "__main__":
    main()
