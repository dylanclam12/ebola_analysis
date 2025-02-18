# Getting Started

This repository reproduces segmentation and subcellular localization analysis of VP35 condensates in a HEK cell model, profiled with MERFISH (Vizgen).

## Organization

### Main Analysis Notebooks
 
1. `*_segmentation_run_all_{24,48}hpt.ipynb`: Perform image preprocessing and segmentation on the 24hpt and 48hpt datasets respectively.

2. `*_condensate_stats_{24,48}hpt.ipynb`: Perform quality control and filtering, compute condensate, and cytoplasmic summary statistics and gene-specific statistics.

### Supporting Notebooks
- `*_segmentation_preprocess.ipynb`: Initial testing of image preprocessing steps.
- `*_segmentation_routine_{24,48}hpt.ipynb`: Testing on a small subset of the dataset to tune parameters.
- `*_condensate_diff_density_{24,48}hpt.ipynb` **[Experimental]**: Single-cell analysis of differential density of condensates between 24hpt and 48hpt. The accuracy for this analysis needs to be verified.

## Installation

1. Clone this repository.
2. Create a virtual environment using `conda` or `venv`. Install Python>=3.10.
3. Install the required packages using the provided `requirements.txt` file.
    ```bash
    pip install -r requirements.txt
    ```

To fully utilize the GPU for segmentation, check proper installation for the following dependencies:
- Cellpose: https://cellpose.readthedocs.io/en/latest/installation.html
- Sopa: https://gustaveroussy.github.io/sopa/getting_started/