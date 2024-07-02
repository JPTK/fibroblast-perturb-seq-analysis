"""
Lara mouse fibroblast resequencing data 2023
Unified run script for Rd1 and Rd3 assays
"""

from typing import Any
from paths import DIR_DATA_GENERATED, DIR_DATA_RNA
from pathlib import Path
import yaml
from parse import parse  # Adjusted to a unified parse script

# Load configurations from the YAML file
# with open("config.yaml", "r") as file:
#     config = yaml.safe_load(file)

LOCAL_DATA_DIR = "/home/jovyan/data/lara_mouse_fibroblast_reseq_2023/reseqFibro2/"
ASSAYS = ["Rd1", "Rd3"]


def load_config(assay: str) -> tuple[list[Path], list[Path]]:
    """
    Loads configuration for a specific assay from the config YAML file.
    It generates file paths for call files and count files based on the processed folders
    specified for each assay in the configuration.

    :param assay: The assay identifier (e.g., 'Rd1', 'Rd3') for which to load the configuration.
    :type assay: str
    :return: A tuple containing two lists: call_files, and count_files.
             - call_files: List of paths to CRISPR call files.
             - count_files: List of paths to count data files.
    :rtype: tuple
    """
    assay_to_prefixes = {
        "Rd1":["Activated_OP1L_NM_NA_Rep1","Activated_OP1L_NM_NA_Rep2","Quiescent_OP1L_NM_NA_Rep1","Quiescent_OP1L_NM_NA_Rep2"],
        "Rd3":["exVivo_OP2_IL1b_1","exVivo_OP2_IL1b_2","exVivo_OP2_resting_1","exVivo_OP2_resting_2","exVivo_OP2_TGFb_1","exVivo_OP2_TGFb_2"],
        }
    prefixes = assay_to_prefixes[assay]
    call_files = [
        Path(DIR_DATA_RNA / f"{prefix}_protospacer_calls_per_cell.csv.gz")
        for prefix in prefixes
    ]
    count_files = [
        Path(DIR_DATA_RNA / f"{prefix}_filtered_feature_bc_matrix.h5")
        for prefix in prefixes
    ]
    return call_files, count_files




def create_h5ad(assay:str, call_files:list[Path], count_files:list[Path]) -> None:
    # generate parsed files
    
    parsed_file = DIR_DATA_GENERATED / f"{assay}_Data_post_sgRNA_rescue-reseq.h5ad"

    # make input data
    input_dict = {
        "call_files": {call_file.name: call_file for call_file in call_files},
        "count_files": {count_file.name: count_file for count_file in count_files}
    }

    # parse
    adata = parse(input_dict, assay)
    # save anndata
    adata.write_h5ad(parsed_file, compression="lzf")


def run(assay: str) -> None:
    call_files, count_files = load_config(assay)
    create_h5ad(assay, call_files, count_files)


if __name__ == "__main__":
    for assay in ASSAYS:
        run(assay)