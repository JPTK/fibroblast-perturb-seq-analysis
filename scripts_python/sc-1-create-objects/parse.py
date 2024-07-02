from typing import Any

import anndata as ad
import pandas as pd
import scanpy as sc
import yaml
from pathlib import Path 

# # Load configurations from the YAML file
# with open("config.yaml", "r") as file:
#     config = yaml.safe_load(file)

STUDY_CODE = "LaraMouseFibroblastReseq2023"



def assign_unique_barcode(adata: ad.AnnData, study_code: str) -> None:
    """
    Generates unique barcodes by adding the paper/experiment context to the existing adata.obs_names

    :param adata: adata object to process with desired cell ids as index
    :type adata: adata object
    :param study_code: unique string to be appended to cell barcode
    :type study_code: str
    """

    unique_ids = [x + "_%s" % study_code for x in adata.obs_names]

    adata.obs_names = pd.Index(unique_ids, name="cell_id")

    return


def get_assay_config(assay: str) -> None | dict[str, Any]:
    """
    Retrieves the configuration for a specific assay from the loaded YAML configuration.

    This function looks up the assay-specific configuration (like processed folders, cell suffix mappings)
    within the global configuration dictionary.

    :param assay: The assay identifier for which to fetch the configuration.
    :type assay: str
    :return: A dictionary containing the assay-specific configuration if found;
             None if the assay is not present in the configuration.
    :rtype: None | dict[str, Any]
    """
    result: None | dict[str, Any] = config.get(assay)
    return result


def load_call_files(key_to_file_path: dict[str,Path]) -> pd.DataFrame:
    """
    Load guide calling csv files in a dataframe.

    :param files: list of RawFiles corresponding to csv files containing guide calling information.
    :return: DataFrame containing information about all guide calls
    """
    guide_counts = []
    for key, file_path in key_to_file_path.items():
        df = pd.read_csv(file_path, index_col="cell_barcode").drop(
            columns=["num_features"]
        )
        df = df.apply(lambda col: col.str.split("|"))
        df = df.apply(pd.Series.explode)
        df = df.pivot(columns="feature_call", values="num_umis").fillna(0)
        df.index = df.index + "-" + str(key).split("/")[1]
        guide_counts.append(df)
    guide_count = pd.concat(guide_counts)
    guide_count = guide_count.fillna(0)
    return guide_count.astype(int)


def combine_columns_by_prefix(df: pd.DataFrame) -> pd.DataFrame:
    """
    Combine columns in a DataFrame based on the prefix of their names.

    This function groups columns that have the same prefix (before the first "_")
    and sums their values. The resulting DataFrame will have combined counts
    for columns with the same prefix.

    :param df: DataFrame with columns to be combined
    :type df: pd.DataFrame
    :return: DataFrame with combined columns
    :rtype: pd.DataFrame
    """

    df_T = df.T
    combined = df_T.groupby(df_T.index.str.split("_").str[0]).sum().T

    return combined


def filter_on_moi(
    guide_count_matrix: pd.DataFrame, threshold: float = 1.75
) -> pd.DataFrame:
    """
    Filter cells of MOI lower or greater than 1.

    This is based off Jake TK's analysis.
    It first filters out some guide counts based on frequency of guides and counts.
    The argument is that low counts paired with low frequencies imply that the
    guide is unlikely to have been integrated in the cell. Following this thresholding,
    the function only keeps cells where the multiplicity of targeting (MOT) is 1. The multiplicity of
    infection (MOI) might be > 1 if all guides target the same gene (thus MOT = 1). Furthermore, cells where MOT = 2
    but one of the guides is a non-targeting control are also kept.
    :param guide_count_matrix: count matrix of size number of cells by number of guides
    :type guide_count_matrix: pd.DataFrame
    :param threshold: threshold applied to count matrix, defaults to 1.75
    :type threshold: float, optional
    :return: dataframe with cell id index and target sgRNA column
    :rtype: pd.DataFrame
    """
    guide_freq_matrix = (
        guide_count_matrix / guide_count_matrix.sum(axis=1).values[:, None]
    )
    guide_count_matrix[
        (guide_freq_matrix <= threshold / 10) & (guide_count_matrix <= 10**threshold)
    ] = 0
    guide_count_matrix_by_target = combine_columns_by_prefix(guide_count_matrix)
    mot = guide_count_matrix_by_target.gt(0).sum(axis=1)
    df1 = guide_count_matrix[mot == 1]
    obs1 = df1.apply(lambda row: "|".join(row.index[row != 0]), axis=1).to_frame(
        name="sgRNA"
    )
    df2 = guide_count_matrix[(mot == 2) & (guide_count_matrix_by_target["NTC"] > 0)]
    obs2 = df2.apply(lambda row: "|".join(row.index[row != 0]), axis=1).to_frame(
        name="sgRNA"
    )
    obs = pd.concat([obs1, obs2])
    return obs


def load_count_data(key_to_filepath: dict[str,Path]) -> ad.AnnData:
    """
    Load scRNA-seq count data in an AnnData.

    :param files: list of paths to files containing the count data
    :return: object containing scRNA-seq count data.
    :rtype: anndata.AnnData
    """
    adatas = []
    for key, file_path in key_to_filepath.items():
        adata = sc.read_10x_h5(file_path)
        stamp = str(key).split("/")[1]
        adata.obs.index = adata.obs.index + "-" + stamp
        adata.var.index = adata.var.index.str.lower()
        adata.var_names_make_unique()  # ensured that none of the targets are affected by this
        adatas.append(adata)
    adata = ad.concat(adatas, axis=0, merge="same")
    return adata


def filter_cells_from_counts(
    counts_adata: ad.AnnData, cell_df: pd.DataFrame
) -> ad.AnnData:
    """
    Filter count adata based on cell dataframe and include cell metadata.

    :param counts_adata: scRNA-seq count data
    :type counts_adata: anndata.AnnData
    :param cell_df: cell metadata
    :type cell_df: pd.DataFrame
    :return: filtered scRNA-seq count data
    :rtype: anndata.AnnData
    """
    adata = counts_adata[cell_df.index]
    adata.obs = pd.concat([adata.obs, cell_df], axis=1)
    return adata


def simplify_target_value(value: str) -> str:
    """Simplifies target string to remove duplicates

    param value: target value
    return: string simplified target value
    """

    unique_values = list(set(value.split("|")))
    return "|".join(unique_values)


def standardise_obs(
    adata: ad.AnnData,
    cell_suffix_to_condition: dict[str, str],
    cell_suffix_to_replicate: dict[str, str],
) -> ad.AnnData:
    """Reads in the anndata object and standardises cell-level metadata.

    :param adata: anndata object containing cell metadata in .obs
    :type adata: anndata.Anndata
    :param cell_suffix_to_condition: dictionary cell suffixes to condition
    :param cell_suffix_to_replicate: dictionary linking cell suffixes to replicates
    :return: adata object with standardised cell metadata in .obs
    :rtype: anndata.Anndata
    """

    adata.obs["misc:umi_count"] = adata.X.toarray().sum(1)

    # Add condition column
    adata.obs["condition"] = (
        adata.obs.index.str.split("-").str[2].map(cell_suffix_to_condition)
    )
    # Add replicate column
    adata.obs["replicate"] = (
        adata.obs.index.str.split("-").str[2].map(cell_suffix_to_replicate)
    )
    # Add CRISPR metadata
    adata.obs["target"] = adata.obs.sgRNA.apply(
        lambda x: "|".join([s.split("_")[0] for s in x.split("|")])
    )
    adata.obs["target"] = adata.obs["target"].str.replace("NTC", "non-targeting")
    adata.obs["misc::target_simple"] = (
        adata.obs["target"]
        .apply(simplify_target_value)
        .str.replace("non-targeting|", "", regex=False)
        .str.replace("|non-targeting", "", regex=False)
    )
    adata.obs["is_nt"] = [
        True if x.split("_")[0] == "non-targeting" else False
        for x in adata.obs["target"]
    ]

    adata.obs["moi"] = adata.obs["sgRNA"].str.split("|").str.len()

    # rename cell ids
    assign_unique_barcode(adata, f"{STUDY_CODE}")

    return adata


def standardise_var(adata: ad.AnnData) -> ad.AnnData:
    """Reads in the anndata object and standardises gene-level metadata.

    :param adata: anndata object containing gene metadata in .var
    :type adata: anndata.Anndata
    :return: adata object with standardised gene metadata in .var
    :rtype: anndata.Anndata
    """

    # Retain only relevant gene id column
    # Capitalize first letters of gene names according to mouse gene name conventions
    adata.var = adata.var.index.to_frame().applymap(lambda x: x.capitalize())
    adata.var.index = adata.var.index.map(lambda x: x.capitalize())
    adata.var.index.name = "MGI_symbol"
    adata.var.columns = ["original_id"]

    return adata


def parse(data_dict: dict[str, dict[str,Path]], assay: str) -> ad.AnnData:
    """
    Reads in raw 10x data and call files, applies light correction of calls based on frequency and total count,
    filters on MOI, and returns an adata object with standardised relation convention.
    Now also takes the assay as a parameter to fetch specific configurations.

    :param data_dict: dict of data : path to that data
    :param assay: Assay identifier, used to load specific configurations
    :return: standardised anndata object
    :rtype: anndata.Anndata
    """
    assay_config = get_assay_config(assay)
    if not assay_config:
        raise ValueError(f"No configuration found for assay: {assay}")

    # Retrieve mappings from configuration
    cell_suffix_to_condition: dict[str, str] = assay_config.get(
        "cell_suffix_to_condition", {}
    )
    cell_suffix_to_replicate: dict[str, str] = assay_config.get(
        "cell_suffix_to_replicate", {}
    )
    assay_description = assay_config.get("assay_description", "")

    guide_counts = load_call_files(data_dict["call_files"])
    cell_df = filter_on_moi(guide_counts)
    adata = load_count_data(data_dict["count_files"])
    adata = filter_cells_from_counts(adata, cell_df)
    adata = standardise_obs(adata, cell_suffix_to_condition, cell_suffix_to_replicate)
    adata = standardise_var(adata)

    # Add study info
    adata.uns["single_cell_technique"] = "10x Genomics Chromium"
    adata.uns[
        "source"
    ] = f"Lara mouse fibroblast resequencing data 2023 {assay_description}"

    return adata