#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Docstring____________________________________________________________________

"""**The module *analysis_helper* provides helper functions used by
screen_analysis.py in the ProtBindScreen package.**

Notes:
    The module analysis_helper includes functions for:

        - Filtering subject proteome dictionaries for completed runs.
        - Parsing JSON files to compute prediction scores.
        - Running Principal Component Analysis (PCA) on prediction scores.
        - Calculating Euclidean distances between native and scrambled
          protein predictions in PCA space.
        - Generating static plots for PCA and volcano plots.

Imports:
    - Standard libraries: fnmatch, json, pathlib, typing
    - Third-party libraries: matplotlib, numpy, pandas, seaborn, scikit-learn,
      statsmodels
"""


# Imports______________________________________________________________________

import fnmatch
import json
from pathlib import Path
from typing import Tuple


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.multitest import multipletests

# Filter subject proteome utilities___________________________________________


def filter_subject_proteome_complete(
    spd_path: Path, all_fasta_files_and_the_prediction_outputs: list[Path]
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Helper function to filter the subject proteome dictionary for completed
    runs.

    Args:
        spd_path (Path): Path to the subject proteome dictionary TSV file.
        all_fasta_files_and_the_prediction_outputs (list[Path]): List of
            paths to all fasta files and the prediction outputs.

    Raises:
        ValueError: If the subject proteome dictionary cannot be loaded or if
            no completed runs are found.
        IOError: If the filtered subject proteome dictionaries TSV files
            cannot be saved.

    Returns:
        Tuple:
            - pd.DataFrame: A pandas DataFrame containing the filtered
            subject proteome dictionary for completed runs.
            - pd.DataFrame: A pandas DataFrame containing the filtered
            subject proteome dictionary for incomplete runs.
    """

    # Load the subject proteome dictionary
    try:
        spd_df = pd.read_csv(spd_path, sep="\t")
    except Exception as e:
        raise ValueError(
            f"Failed to load subject proteome dictionary: {e}"
        ) from e

    spd_completion_checks_df = spd_df.copy()

    # Add columns for completion checks
    spd_completion_checks_df = spd_completion_checks_df.assign(
        fasta_file_n_count=pd.Series(dtype="Int64"),
        done_txt_file_n_count=pd.Series(dtype="Int64"),
        fasta_file_s_count=pd.Series(dtype="Int64"),
        done_txt_file_s_count=pd.Series(dtype="Int64"),
        run_completion=pd.Series(dtype="boolean"),
    )

    # Loop through each subject ID and check for completion
    for idx, subject_id in spd_completion_checks_df["subject_id"].items():
        fasta_n_matches = [f for f in all_fasta_files_and_the_prediction_outputs if f.name.endswith(f"_{subject_id}n.fasta")]
        fasta_s_matches = [f for f in all_fasta_files_and_the_prediction_outputs if f.name.endswith(f"_{subject_id}s.fasta")]
        done_n_matches  = [f for f in all_fasta_files_and_the_prediction_outputs if f.name.endswith(f"_{subject_id}n.done.txt")]
        done_s_matches  = [f for f in all_fasta_files_and_the_prediction_outputs if f.name.endswith(f"_{subject_id}s.done.txt")]

        fasta_n_count = len(fasta_n_matches)
        fasta_s_count = len(fasta_s_matches)
        done_n_count = len(done_n_matches)
        done_s_count = len(done_s_matches)

        run_complete = (
            fasta_n_count == 1
            and done_n_count == 1
            and fasta_s_count == 1
            and done_s_count == 1
        )

        # Write to dataframe
        spd_completion_checks_df.at[idx, "fasta_file_n_count"] = fasta_n_count
        spd_completion_checks_df.at[idx, "done_txt_file_n_count"] = (
            done_n_count
        )
        spd_completion_checks_df.at[idx, "fasta_file_s_count"] = fasta_s_count
        spd_completion_checks_df.at[idx, "done_txt_file_s_count"] = (
            done_s_count
        )
        spd_completion_checks_df.at[idx, "run_completion"] = run_complete

    # Filter subject proteome dictionary based on completion checks
    if not spd_completion_checks_df["run_completion"].any():
        raise ValueError(
            "No completed runs found."
            " Please check the prediction results directory."
        )
    spd_df_complete_runs = spd_completion_checks_df[
        spd_completion_checks_df["run_completion"].eq(True)
    ].copy()
    spd_df_incomplete_runs = spd_completion_checks_df[
        spd_completion_checks_df["run_completion"].eq(False)
    ].copy()

    return spd_df_complete_runs, spd_df_incomplete_runs


# Json parser utilities_______________________________________________________


def parse_json_compute_scores(
    subject_ids_list: list,
    all_fasta_files_and_the_prediction_outputs: list[Path],
    included_matrices_dict: dict,
    query_len: int,
) -> pd.DataFrame:
    """Helper function to parse JSON files and compute all prediction scores.

    Args:
        subject_ids_list (list): List of subject IDs to process.
        all_fasta_files_and_the_prediction_outputs (list[Path]): List of paths to all FASTA files and prediction outputs.
        included_matrices_dict (dict): Dictionary of included matrices with
            amino acid ranges.
        query_len (int): Length of the query protein.

    Returns:
        pd.DataFrame: DataFrame containing the computed prediction scores.
    """

    # Set constants
    MODEL_RANKS = ["1", "2", "3", "4", "5"]
    CONFIDENCE_THRESHOLD = 1 / (
        1 + 20
    )  # corresponds to pAE of 20, i.e. 0.047619

    matrices = list(included_matrices_dict.keys())
    query_coords = []  # The values stored in query_coords are 0-based indices
    for matrix in matrices:
        matrix_data = included_matrices_dict[matrix]
        aa_ranges_i = matrix_data["aa_ranges_i"]
        for part in aa_ranges_i.split(","):
            if "-" in part:
                start, end = map(int, part.split("-"))
                query_coords.extend(range(start - 1, end))
            else:
                query_coords.append(int(part) - 1)

    # Format a seed DataFrame for PCA input
    # Initialize the seed DataFrame with the required columns
    empty_df = pd.DataFrame()
    seed_df = empty_df.assign(
        full_subject_id=pd.Series(dtype="string"),
        subject_id=pd.Series(dtype="string"),
        nstag=pd.Series(dtype="string"),
    )

    for model in MODEL_RANKS:
        subject_len_col = f"m{model}_subject_len"
        seed_df = seed_df.assign(**{subject_len_col: pd.Series(dtype="Int64")})
        plddt_col = f"m{model}_plddt"
        ptm_col = f"m{model}_ptm"
        iptm_col = f"m{model}_iptm"

        std_scores_columns = [plddt_col, ptm_col, iptm_col]
        seed_df = seed_df.assign(
            **{col: pd.Series(dtype="float64") for col in std_scores_columns}
        )

    for model in MODEL_RANKS:
        for coord in query_coords:
            coord1 = coord + 1  # Convert to 1-based index for column names
            for strip_type in ["a", "b"]:
                strip_name = f"m{model}_s{coord1}{strip_type}"
                seed_df = seed_df.assign(
                    **{strip_name: pd.Series(dtype="float64")}
                )

    # Populate the identification columns of the seed DataFrame
    identification_rows = []
    for subject_id in subject_ids_list:
        for nstag in ["n", "s"]:
            identification_rows.append(
                {
                    "full_subject_id": f"{subject_id}{nstag}",
                    "subject_id": subject_id,
                    "nstag": nstag,
                }
            )
    seed_df = pd.concat(
        [seed_df, pd.DataFrame(identification_rows)], ignore_index=True
    )

    # Process each JSON file to get the I scores per strip.
    # The confidence matrix is derived from the pAE matrix as
    # confidence = 1/(1 + pAE)
    # I scores are a sum of filtered confidence scores for each strip
    # THRESHOLD = 1/(1 + 20) ~ 0.047619, corresponds to a pAE of 20.

    data_rows = []
    full_subject_ids_list = seed_df["full_subject_id"].unique().tolist()
    for full_subject_id in full_subject_ids_list:
        skipped_subjects = {}
        data_row = {}
        data_row["full_subject_id"] = full_subject_id
        for model in MODEL_RANKS:
            # Load the JSON file for the current subject and model rank
            try:
                json_file_path = next(
                    (
                        f for f in all_fasta_files_and_the_prediction_outputs
                        if fnmatch.fnmatch(f.name, f"*_{full_subject_id}_*_rank_00{model}*.json")
                    )
                )
            except StopIteration:
                skipped_subjects.setdefault(full_subject_id, []).append(model)
                continue

            try:
                with open(json_file_path, "r") as f:
                    json_data = json.load(f)
            except (json.JSONDecodeError, OSError):
                skipped_subjects.setdefault(full_subject_id, []).append(model)
                continue

            # Load data from JSON
            try:
                pae_matrix = np.array(json_data["pae"])
            except (KeyError, ValueError):
                skipped_subjects.setdefault(full_subject_id, []).append(model)
                continue

            confidence_matrix = 1 / (1 + pae_matrix)

            total_len = int(confidence_matrix.shape[0])
            query_len = int(query_len)
            subject_len = int(total_len - query_len)
            data_row[f"m{model}_subject_len"] = subject_len

            try:
                plddt = float(round(np.mean(json_data["plddt"]), 2))
            except (KeyError, ValueError):
                skipped_subjects.setdefault(full_subject_id, []).append(model)
                continue
            data_row[f"m{model}_plddt"] = plddt

            try:
                ptm = float(round(json_data["ptm"], 2))
            except (KeyError, ValueError):
                skipped_subjects.setdefault(full_subject_id, []).append(model)
                continue
            data_row[f"m{model}_ptm"] = ptm

            try:
                iptm = float(round(json_data["iptm"], 2))
            except (KeyError, ValueError):
                skipped_subjects.setdefault(full_subject_id, []).append(model)
                continue
            data_row[f"m{model}_iptm"] = iptm

            # Keep the original list of matrices and query coords (see above)
            matrices = matrices
            query_coords = query_coords

            # Get the subject coords for the analyzed matrices
            # and compute the interaction scores
            for matrix in matrices:
                matrix_data = included_matrices_dict[matrix]
                aa_ranges_j = matrix_data["aa_ranges_j"]

                subject_coords = []
                # The values stored in subject_coords are 0-based indices
                for i in aa_ranges_j.split(","):
                    i = i.strip()
                    if "-" in i:
                        start, end = map(int, i.split("-"))
                        subject_coords.extend(range(start - 1, end))
                    elif ":" in i:
                        start = int(i[:-1])
                        subject_coords.extend(range(start - 1, total_len))
                    else:
                        subject_coords.append(int(i) - 1)

                for coord in query_coords:
                    # Extract the confidence scores for strips a and b
                    strip_a = [
                        confidence_matrix[coord, sc] for sc in subject_coords
                    ]
                    strip_b = [
                        confidence_matrix[sc, coord] for sc in subject_coords
                    ]

                    # Calculate the interaction scores for strips a and b
                    # by summing confidence scores above the threshold
                    try:
                        i_score_a = float(
                            round(
                                sum(
                                    confidence_score
                                    for confidence_score in strip_a
                                    if confidence_score >= CONFIDENCE_THRESHOLD
                                ),
                                3,
                            )
                        )
                    except (TypeError, ValueError):
                        skipped_subjects.setdefault(
                            full_subject_id, []
                        ).append(model)
                        continue
                    try:
                        i_score_b = float(
                            round(
                                sum(
                                    confidence_score
                                    for confidence_score in strip_b
                                    if confidence_score >= CONFIDENCE_THRESHOLD
                                ),
                                3,
                            )
                        )
                    except (TypeError, ValueError):
                        skipped_subjects.setdefault(
                            full_subject_id, []
                        ).append(model)
                        continue
                    coord1 = (
                        coord + 1
                    )  # Convert to 1-based index for column names
                    data_row[f"m{model}_s{coord1}a"] = i_score_a
                    data_row[f"m{model}_s{coord1}b"] = i_score_b
        if skipped_subjects.get(full_subject_id):
            continue
        data_rows.append(data_row)

    data_df = pd.DataFrame(data_rows)
    seed_df = seed_df.set_index("full_subject_id")
    data_df = data_df.set_index("full_subject_id")
    seed_df.update(data_df)
    seed_df = seed_df.reset_index()
    seed_df = seed_df.dropna(
        subset=[col for col in seed_df.columns if col != "full_subject_id"],
        how="any",
    )
    seed_df = seed_df.sort_values(by=["full_subject_id"]).reset_index(
        drop=True
    )
    scores_complete_runs_df = seed_df

    return scores_complete_runs_df


# PCA utilities_______________________________________________________________


def run_pca_on_completed_runs(
    scores_complete_runs_df: pd.DataFrame,
) -> Tuple[dict, dict]:
    """Helper function to run Principal Component Analysis (PCA) on the
    prediction scores data for completed runs.

    Args:
        scores_complete_runs_df (pd.DataFrame): DataFrame containing
            prediction scores.

    Raises:
        ValueError: If PCA computation fails due to empty DataFrame or
            missing required columns.

    Returns:
        Tuple[dict, dict]: PCA results for model rank 1 and mean of
        model ranks 1-5.

    Notes:
        - PCA is performed separately on scores derived from model rank 1
        (m1) and the mean of model ranks 1-5 (m15mean).
        - The function returns dictionaries containing PCA objects,
        explained variance, PCA results DataFrames, and standardized input
        DataFrames for both m1 and m15mean datasets.
    """

    scores_complete_runs_df = scores_complete_runs_df.copy()
    #  Dataset with scores derived from model rank 1 (m1) only
    m1_scores_df = scores_complete_runs_df[["full_subject_id"]].join(
        scores_complete_runs_df.filter(like="m1_s")
    )

    #  Dataset with scores derived from mean of model ranks 1-5 (m15mean)
    id_only_df = scores_complete_runs_df[["full_subject_id"]].copy()
    # Calculate mean scores across models 1-5 for each strip
    scores_columns = scores_complete_runs_df.columns[
        scores_complete_runs_df.columns.str.match(r"m[1-5]_s\d+[ab]")
    ]
    scores_only_df = scores_complete_runs_df[scores_columns].copy()
    mean_scores_df = (
        scores_only_df.T.groupby(lambda idx: idx.split("_", 1)[1])
        .mean()
        .T.add_prefix("m15mean_")
    )
    m15mean_scores_df = pd.concat([id_only_df, mean_scores_df], axis=1)

    # Run PCA for each dataset and return results in a dictionary
    def run_pca(data: pd.DataFrame, label: str) -> dict:
        """Run PCA on the provided DataFrame and return results.

        Args:
            data (pd.DataFrame): Input DataFrame with 'full_subject_id' and
                numerical columns.
            label (str): Label for the dataset (e.g., 'm1' or 'm15mean').

        Raises:
            ValueError: If no numerical columns with variance > 0 are found or
                if data contains NaN or infinite values.

        Returns:
            dict:
                - pca: Fitted PCA object.
                - total_var: Total explained variance by selected components.
                - var_details: Explained variance ratio for each component.
                - pca_results_df: DataFrame with PCA results including
                    'full_subject_id'.
                - pca_input_xscaled_df: DataFrame with standardized input data.
                - pca_input_xreconstructed_df: DataFrame with inverse
                    transformed data from PCA components.

        Notes:

            Principal Component Analysis (PCA) is a statistical technique used
            to reduce the dimensionality of a dataset while preserving as much
            variance as possible. It transforms the original features into a
            new set of uncorrelated features (principal components) ordered by
            the amount of variance they capture.

            - In this function, PCA is applied to the numerical columns of the
            input DataFrame.
            - The data is first standardized to have a mean of 0 and a standard
            deviation of 1.
            - PCA is then performed to retain components that explain 95% of
            the variance, using a full singular Value Decomposition (SVD)
            solver.
            - The data is inverse transformed to reconstruct it in the
            original (scaled) feature space.
            - The function returns a dictionary containing the PCA object,
            total explained variance, variance details, PCA results DataFrame,
            standardized input DataFrame, and inverse transformed DataFrame.


        """
        # Get only numerical columns for PCA
        if data["full_subject_id"].duplicated().any():
            raise ValueError(
                "Duplicate full_subject_id values found in input data."
            )

        data = data.set_index("full_subject_id")
        ids_index = data.index

        variables_df = data.select_dtypes(include=["float64"])

        # Remove columns with zero variance
        variables_df = variables_df.loc[:, variables_df.std() > 0]

        # Fail if no numerical columns with variance > 0
        nonzero_mask = variables_df.std(ddof=0) > 0
        variables_df = variables_df.loc[:, nonzero_mask]

        if variables_df.shape[1] == 0:
            raise ValueError(
                "No numerical columns with variance > 0 found for PCA."
            )

        # Fail if there are any missing values
        if not np.isfinite(variables_df.to_numpy()).all():
            raise ValueError("Data for PCA contains NaN or infinite values.")

        # Standardize the data before PCA
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(variables_df)
        assert X_scaled.shape[0] == len(ids_index)

        # Run PCA to retain 95% of variance
        pca = PCA(n_components=0.95, svd_solver="full")
        principal_components = pca.fit_transform(X_scaled)
        assert principal_components.shape[0] == len(ids_index)

        # Inverse transform to reconstruct data
        # in the original (scaled) feature space
        X_reconstructed = pca.inverse_transform(principal_components)
        assert X_reconstructed.shape == X_scaled.shape

        # Create a DataFrame for PCA results
        pc_cols = [
            f"{label}_PC{i+1}" for i in range(principal_components.shape[1])
        ]
        pc_data_only_df = pd.DataFrame(
            principal_components, index=ids_index, columns=pc_cols
        )
        pca_results_df = pc_data_only_df.reset_index().rename(
            columns={"index": "full_subject_id"}
        )

        # Calculate explained variance details
        total_var = float(pca.explained_variance_ratio_.sum())
        var_details = {
            f"{label}_PC{i+1}": float(var)
            for i, var in enumerate(pca.explained_variance_ratio_)
        }

        # Capture scaled data and inverse transformed data for documentation
        nonzerovar_variables_names = variables_df.columns.tolist()
        scaled_data_only_df = pd.DataFrame(
            X_scaled, index=ids_index, columns=nonzerovar_variables_names
        )
        pca_input_xscaled_df = scaled_data_only_df.reset_index().rename(
            columns={"index": "full_subject_id"}
        )

        reconstructed_data_only_df = pd.DataFrame(
            X_reconstructed,
            index=ids_index,
            columns=nonzerovar_variables_names,
        )
        pca_input_xreconstructed_df = (
            reconstructed_data_only_df.reset_index().rename(
                columns={"index": "full_subject_id"}
            )
        )

        assert (
            pca_results_df["full_subject_id"].values
            == pca_input_xscaled_df["full_subject_id"].values
        ).all()
        assert (
            pca_results_df["full_subject_id"].values
            == pca_input_xreconstructed_df["full_subject_id"].values
        ).all()

        pca_results_dict = {
            "model": pca,
            "total_var": total_var,
            "var_details": var_details,
            "results_df": pca_results_df,
            "input_xscaled_df": pca_input_xscaled_df,
            "input_xreconstructed_df": pca_input_xreconstructed_df,
        }

        return pca_results_dict

    try:
        m1_pca_results_dict = run_pca(m1_scores_df, label="m1")
    except Exception as e:
        raise ValueError(f"PCA computation failed for m1 dataset: {e}") from e

    try:
        m15mean_pca_results_dict = run_pca(m15mean_scores_df, label="m15mean")
    except Exception as e:
        raise ValueError(
            f"PCA computation failed for m15mean dataset: {e}"
        ) from e

    return m1_pca_results_dict, m15mean_pca_results_dict


# N-S distance utilities______________________________________________________


def run_euclidean_distance_analysis_on_pca_results(
    pca_results_df: pd.DataFrame,
) -> Tuple[pd.DataFrame, dict]:
    """Helper function to run Euclidean distance analysis on PCA results.

    Args:
        pca_results_df (pd.DataFrame): DataFrame containing PCA results.

    Returns:
        Tuple:
            - pd.DataFrame: DataFrames containing NS distances
            - dict: Dictionary containing hits summary.

    Notes:
        - n = native
        - s = scrambled
        - PC1, PC2, ... PCx = principal components from PCA analysis,
        where x is the number of components that explain 95% of the variance.
        - The input Dataframe is expected to have columns:
            - full_subject_id
            - PC1
            - PC2
            - ...
            - PCx
        - full_subject_id is formatted as {subject_id}{nstag}, where nstag is
        either 'n' or 's'.
        - The function calculates Euclidean distance between the PCA
        coordinates of native (n) and scrambled (s) subjects for each
        subject ID in the multi-dimensional PCA space defined by the
        principal components (PC1 o PCx).

        - The function also computes z-scores, min-max scaled distances,
        p-values, FDRs, negative log10 p-values, and negative log10 FDRs for the NS distances.

        - A null distribution for p-values is generated by randomly
        shuffling the coordinates and calculating distances for the
        shuffled coordinates.

        - The hits summary includes total hits, significant hits,
        and their respective p-values.
    """

    # Reshape the PCA results DataFrame for distance stats calculation
    pc_df = pca_results_df.copy()

    # Split full_subject_id into subject_id and nstag columns
    # Drop full_subject_id column
    ids = pc_df["full_subject_id"].astype(str)
    pc_df.insert(1, "subject_id", ids.str[:-1])
    pc_df.insert(2, "nstag", ids.str[-1])
    pc_df = pc_df.drop(columns=["full_subject_id"])

    # Reshape from wide to long format
    # Then back to wide with multi-index
    # to have columns PC1_n, PC1_s, PC2_n, PC2_s, ...
    # Converty pd.DataFrame to Numpy array for distance calculations

    pc_cols = [
        c
        for c in pc_df.columns
        if c.startswith("m1_PC") or c.startswith("m15mean_PC")
    ]
    if not pc_cols:
        raise ValueError("No PCA columns found in the input DataFrame.")

    pc_cols = sorted(
        pc_cols,
        key=lambda s: int(s.split("_", 1)[1][2:]),
    )

    pc_long_df = pc_df.melt(
        id_vars=["subject_id", "nstag"],
        value_vars=pc_cols,
        var_name="pc_axis",
        value_name="pc_coord",
    )
    pc_wide_midx_df = pc_long_df.pivot(
        index=["subject_id"], columns=["pc_axis", "nstag"], values="pc_coord"
    )
    pc_wide_midx_df = pc_wide_midx_df.reindex(
        columns=pd.MultiIndex.from_product([pc_cols, ["n", "s"]])
    )

    pc_wide_midx_df.columns = [
        f"{col[0]}_{col[1]}" for col in pc_wide_midx_df.columns
    ]

    pc_wide_midx_df = pc_wide_midx_df.reset_index()

    if pc_wide_midx_df.isna().any().any():
        raise ValueError("Missing PC_nstag values after reshape.")

    # Calculate Euclidean distances between n and s in all dimensions
    # for each subject
    eud_stats_df = pc_wide_midx_df.copy()
    hits_summary = {}
    n_cols = [f"{col}_n" for col in pc_cols]
    s_cols = [f"{col}_s" for col in pc_cols]

    if len(n_cols) != len(s_cols):
        raise ValueError("Mismatch in number of n and s PC columns.")

    diff = eud_stats_df[n_cols].to_numpy(dtype="float64") - eud_stats_df[
        s_cols
    ].to_numpy(dtype="float64")

    eud_stats_df["ns_euclidean_distance"] = np.sqrt(np.sum(diff**2, axis=1))

    # Block permute the s coordinates and calculate null distances
    # to generate a null distribution for p-values
    n_only_ndarray = eud_stats_df[n_cols].to_numpy(dtype="float64")
    s_only_ndarray = eud_stats_df[s_cols].to_numpy(dtype="float64")

    null_distances = []
    rng = np.random.default_rng(seed=42)
    for _ in range(10000):
        # Randomly shuffle the coordinates
        shuffled_s_only_ndarray = rng.permutation(s_only_ndarray)

        # Calculate the distance for the shuffled coordinates
        null_diff = n_only_ndarray - shuffled_s_only_ndarray
        null_distance = np.sqrt(np.sum(null_diff**2, axis=1))
        null_distances.append(null_distance)

    null_distances = np.concatenate(null_distances, axis=None)
    len_null = len(null_distances)

    observed_distances = eud_stats_df["ns_euclidean_distance"].to_numpy(
        dtype="float64"
    )

    # Calculate p-values, FDRs, and negative log10 FDRs
    p_values = np.array(
        [
            (np.sum(null_distances >= d) + 1) / (len_null + 1)
            for d in observed_distances
        ]
    )

    eud_stats_df["ns_euclidean_distance_pvalue"] = p_values

    _, fdrs, _, _ = multipletests(p_values, method="fdr_bh")
    eud_stats_df["ns_euclidean_distance_fdr"] = fdrs

    neglog10p = -np.log10(eud_stats_df["ns_euclidean_distance_pvalue"])
    eud_stats_df["ns_euclidean_distance_neglog10pvalue"] = neglog10p

    neglog10fdr = -np.log10(eud_stats_df["ns_euclidean_distance_fdr"])
    eud_stats_df["ns_euclidean_distance_neglog10fdr"] = neglog10fdr

    # Prepare hits summary
    total_hits = int(eud_stats_df.shape[0])
    significant_hits = int(
        (eud_stats_df["ns_euclidean_distance_fdr"] < 0.05).sum()
    )
    significant_rate = (
        float((significant_hits / total_hits) * 100)
        if total_hits > 0
        else float("nan")
    )
    hits_summary = {
        "total_hits": total_hits,
        "significant_hits": significant_hits,
        "significant_rate": significant_rate,
    }
    return eud_stats_df, hits_summary


# Graph generation utilities__________________________________________________


def make_static_plots(
    viewer_pca_source_df: pd.DataFrame, viewer_volcano_source_df: pd.DataFrame
) -> plt.Figure:
    """Helper function to create static plots for PCA and volcano plots.

    Args:
        viewer_pca_source_df (pd.DataFrame): DataFrame containing PCA results
            for visualization.
        viewer_volcano_source_df (pd.DataFrame): DataFrame containing volcano
            plot data for visualization.

    Returns:
        plt.Figure: Matplotlib Figure object containing the plots.

    Notes:
        - The function creates a 2x2 grid of scatter plots using seaborn and
        matplotlib.
        - The top row contains volcano plots for model rank 1 and the mean of
        the top 5 models.
        - The bottom row contains PCA plots for model rank 1 and the mean of
        the top 5 models.
        - Each plot includes appropriate titles, axis labels, and styling.
    """
    # Create a 2x2 grid of scatter plots
    sns.set_theme(style="whitegrid", font_scale=1.2)
    fig, ax = plt.subplots(2, 2, figsize=(12, 12))

    # Volcano plots
    sns.scatterplot(
        data=viewer_volcano_source_df,
        x="m1_ns_euclidean_distance",
        y="m1_ns_euclidean_distance_neglog10pvalue",
        ax=ax[0, 0],
        s=100,
        legend=False,
        color="blue",
        alpha=0.6,
    )
    ax[0, 0].set_title("Volcano Plot - Model Rank 1")
    ax[0, 0].set_xlabel("Native to Scrambled Distance")
    ax[0, 0].set_ylabel("-log10(p-value)")

    sns.scatterplot(
        data=viewer_volcano_source_df,
        x="m15mean_ns_euclidean_distance",
        y="m15mean_ns_euclidean_distance_neglog10pvalue",
        ax=ax[0, 1],
        s=100,
        legend=False,
        color="orange",
        alpha=0.6,
    )
    ax[0, 1].set_title("Volcano Plot - Mean of Top 5 Models")
    ax[0, 1].set_xlabel("Native to Scrambled Distance")
    ax[0, 1].set_ylabel("-log10(p-value)")

    # PCA plots
    sns.scatterplot(
        data=viewer_pca_source_df,
        x="m1_PC1",
        y="m1_PC2",
        ax=ax[1, 0],
        s=100,
        legend=False,
        color="blue",
        alpha=0.6,
    )
    ax[1, 0].set_title("PCA Plot - Model Rank 1")
    ax[1, 0].set_xlabel("PC1")
    ax[1, 0].set_ylabel("PC2")

    sns.scatterplot(
        data=viewer_pca_source_df,
        x="m15mean_PC1",
        y="m15mean_PC2",
        ax=ax[1, 1],
        s=100,
        legend=False,
        color="orange",
        alpha=0.6,
    )
    ax[1, 1].set_title("PCA Plot - Mean of Top 5 Models")
    ax[1, 1].set_xlabel("PC1")
    ax[1, 1].set_ylabel("PC2")

    fig.tight_layout()

    return fig
