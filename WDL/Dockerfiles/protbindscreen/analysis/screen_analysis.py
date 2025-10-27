#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Docstring____________________________________________________________________

"""**The script *screen_analysis.py* runs the analysis mode of the
ProtBindScreen tool.**

Notes:
    - THIS .PY WAS ADAPTED FOR WDL. THE INPUT ARGUMENTS WERE CHANGED TO FIT THE WDL TASK

Imports:
        - Standard libraries: argparse, datetime, joblib, json, logging,
        pathlib, sys, time
        - Third-party libraries: matplotlib, pandas
        - ProtBindScreen package:

            - **analysis_helper**, from protbindscreen.analysis
            - **housekeeping_helper**, from protbindscreen.general

"""

# Imports______________________________________________________________________

import argparse
from datetime import datetime
from joblib import dump
import json
import logging
from pathlib import Path
import sys
import time
import shutil
import matplotlib.pyplot as plt
import pandas as pd

from protbindscreen.general import housekeeping_helper as hh
from protbindscreen.analysis import analysis_helper as ah


# Main logic__________________________________________________________________
def screen_analysis(
        screen_dir: Path, 
        analysis_dir: Path, 
        analysis_name: str, 
        query_len: int, 
        all_fasta_files_and_the_prediction_outputs: list[Path], 
        subject_proteome_dictionary: Path, 
        analysis_matrices: dict) -> None:
    """Wrapper function for running the ProtBindScreen tool in analysis mode.

    Args:
        analysis_dir (Path): Path to the analysis directory.
        analysis_name (str): Name of the analysis.
        query_len (int): Length of the query protein.
        all_fasta_files_and_the_prediction_outputs (list[Path]):
            List of paths to all fasta files and their prediction outputs.
        subject_proteome_dictionary (Path):
            Path to the subject proteome dictionary file.
        analysis_matrices (dict):
            Dictionary specifying the analysis matrices.

    Side Effects:
        1. **Sets up logging**:
            - logs are directed to the parent directory of the input file
        2. **Writes files and directories**:
            - writes included matrices file
            - writes filtered subject proteome dictionary file
            - writes complete data DataFrame file
            - writes PCA results files
            - writes NS distance results files
            - writes final results file
            - writes graphs files (SVG and PDF)
        3. **sys.exit(1)**: Exits the program with an error message.
    """

    start_time = time.time()

    SCRIPT_TAG = Path(__file__).stem
    LOGGING_LEVEL = logging.INFO

    timestamp = datetime.now().strftime("%Y%m%d_%H%M")
    out_log_path = analysis_dir / f"{SCRIPT_TAG}_{timestamp}.out"
    err_log_path = analysis_dir / f"{SCRIPT_TAG}_{timestamp}.err"

    hh.setup_central_logging(out_log_path, err_log_path, level=LOGGING_LEVEL)

    logging.info("Welcome to ProtBindScreen!")
    logging.info("Running a new analysis...")
    logging.info("=" * 79)

    logging.info("Logging")
    logging.info("-" * 79)
    logging.info("Logging was initialized.")
    logging.info(f"stdout log → {out_log_path}")
    logging.info(f"stderr log → {err_log_path}")
    logging.info(
        f"The logging level was set to: {logging.getLevelName(LOGGING_LEVEL)}"
    )
    logging.info("=" * 79)

    logging.info("CLI Information")
    logging.info("-" * 79)
    logging.info(f"Running: {Path(__file__)}")

    logging.info("=" * 79)

    logging.info("Analysis Input")
    logging.info("-" * 79)

    logging.info(f"Screen directory: {screen_dir}")

    logging.info(f"Analysis directory: {analysis_dir}")

    logging.info(f"Analysis name: {analysis_name}")

    logging.info(f"Query length: {query_len}")

    analysis_matrices_dict = analysis_matrices
    included_matrices_dict = {}
    required_matrix = analysis_matrices_dict["matrix_required"]
    included_matrices_dict["matrix_required"] = required_matrix

    optional_matrices_list = [
        k for k in analysis_matrices_dict if k.startswith("matrix_optional_")
    ]
    for matrix in optional_matrices_list:
        if analysis_matrices_dict[matrix]["include"] == "true":
            included_matrices_dict[matrix] = analysis_matrices_dict[matrix]
    logging.info(
        f"Included matrices for analysis:\n"
        f"{json.dumps(included_matrices_dict, indent=4)}"
    )

    included_matrices_path = analysis_dir / "included_matrices.json"
    try:
        with open(included_matrices_path, "w") as f:
            json.dump(included_matrices_dict, f, indent=4)
    except Exception as e:
        logging.error(
            f"Failed to write included matrices file: {e}", exc_info=True
        )
        sys.exit(1)

    logging.info(f"Included matrices saved to {included_matrices_path}")
    logging.info("=" * 79)

    # Filter the subject proteome dictionary for completed runs
    logging.info("Filtering subject proteome dictionary for completed runs...")
    logging.info("-" * 79)

    ### 6.22.25 Copy subject proteome dictionary to screen directory in the Results dir of wdl task
    shutil.copy(subject_proteome_dictionary, screen_dir)

    spd_path = screen_dir / "subject_proteome_dictionary.tsv"
    if not spd_path.exists():
        logging.error(
            f"Subject proteome dictionary file '{spd_path}' does not exist."
        )
        sys.exit(1)
    logging.info(f"Subject proteome dictionary path: {spd_path}")


    try:
        spd_df_complete_runs, spd_df_incomplete_runs = (
            ah.filter_subject_proteome_complete(spd_path, all_fasta_files_and_the_prediction_outputs)
        )
    except Exception as e:
        logging.error(f"Data preparation failed: {e}", exc_info=True)
        sys.exit(1)

    # Save the filtered subject proteome dictionary
    spd_complete_runs_path = (
        analysis_dir / "subject_proteome_dictionary_complete_runs.tsv"
    )
    spd_incomplete_runs_path = (
        analysis_dir / "subject_proteome_dictionary_incomplete_runs.tsv"
    )

    try:    
        spd_df_complete_runs.to_csv(
            spd_complete_runs_path,
            sep="\t",
            index=False,
        )
        spd_df_incomplete_runs.to_csv(
            spd_incomplete_runs_path,
            sep="\t",
            index=False,
        )
    except Exception as e:
        logging.error(
            f"Failed to save filtered subject proteome dictionaries: {e}",
            exc_info=True,
        )
        sys.exit(1)

    # Rename 'native' to 'n' for merging because:
    spd_df_complete_runs["nstag"] = spd_df_complete_runs["nstag"].replace(
        {"native": "n"}
    )

    # Add filtered subject proteome dictionary to final results dataset
    final_results_df = spd_df_complete_runs.copy()

    # Add filtered subject proteome dictionary to volcano source dataset
    viewer_volcano_source_df = spd_df_complete_runs.copy()
    viewer_volcano_source_df["full_subject_id"] = (
        viewer_volcano_source_df["subject_id"]
        + viewer_volcano_source_df["nstag"]
    )
    col = viewer_volcano_source_df.pop("full_subject_id")
    viewer_volcano_source_df.insert(0, col.name, col)

    logging.info(
        "Filtered subject proteome dictionary with completed runs"
        f" saved to {spd_complete_runs_path}."
        " This file will be used for the analysis."
    )
    logging.info(
        "Filtered subject proteome dictionary with incomplete runs"
        f" saved to {spd_incomplete_runs_path}."
        " This file is for documentation only."
    )

    logging.info("Filtering completed.")
    logging.info("=" * 79)

    # Parse JSON files and compute scores for completed runs
    logging.info(
        "Parsing JSON files and computing scores for completed runs..."
    )
    logging.info("-" * 79)

    subject_ids_list = spd_df_complete_runs["subject_id"].unique().tolist()
    logging.info(
        "Number of subject proteins with completed runs:"
        f" {len(subject_ids_list)}"
    )

    try:
        scores_complete_runs_df = ah.parse_json_compute_scores(
            subject_ids_list, all_fasta_files_and_the_prediction_outputs, included_matrices_dict, query_len
        )
    except Exception as e:
        logging.error(
            f"JSON parsing and score computation failed: {e}", exc_info=True
        )
        sys.exit(1)

    scores_complete_runs_path = analysis_dir / "scores_complete_runs.tsv"
    try:
        scores_complete_runs_df.to_csv(
            scores_complete_runs_path, sep="\t", index=False, encoding="utf-8"
        )
    except Exception as e:
        logging.error(f"Failed to save complete data file: {e}", exc_info=True)
        sys.exit(1)

    # Add standard scores for native subjects to final results dataset

    std_scores_columns = scores_complete_runs_df.filter(
        regex=r"m\d+_(subject_len|plddt|ptm|iptm)$"
    ).columns.tolist()
    filtered_native_standard_scores_df = scores_complete_runs_df.loc[
        scores_complete_runs_df["nstag"] == "n",
        ["full_subject_id", "subject_id", "nstag"] + std_scores_columns,
    ]

    final_results_df = final_results_df.merge(
        filtered_native_standard_scores_df,
        on=["subject_id", "nstag"],
        how="inner",
    )

    # Reorder columns so that full_subject_id is first
    col = final_results_df.pop("full_subject_id")
    final_results_df.insert(0, col.name, col)

    logging.info(f"Complete data saved to {scores_complete_runs_path}")
    logging.info("This file will be used for further analysis.")
    logging.info("JSON parsing and score computation completed.")
    logging.info("=" * 79)

    # Run PCA on the data
    logging.info(
        "Running Principal Component Analysis (PCA)"
        " on the prediction scores data for completed runs..."
    )
    logging.info("-" * 79)

    try:
        m1_pca_results_dict, m15mean_pca_results_dict = (
            ah.run_pca_on_completed_runs(scores_complete_runs_df)
        )
    except Exception as e:
        logging.error(f"PCA computation failed: {e}", exc_info=True)
        sys.exit(1)

    # Save PCA results
    pca_metadata = {}
    pca_results_dicts = {
        "m1_pca": m1_pca_results_dict,
        "m15mean_pca": m15mean_pca_results_dict,
    }
    for name, entry in pca_results_dicts.items():
        variance_summary = {}
        for key, value in entry.items():
            # Save PCA object to a joblib file and update metadata
            if key == "model":
                path = analysis_dir / f"{name}_model.joblib"
                dump(value, path)
                pca_metadata[f"{name}_model_path"] = str(path)

            # Save the dataframes to TSV files and update metadata
            dfs = {
                "results": "results_df",
                "input_xscaled": "input_xscaled_df",
                "input_xreconstructed": "input_xreconstructed_df",
            }
            for df_name, df in dfs.items():
                if key == df:
                    path = analysis_dir / f"{name}_{df_name}.tsv"
                    df = entry[df]
                    df.to_csv(path, sep="\t", index=False, encoding="utf-8")
                    pca_metadata[f"{name}_{df_name}_path"] = str(path)

            # Save variance details to a JSON file and update metadata
            if key == "total_var":
                variance_summary["total_var"] = float(value)
            if key == "var_details":
                variance_summary["var_details"] = value
        path = analysis_dir / f"{name}_var_details.json"
        with open(path, "w") as f:
            json.dump(variance_summary, f, indent=4)
        pca_metadata[f"{name}_var_details"] = str(path)

    m1_pca_results_df = m1_pca_results_dict["results_df"]
    m15mean_pca_results_df = m15mean_pca_results_dict["results_df"]

    # Prepare PCA graph input for viewer
    pca_results_df = pd.merge(
        m1_pca_results_df,
        m15mean_pca_results_df,
        on="full_subject_id",
    )
    id_df = pca_results_df[["full_subject_id"]].copy()
    full_subject_id = id_df["full_subject_id"].astype(str)
    id_df["subject_id"] = full_subject_id.str[:-1]
    id_df["nstag"] = full_subject_id.str[-1]
    id_df = id_df.merge(
        spd_df_complete_runs[["subject_id", "identifier", "gene_name"]],
        on=["subject_id"],
        how="inner",
    )
    viewer_pca_source_df = id_df.merge(
        pca_results_df,
        on=["full_subject_id"],
        how="inner",
    )

    logging.info(
        f"PCA model for m1 saved to {pca_metadata['m1_pca_model_path']}"
    )
    logging.info(
        "PCA model for m15mean saved to"
        f" {pca_metadata['m15mean_pca_model_path']}"
    )
    logging.info(
        "PCA variance details for m1 saved to"
        f" {pca_metadata['m1_pca_var_details']}"
    )
    logging.info(
        "PCA variance details for m15mean saved to"
        f" {pca_metadata['m15mean_pca_var_details']}"
    )
    logging.info(
        "PCA results for m1 saved to" f" {pca_metadata['m1_pca_results_path']}"
    )
    logging.info(
        "PCA results for m15mean saved to"
        f" {pca_metadata['m15mean_pca_results_path']}"
    )
    logging.info(
        "PCA input xscaled for m1 saved to"
        f" {pca_metadata['m1_pca_input_xscaled_path']}"
    )
    logging.info(
        "PCA input xscaled for m15mean saved to"
        f" {pca_metadata['m15mean_pca_input_xscaled_path']}"
    )
    logging.info(
        "PCA input xreconstructed for m1 saved to"
        f" {pca_metadata['m1_pca_input_xreconstructed_path']}"
    )
    logging.info(
        "PCA input xreconstructed for m15mean saved to"
        f" {pca_metadata['m15mean_pca_input_xreconstructed_path']}"
    )

    logging.info("PCA completed.")
    logging.info("=" * 79)

    # Compute N-S distances
    logging.info("Computing N-S euclidean distances statistics...")
    logging.info("-" * 79)

    m1_ns_eud_stats_df, m1_ns_eud_hits_summary = (
        ah.run_euclidean_distance_analysis_on_pca_results(m1_pca_results_df)
    )
    m15mean_ns_eud_stats_df, m15mean_ns_eud_hits_summary = (
        ah.run_euclidean_distance_analysis_on_pca_results(
            m15mean_pca_results_df
        )
    )

    # Save N-S distance results
    m1_ns_eud_stats_df_path = analysis_dir / "m1_ns_eud_stats.tsv"
    with open(m1_ns_eud_stats_df_path, "w") as f:
        m1_ns_eud_stats_df.to_csv(f, sep="\t", index=False, encoding="utf-8")

    m15mean_ns_eud_stats_df_path = analysis_dir / "m15mean_ns_eud_stats.tsv"
    with open(m15mean_ns_eud_stats_df_path, "w") as f:
        m15mean_ns_eud_stats_df.to_csv(
            f, sep="\t", index=False, encoding="utf-8"
        )

    m1_ns_eud_hits_summary_path = analysis_dir / "m1_ns_eud_hits_summary.json"
    with open(m1_ns_eud_hits_summary_path, "w") as f:
        json.dump(m1_ns_eud_hits_summary, f, indent=4)

    m15mean_ns_eud_hits_summary_path = (
        analysis_dir / "m15mean_ns_eud_hits_summary.json"
    )
    with open(m15mean_ns_eud_hits_summary_path, "w") as f:
        json.dump(m15mean_ns_eud_hits_summary, f, indent=4)

    # Add N-S distance results to final results dataset
    m1_ns_eud_stats_df = m1_ns_eud_stats_df.rename(
        columns={
            c: f"m1_{c}"
            for c in m1_ns_eud_stats_df.columns
            if c not in ["subject_id", "nstag"]
        }
    )
    final_results_df = final_results_df.merge(
        m1_ns_eud_stats_df,
        on=["subject_id"],
        how="inner",
    )
    m15mean_ns_eud_stats_df = m15mean_ns_eud_stats_df.rename(
        columns={
            c: f"m15mean_{c}"
            for c in m15mean_ns_eud_stats_df.columns
            if c not in ["subject_id", "nstag"]
        }
    )
    final_results_df = final_results_df.merge(
        m15mean_ns_eud_stats_df,
        on=["subject_id"],
        how="inner",
    )

    # Add N-S distance results to volcano source dataset
    viewer_volcano_source_df = viewer_volcano_source_df.merge(
        m1_ns_eud_stats_df[
            [
                "subject_id",
                "m1_ns_euclidean_distance",
                "m1_ns_euclidean_distance_pvalue",
                "m1_ns_euclidean_distance_fdr",
                "m1_ns_euclidean_distance_neglog10pvalue",
                "m1_ns_euclidean_distance_neglog10fdr",
            ]
        ],
        on=["subject_id"],
        how="inner",
    )
    viewer_volcano_source_df = viewer_volcano_source_df.merge(
        m15mean_ns_eud_stats_df[
            [
                "subject_id",
                "m15mean_ns_euclidean_distance",
                "m15mean_ns_euclidean_distance_pvalue",
                "m15mean_ns_euclidean_distance_fdr",
                "m15mean_ns_euclidean_distance_neglog10pvalue",
                "m15mean_ns_euclidean_distance_neglog10fdr",
            ]
        ],
        on=["subject_id"],
        how="inner",
    )
    logging.info(f"N-S distances for m1 saved to {m1_ns_eud_stats_df_path}.")
    logging.info(
        f"N-S distances for m15mean saved to {m15mean_ns_eud_stats_df_path}."
    )
    logging.info(
        f"N-S distance summary for m1 saved to {m1_ns_eud_hits_summary_path}."
    )
    logging.info(
        "N-S distance summary for m15mean saved to"
        f" {m15mean_ns_eud_hits_summary_path}."
    )
    logging.info("NS distances computed and saved.")
    logging.info("=" * 79)

    # Make final results
    logging.info(
        "Generating final results dataset and viewer source datasets..."
    )
    logging.info("-" * 79)
    screen_name = screen_dir.name
    datetime_str = datetime.now().strftime("%Y%m%d_%H%M%S")

    # Save final results dataset
    screen_analysis_final_results_path = (
        analysis_dir
        / f"final_results_{screen_name}_{analysis_name}_{datetime_str}.tsv"
    )
    with open(screen_analysis_final_results_path, "w") as f:
        final_results_df.to_csv(f, sep="\t", index=False, encoding="utf-8")

    logging.info(
        f"Final results saved to {screen_analysis_final_results_path}."
    )

    # Save PCA graph input for viewer
    viewer_pca_source_path = analysis_dir / "viewer_pca_source.tsv"
    with open(viewer_pca_source_path, "w") as f:
        viewer_pca_source_df.to_csv(f, sep="\t", index=False, encoding="utf-8")

    logging.info(
        f"PCA graph input for viewer saved to {viewer_pca_source_path}."
    )

    viewer_volcano_source_path = analysis_dir / "viewer_volcano_source.tsv"
    with open(viewer_volcano_source_path, "w") as f:
        viewer_volcano_source_df.to_csv(
            f, sep="\t", index=False, encoding="utf-8"
        )

    logging.info(
        "Volcano graph input for viewer saved to"
        f" {viewer_volcano_source_path}."
    )
    logging.info("Final results generated.")
    logging.info("=" * 79)

    # Prepare static plots
    logging.info("Preparing static plots...")
    logging.info("-" * 79)

    try:
        fig = ah.make_static_plots(
            viewer_pca_source_df, viewer_volcano_source_df
        )
    except Exception as e:
        logging.error(f"Graph preparation failed: {e}", exc_info=True)
        sys.exit(1)

    graph_svg_path = (
        analysis_dir
        / f"static_plots_{screen_name}_{analysis_name}_{datetime_str}.svg"
    )
    graph_pdf_path = (
        analysis_dir
        / f"static_plots_{screen_name}_{analysis_name}_{datetime_str}.pdf"
    )

    fig.savefig(graph_svg_path, format="svg", bbox_inches="tight", dpi=300)
    fig.savefig(graph_pdf_path, format="pdf", bbox_inches="tight", dpi=300)
    plt.close(fig)

    logging.info(f"Graph SVG saved to {graph_svg_path}.")
    logging.info(f"Graph PDF saved to {graph_pdf_path}.")
    logging.info("Static graphs and interactive visualization data prepared.")
    logging.info("=" * 79)

    done_path = analysis_dir / f"{analysis_name}_done.txt"
    with open(done_path, "w") as f:
        f.write("ProtBindScreen analysis completed successfully.\n")
    logging.info(f"Done file written to: {done_path}")

    elapsed = time.time() - start_time
    days, rem = divmod(elapsed, 86400)
    hours, rem = divmod(rem, 3600)
    minutes, seconds = divmod(rem, 60)

    logging.info(
        "Screen analysis completed successfully in:"
        f" {days:.0f} days, {hours:.0f} hours,"
        f" {minutes:.0f} minutes, and {seconds:.0f} seconds."
    )
    logging.info("You can now proceed to the visualization step.")
    logging.info("Follow the instructions in README file.")
    logging.info("Thank you for using ProtBindScreen!")
    logging.info("Goodbye!")
    logging.info("=" * 79)


# CLI entry point_____________________________________________________________
if __name__ == "__main__":
    # Parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--screen_dir",
        type=Path,
        required=True,
        help="Path to the screen directory.",
    )
    parser.add_argument(
        "--analysis_dir",
        type=Path,
        required=True,
        help="Path to the analysis directory.",
    )
    parser.add_argument(
        "--analysis_name",
        type=str,
        required=True,
        help="Name of the analysis.",
    )
    parser.add_argument(
        "--query_len",
        type=int,
        required=True,
        help="Length of the query protein.",
    )
    parser.add_argument(
        "--all_fasta_files_and_the_prediction_outputs",
        type=list[Path],
        required=True,
        help="List of paths to all fasta files and their prediction outputs.",
    )
    parser.add_argument(
        "--subject_proteome_dictionary",
        type=Path,
        required=True,
        help="Path to the subject proteome dictionary file.",
    )
    parser.add_argument(
        "--analysis_matrices",
        type=dict,
        required=True,
        help="Dictionary specifying the analysis matrices.",
    )
    args = parser.parse_args()

    # Run the main function with parsed arguments
    screen_analysis(
        screen_dir=args.screen_dir,
        analysis_dir=args.analysis_dir,
        analysis_name=args.analysis_name,
        query_len=args.query_len,
        all_fasta_files_and_the_prediction_outputs=args.all_fasta_files_and_the_prediction_outputs,
        subject_proteome_dictionary=args.subject_proteome_dictionary,
        analysis_matrices=args.analysis_matrices,
    )
