#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Docstring___________________________________________________________________

"""**The module *submission_helper* provides helper functions used to generate FASTA files.**

Notes:
    THIS .PY WAS ADAPTED FOR WDL. THE INPUT ARGUMENTS WERE CHANGED TO FIT THE WDL TASK

Notes:
    The module submission_helper includes functions for:

    - Making FASTA files

Imports:
    - Standard libraries: pathlib
    - Third-party libraries: pandas
"""


# Imports_____________________________________________________________________

from pathlib import Path

import pandas as pd


# Screen input processing utilities___________________________________________


def make_fasta_files(
    query_protein_path: Path,
    subject_native_sequences_path: Path,
    subject_scrambled_sequences_path: Path,
    predictions_dir: Path,
) -> int:
    """Helper function to generate FASTA files.

    Args:
        query_protein_path (Path): Path to the query protein FASTA file.
        subject_native_sequences_path (Path): Path to the subject native
            sequences FASTA file.
        subject_scrambled_sequences_path (Path): Path to the subject scrambled
            sequences FASTA file.
        predictions_dir (Path): Directory where FASTA files will be saved.

    Raises:
        FileNotFoundError: If any of the input files are not found.
        ValueError: For any errors that occur during FASTA file generation.

    Side Effects:
        - Writes FASTA files in the specified predictions directory.

    Returns:
        int: The number of generated FASTA files.

    Notes:
        - The function will create FASTA files in the
        specified predictions directory.
        - The naming convention for FASTA files is:
        {query_name}_{full_subject_id}.fasta
        - Each FASTA file will contain the query sequence
        and one subject sequence.
        - The total number of FASTA files equals the number of
        subject sequences.
    """

    # Check if input files and predictions directory exist
    if not query_protein_path.exists():
        raise FileNotFoundError(
            "Query protein file not found:" f" {query_protein_path}"
        )
    if not subject_native_sequences_path.exists():
        raise FileNotFoundError(
            "Subject native sequences file not found:"
            f" {subject_native_sequences_path}"
        )
    if not subject_scrambled_sequences_path.exists():
        raise FileNotFoundError(
            "Subject scrambled sequences file not found:"
            f" {subject_scrambled_sequences_path}"
        )
    if not predictions_dir.exists():
        predictions_dir.mkdir(parents=True, exist_ok=True)

    # Load input files into DataFrames
    try:
        query_df = pd.read_csv(query_protein_path, sep="\t")
        spns_df = pd.read_csv(subject_native_sequences_path, sep="\t")
        spss_df = pd.read_csv(subject_scrambled_sequences_path, sep="\t")

    except Exception as e:
        raise ValueError(
            f"Error occurred while reading input files: {e}"
        ) from e

    # Create a sorted DataFrame for allsubject sequences
    try:
        subject_df = pd.concat([spns_df, spss_df], ignore_index=True)
    except Exception as e:
        raise ValueError(
            f"Error occurred while processing subject sequences: {e}"
        ) from e

    try:
        subject_df = subject_df.sort_values(by="full_subject_id").reset_index(
            drop=True
        )
    except Exception as e:
        raise ValueError(
            f"Error occurred while sorting subject sequences: {e}"
        ) from e

    # Create the FASTA files
    try:
        for query_row in query_df.itertuples(index=False):
            query_name = query_row.query_name
            query_seq = query_row.query_sequence

            for subject_row in subject_df.itertuples(index=False):
                full_subject_id = subject_row.full_subject_id
                subject_seq = subject_row.sequence

                fasta_filename = f"{query_name}_{full_subject_id}.fasta"
                with open(predictions_dir / fasta_filename, "w") as f:
                    f.write(f">{query_name}_{full_subject_id}\n")
                    f.write(f"{query_seq}:{subject_seq}\n")
    except Exception as e:
        raise ValueError(
            f"Error occurred while creating FASTA files: {e}"
        ) from e

    count = sum(
        1 for fasta in predictions_dir.glob("*.fasta") if fasta.is_file()
    )

    return count


