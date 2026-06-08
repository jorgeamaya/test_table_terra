version 1.0

import "../Tasks/prepare_fasta_files/prepare_fasta_files.wdl" as test_prepare_fasta_files_t

workflow TestEnvironment {
	input {
		String bucket_name
		String screen_id
		String query_name
		String query_sequence
		Array[File] subject_proteome_datasets
	}

	call test_prepare_fasta_files_t.TestPrepareFastaFiles as t_001_prepare_fasta_files {
		input:
			bucket_name = bucket_name,
			screen_id = screen_id,
			query_name = query_name,
			query_sequence = query_sequence,
			subject_proteome_datasets = subject_proteome_datasets
	}

	output {
		File fasta_input_inventory = t_001_prepare_fasta_files.fasta_input_inventory
	}
}
