version 1.0

import "../Tasks/prepare_fasta_files/prepare_fasta_files.wdl" as prepare_fasta_files_t
import "../Tasks/local_msa_colabfold_search/local_msa_colabfold_search_test.wdl" as local_msa_colabfold_search_t
import "../Tasks/predict_with_colabfold/predict_with_colabfold_test.wdl" as predict_with_colabfold_t

workflow ProtBindScreenSubmitMode {
	input {
		String query_name
		String query_sequence
		File subject_native_sequences_file 
		File subject_scrambled_sequences_file
		String colabfold_db_path	
	}

	call prepare_fasta_files_t.PrepareFastaFiles as t_001_prepare_fasta_files {
		input:
			query_name = query_name,
			query_sequence = query_sequence,
			subject_native_sequences_file = subject_native_sequences_file,
			subject_scrambled_sequences_file = subject_scrambled_sequences_file
	}

	call local_msa_colabfold_search_t.LocalMSAColabfoldSearch as t_002_local_msa_colabfold_search {
		input:
			query_name = query_name,
			colabfold_input_files = t_001_prepare_fasta_files.colabfold_input_files,
			colabfold_db_path = colabfold_db_path
	}

	call predict_with_colabfold_t.PredictWithColabfold as t_003_predict_with_colabfold {
		input:
			query_name = query_name,
			a3m_files = t_002_local_msa_colabfold_search.colabfold_a3m_files
	}

	output {
		Array[File] colabfold_input_files = t_001_prepare_fasta_files.colabfold_input_files
		Array[File] colabfold_output_files = t_003_predict_with_colabfold.colabfold_output_files
	}
}
