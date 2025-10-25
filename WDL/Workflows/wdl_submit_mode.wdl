version 1.0

import "../Tasks/prepare_fasta_files/prepare_fasta_files.wdl" as prepare_fasta_files_t
#import "../Tasks/predict_with_colabfold/predict_with_colabfold.wdl" as predict_with_colabfold_t

workflow ProtBindScreenSubmitMode {
	input {
		String query_name
		String query_sequence

		File subject_native_sequences_file 
		File subject_scrambled_sequences_file
	}

	call prepare_fasta_files_t.PrepareFastaFiles as t_001_prepare_fasta_files {
		input:
			query_name = query_name,
			query_sequence = query_sequence,
			subject_native_sequences_file = subject_native_sequences_file,
			subject_scrambled_sequences_file = subject_scrambled_sequences_file
	}
	
#	scatter (fasta_file in t_001_prepare_fasta_files.colabfold_input_files) {
#		call predict_with_colabfold_t.PredictWithColabfold as t_002_predict_with_colabfold {
#			input:
#				fasta_file = fasta_file
#		}
#	}

	output {
		Array[File] colabfold_input_files = t_001_prepare_fasta_files.colabfold_input_files
		#Array[File] colabfold_output_files = flatten(t_002_predict_with_colabfold.colabfold_output_files)
	}
}
