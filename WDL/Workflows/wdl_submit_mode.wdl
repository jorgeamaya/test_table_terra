version 1.0

import "../Tasks/prepare_fasta_groups/prepare_fasta_groups.wdl" as prepare_fasta_groups_t
#import "../tasks/mock_task/mock_task.wdl" as mock_task_t
import "../Tasks/predict_with_colabfold/predict_with_colabfold_locally.wdl" as predict_with_colabfold_locally_t
#import "../tasks/predict_with_colabfold/predict_with_colabfold.wdl" as predict_with_colabfold_t

workflow ProtBindScreenSubmitMode {
	input {
		String query_name
		String query_sequence
		String subject_proteome

		File subject_native_sequences_file 
		File subject_scrambled_sequences_file 
		File subject_proteome_dict 
	}


	call prepare_fasta_groups_t.PrepareFastaGroups as t_001_prepare_fasta_groups {
		input:
			query_name = query_name,
			query_sequence = query_sequence,
			subject_proteome = subject_proteome,
			subject_native_sequences_file = subject_native_sequences_file,
			subject_scrambled_sequences_file = subject_scrambled_sequences_file,
			subject_proteome_dict = subject_proteome_dict
	}

	#call mock_task_t.SeeHowInputVarLook as t_002_see_how_input_var_look {
	#	input:
	#		predictions_input_files = t_001_prepare_fasta_groups.predictions_input_files,
	#		screen_dir_string_path = t_001_prepare_fasta_groups.screen_dir_string_path,
	#		size_tags = size_tags
	#}

	
	scatter (fasta_file in t_001_prepare_fasta_groups.predictions_input_files) {
		call predict_with_colabfold_locally_t.PredictWithColabfoldLocally as t_002_predict_with_colabfold_locally {
			input:
				fasta_file = fasta_file
		}
	}

	#call predict_with_colabfold_locally_t.PredictWithColabfoldLocally as t_002_predict_with_colabfold_locally {
	#	input:
	#		predictions_input_files = t_001_prepare_fasta_groups.predictions_input_files
	#}

	output {
		Array[File] predictions_input_files = t_001_prepare_fasta_groups.predictions_input_files
		Array[File] predictions_input_counts = t_001_prepare_fasta_groups.predictions_input_counts
		Array[File] placement_data = t_001_prepare_fasta_groups.placement_data
		String screen_dir_string_path = t_001_prepare_fasta_groups.screen_dir_string_path
		Array[File] colabfold_outputs = flatten(t_002_predict_with_colabfold_locally.colabfold_outputs)
		#Array[File] colabfold_outputs = t_002_predict_with_colabfold_locally.colabfold_outputs
		#File mock_report = t_002_see_how_input_var_look.mock_report
	}
}
