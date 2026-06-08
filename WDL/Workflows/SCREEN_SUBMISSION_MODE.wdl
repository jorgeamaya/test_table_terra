version 1.0

import "../Tasks/prepare_fasta_files/prepare_fasta_files.wdl" as prepare_fasta_files_t
import "../Tasks/predict_with_colabfold/predict_with_colabfold.wdl" as predict_with_colabfold_t
import "../Tasks/inventory_colabfold_outputs/inventory_colabfold_outputs.wdl" as inventory_colabfold_outputs_t

workflow ProtBindScreenSubmitMode {
	input {
		String bucket_name
		String screen_id
		String query_name
		String query_sequence
		Array[File] subject_proteome_datasets
	}

	call prepare_fasta_files_t.PrepareFastaFiles as t_001_prepare_fasta_files {
		input:
			bucket_name = bucket_name,
			screen_id = screen_id,
			query_name = query_name,
			query_sequence = query_sequence,
			subject_proteome_datasets = subject_proteome_datasets
	}

	scatter (fasta_input_group_inventory in t_001_prepare_fasta_files.fasta_input_group_inventories) {
		call predict_with_colabfold_t.PredictWithColabfold as t_002_predict_with_colabfold {
			input:
				bucket_name = bucket_name,
				screen_id = screen_id,
				fasta_input_group_inventory = fasta_input_group_inventory
		}
	}

	call inventory_colabfold_outputs_t.InventoryColabfoldOutputs as t_003_inventory_colabfold_outputs {
		input:
			bucket_name = bucket_name,
			screen_id = screen_id,
			colabfold_output_group_inventories = t_002_predict_with_colabfold.colabfold_output_group_inventory
	}

	output {
		File fasta_input_inventory = t_001_prepare_fasta_files.fasta_input_inventory
		File colabfold_output_inventory = t_003_inventory_colabfold_outputs.colabfold_output_inventory
	}
}
