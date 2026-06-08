version 1.0

import "../Tasks/analyze_screen/analyze_screen.wdl" as analyze_screen_t

workflow ProtBindScreenAnalyzeMode {
	input {
		String bucket_name
		String screen_id
		File fasta_input_inventory 
		File colabfold_output_inventory
		Array[File] subject_proteome_datasets 
		String analysis_id
		Int query_len 
		Array[String] aa_ranges_i
		Array[String] aa_ranges_j
		Int pae_threshold
	}

	call analyze_screen_t.AnalyzeScreen as t_001_analyze_screen {
		input:
		    bucket_name = bucket_name,
		    screen_id = screen_id,
		    fasta_input_inventory = fasta_input_inventory,
		    colabfold_output_inventory = colabfold_output_inventory,
		    subject_proteome_datasets = subject_proteome_datasets,
		    analysis_id = analysis_id,
		    query_len = query_len,
		    aa_ranges_i = aa_ranges_i,
		    aa_ranges_j = aa_ranges_j,
		    pae_threshold = pae_threshold
	}
	output {
		File analysis_output_inventory = t_001_analyze_screen.analysis_output_inventory
	}
}
