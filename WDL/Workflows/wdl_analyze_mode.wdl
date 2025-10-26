version 1.0

import "../Tasks/analyze_screen/analyze_screen.wdl" as analyze_screen_t

workflow ProtBindScreenAnalyzeMode {
	input {
		String screen_name 
		String analysis_name 
		Array[File] colabfold_output_files 
		Array[File] colabfold_input_files
		File subject_proteome_dictionary_file 
		Int query_len 
		Array[String] aa_ranges_i #=["1-50","51-100"] example
		Array[String] aa_ranges_j #=["1-50","51-100"] example
	}

	Array[File] all_fasta_files_and_the_prediction_outputs = flatten([colabfold_input_files, colabfold_output_files])

	call analyze_screen_t.AnalyzeScreen as t_001_analyze_screen {
		input:
			screen_name = screen_name,
			analysis_name = analysis_name,
			all_fasta_files_and_the_prediction_outputs = all_fasta_files_and_the_prediction_outputs, 
			subject_proteome_dictionary_file = subject_proteome_dictionary_file,
			query_len = query_len,
			aa_ranges_i = aa_ranges_i,
			aa_ranges_j = aa_ranges_j
	}
	output {
		Array[File] analysis_output_files = t_001_analyze_screen.analysis_output_files
		Array[String] aa_ranges_i_final = aa_ranges_i
		Array[String] aa_ranges_j_final = aa_ranges_j
	}
}
