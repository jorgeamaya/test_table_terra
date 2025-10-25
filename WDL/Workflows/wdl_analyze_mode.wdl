version 1.0

import "../Tasks/analyze_screen/analyze_screen.wdl" as analyze_screen_t

workflow ProtBindScreenAnalyzeMode {
	input {
        String screen_name #=name_of_the_screen_run
        String analysis_name #=name_of_the_analysis_screen_run
	Array[File] prediction_outputs #they can be just a long array, really doesn't matter because we have all the info in the file names
        Array[File] all_fasta_files
	File subject_proteome_dictionary_file #=subject_proteome_dictionary.tsv 
        Int query_len #keep name query_len for compatibility with analysis script
	Array[String] aa_ranges_i #=["1-50","51-100"] example
	Array[String] aa_ranges_j #=["1-50","51-100"] example
	}

	Array[File] all_fasta_files_and_the_prediction_outputs = flatten([all_fasta_files, prediction_outputs])
	
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
	}
}
