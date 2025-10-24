version 1.0

import "../tasks/analyze_screen/analyze_screen.wdl" as analyze_screen_t


workflow ProtBindScreenAnalyzeMode {
	input {
        String screen_name #=name_of_the_screen_run
        String analysis_name #=name_of_the_analysis_screen_run
	Array[File]  all_fasta_files_and_the_prediction_outputs #they can be just a long array, really doesn't matter because we have all the info in the file names
	File subject_proteome_dictionary #=subject_proteome_dictionary.tsv 
        Int query_len #keep name query_len for compatibility with analysis script
        Map[String, Map[String, String]] analysis_matrices #=analysis_matrices.json that mat string mat string string was just suggested by wdl syntax highlighter but maybe it's good I leave it for you
	}

	call analyze_screen_t.AnalyzeScreen as t_001_analyze_screen {
		input:
            screen_name = screen_name,
            analysis_name = analysis_name,
            all_fasta_files_and_the_prediction_outputs = all_fasta_files_and_the_prediction_outputs, 
            subject_proteome_dictionary = subject_proteome_dictionary,
            query_len = query_len,
            analysis_matrices = analysis_matrices
	}
	output {
		Array[File] analysis_output_files = t_001_analyze_screen.analysis_output_files
	}
}
