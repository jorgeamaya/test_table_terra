version 1.0

import "../Tasks/test/test.wdl" as test_t

workflow TestEnvironment {

	input {
		String bucket_name
	}

	call test_t.TestEnvironment as t_001_test_environment {
		input:
			bucket_name = bucket_name
	}

	output {
		String bucket_name_out = t_001_test_environment.bucket_name_out
		String bucket_uri_out = t_001_test_environment.bucket_uri_out
	}
}
