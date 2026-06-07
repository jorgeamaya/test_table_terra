version 1.0

task TestEnvironment {

	input {
		String bucket_name
	}

	command <<<
		set -e

		BUCKET_URI="gs://~{bucket_name}"

		echo "~{bucket_name}" > bucket_name.txt
		echo "${BUCKET_URI}" > bucket_uri.txt
	>>>

	output {
		String bucket_name_out = read_string("bucket_name.txt")
		String bucket_uri_out = read_string("bucket_uri.txt")
	}

	runtime {
		docker: "ubuntu:22.04"
	}
}
