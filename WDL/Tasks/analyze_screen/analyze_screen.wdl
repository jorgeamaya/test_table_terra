version 1.0

task AnalyzeScreen {
    input {
        String screen_name
        String analysis_name
        Array[File]  all_fasta_files_and_the_prediction_outputs
        File subject_proteome_dictionary_file
        Int query_len 
        Array[String] aa_ranges_i #=["1-50","51-100"] example
        Array[String] aa_ranges_j #=["1-50","51-100"] example
    }

    command <<<
    set -euxo pipefail

    conda info --envs
    python -c "import protbindscreen" && echo "Available" || echo "Not available"

    python3 - <<'PYCODE'
import json
from pathlib import Path
from protbindscreen.analysis.screen_analysis import screen_analysis

screen_dir = Path("Results/~{screen_name}")
analysis_dir = Path("Results/~{screen_name}/analysis/~{analysis_name}")
analysis_dir.mkdir(parents=True, exist_ok=True)
print("DEBUG: screen_dir =", screen_dir)
print("DEBUG: analysis_dir =", analysis_dir)

python_list_files=[Path(x) for x in "~{sep=' ' all_fasta_files_and_the_prediction_outputs}".split()]
print("DEBUG: python_list_files =", python_list_files)

# ---- Build analysis_matrices manually ----
aa_ranges_i = "~{sep=' ' aa_ranges_i}".split()
aa_ranges_j = "~{sep=' ' aa_ranges_j}".split()

analysis_matrices = {}
for i, (range_i, range_j) in enumerate(zip(aa_ranges_i, aa_ranges_j), start=0):
    i="required" if i == 0 else f"optional_{i}"
    matrix_name = f"matrix_{i}"
    analysis_matrices[matrix_name] = {
        "aa_ranges_i": range_i,
        "aa_ranges_j": range_j
    }
    if i != "required":
        analysis_matrices[matrix_name]["include"] = "true"

# Pretty print for debug
print("DEBUG: Loaded analysis_matrices (normalized):")
print(json.dumps(analysis_matrices, indent=4))

screen_analysis(
    screen_dir=screen_dir,
    analysis_dir=analysis_dir,
    analysis_name='~{analysis_name}',
    query_len=~{query_len},
    all_fasta_files_and_the_prediction_outputs=python_list_files,
    subject_proteome_dictionary=Path('~{subject_proteome_dictionary_file}'),
    analysis_matrices=analysis_matrices
)
PYCODE
    >>>
  
    output {
        Array[File] analysis_output_files = glob("Results/~{screen_name}/**/**/**")
    }

    runtime {
        cpu: 1 
        memory: "8 GiB" 
        disks: "local-disk 10 HDD" 
        bootDiskSizeGb: 10 
        preemptible: 3
        maxRetries: 1
        docker: 'us-central1-docker.pkg.dev/global-axe-475818-q0/protbindscreen-docker-repo/protbindscreen:v0.0.2'
    }
}
