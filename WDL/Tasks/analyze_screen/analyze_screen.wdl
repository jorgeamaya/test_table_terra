version 1.0

task AnalyzeScreen {
    input {
        String screen_name
        String analysis_name
        Array[File]  all_fasta_files_and_the_prediction_outputs
	File subject_proteome_dictionary
        Int query_len 
        Map[String, Map[String, String]] analysis_matrices 
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

analysis_matrices_path = Path("~{write_json(analysis_matrices)}")
print("DEBUG: analysis_matrices path =", analysis_matrices_path)
with open(analysis_matrices_path, "r") as f:
    analysis_matrices = json.load(f)

# Convert from WDL-style [{"left": k, "right": v}] to normal dict
if isinstance(analysis_matrices, list) and all("left" in item and "right" in item for item in analysis_matrices):
    analysis_matrices = {item["left"]: item["right"] for item in analysis_matrices}

# Pretty print for debug
print("DEBUG: Loaded analysis_matrices (normalized):")
print(json.dumps(analysis_matrices, indent=4))

screen_analysis(
    screen_dir=screen_dir,
    analysis_dir=analysis_dir,
    analysis_name='~{analysis_name}',
    query_len=~{query_len},
    all_fasta_files_and_the_prediction_outputs=python_list_files,
    subject_proteome_dictionary=Path('~{subject_proteome_dictionary}'),
    analysis_matrices=analysis_matrices
)
PYCODE
    >>>
  
    output {
        Array[File] analysis_output_files = glob("Results/~{screen_name}/**")
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
