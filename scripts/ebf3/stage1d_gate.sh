#!/bin/bash
#SBATCH --job-name=ebf3_qc_gate
#SBATCH --account=rc_dphansti_pi
#SBATCH --partition=general
#SBATCH --cpus-per-task=4
#SBATCH --mem=32g
#SBATCH --time=4:00:00
#SBATCH --output=data/logs/ebf3_qc_gate_%j.out
# ============================================================================
# Ebf3 Stage 1d go/no-go gate (submit with --dependency=afterok:<hic job>).
# Exits non-zero, which stops every afterok-dependent job (the smoke test), if:
#   * hic_qc.py fails,
#   * a 135.80 / 137.70 Mb window edge sits in a strong insulation boundary
#     (then the window must be widened to 135,467,000-137,930,000 and
#     everything downstream regenerated: not automated, needs a person),
#   * no resolution meets the pre-registered endpoint rule,
#   * check_inputs.sh reports any FAIL.
# ============================================================================
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
for s in ES ncx_NPC ncx_CN; do
    test -s "data/mcool/${s}_chr7.mcool" || { echo "GATE FAIL: missing data/mcool/${s}_chr7.mcool"; exit 1; }
done

echo "== hic_qc.py"
.conda_envs/ebf3_analysis/bin/python -I scripts/ebf3/hic_qc.py

echo "== gate decision"
.conda_envs/ebf3_analysis/bin/python -I - <<'EOF'
import json, sys
q = json.load(open("results/ebf3/hic_qc.json"))
print("endpoint_resolution_bp:", q["endpoint_resolution_bp"], "| widen_window:", q["widen_window"])
bad = []
if q["widen_window"]:
    bad.append("a window edge sits in a strong boundary -> widen to 135,467,000-137,930,000 (manual)")
if not q["endpoint_resolution_bp"]:
    bad.append("no resolution meets the pre-registered E-P count rule")
if bad:
    print("GATE FAIL:", "; ".join(bad)); sys.exit(1)
print("GATE PASS (Hi-C QC)")
EOF

echo "== check_inputs.sh (ebf3 config)"
bash check_inputs.sh --configfile config/ChromSimConfig.ebf3.yaml | tee results/ebf3/check_inputs_ebf3.txt
if grep -q "FAIL" results/ebf3/check_inputs_ebf3.txt; then echo "GATE FAIL: check_inputs"; exit 1; fi
echo "GATE PASS (all)"
