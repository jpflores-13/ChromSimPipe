#!/bin/bash
#SBATCH --job-name=bdnf_gate
#SBATCH --account=rc_dphansti_pi
#SBATCH --partition=general
#SBATCH --cpus-per-task=4
#SBATCH --mem=32g
#SBATCH --time=2:00:00
#SBATCH --output=data/logs/bdnf_gate_%j.out
# Bdnf go/no-go gate (submit with --dependency=afterok:<chr2 Hi-C job>). Fails, which
# stops the afterok-dependent production run, if a window edge sits in a strong
# boundary or check_inputs.sh reports any FAIL.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
for s in ES ncx_NPC ncx_CN; do
    test -s "data/mcool/${s}_chr2.mcool" || { echo "GATE FAIL: missing data/mcool/${s}_chr2.mcool"; exit 1; }
done
echo "== window edges"
.conda_envs/ebf3_analysis/bin/python -I scripts/ebf3/window_edge_check.py --chrom chr2 \
    --start 108800000 --end 110300000 --mcool 'data/mcool/{state}_chr2.mcool' \
    --out results/bdnf/window_edge_check.json --fig figures/bdnf/02_insulation_window_check.png \
    || { echo "GATE FAIL: a window edge sits in a strong boundary (see results/bdnf/window_edge_check.json)"; exit 1; }
echo "== check_inputs.sh"
bash check_inputs.sh --configfile config/ChromSimConfig.bdnf.yaml | tee results/bdnf/check_inputs_bdnf.txt
if grep -q "FAIL" results/bdnf/check_inputs_bdnf.txt; then echo "GATE FAIL: check_inputs"; exit 1; fi
echo "GATE PASS"
