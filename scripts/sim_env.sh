## Sourced by check_inputs.sh and sweep_1d.sh: sets SIMENV to the simulation
## conda environment's folder, from `condaEnvs: sim:` in the config ($CONFIG).
##   an env file (.yml/.yaml) -> the environment Snakemake built from it
##                               (it keeps a copy of the file next to it)
##   a path                   -> that environment (e.g. one an admin built)
##   a name                   -> looked up with `conda env list`
find_sim_env() {
    local entry
    entry=$(sed -nE "s/^[[:space:]]+sim:[[:space:]]*['\"]?([^'\"#[:space:]]+).*/\1/p" "${CONFIG}" | head -1)
    entry="${entry/#\~/$HOME}"
    if [[ -z "${entry}" || "${entry}" == *.yml || "${entry}" == *.yaml ]]; then
        SIMENV=$(grep -l hictk .snakemake/conda/*.yaml 2>/dev/null | head -1 | sed 's/\.yaml$//')
        if [[ -z "${SIMENV}" || ! -x "${SIMENV}/bin/python" ]]; then
            echo "The conda environments are not built yet. On a compute node, run:" >&2
            echo "  bash ChromSimPipe.sh --conda-create-envs-only" >&2
            return 1
        fi
    elif [[ "${entry}" == */* ]]; then
        SIMENV="${entry}"
    else
        SIMENV=$(conda env list 2>/dev/null | awk -v n="${entry}" '$1 == n {print $NF}')
        if [[ -z "${SIMENV}" ]]; then
            echo "condaEnvs: sim is '${entry}', but 'conda env list' has no environment of that" >&2
            echo "name (or conda is not on PATH). Use the environment's full path in the config." >&2
            return 1
        fi
    fi
    if [[ ! -x "${SIMENV}/bin/python" ]]; then
        echo "No Python in the simulation environment ${SIMENV} (condaEnvs: sim in ${CONFIG})." >&2
        return 1
    fi
}
