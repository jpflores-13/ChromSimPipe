# data/

Input files for a run: the CTCF peak files and `.hic` maps named in
`ChromSimSamplesheet.txt`, and the genome FASTA named by `genome` in
`config/ChromSimConfig.yaml`. Everything in this folder except this README is
git-ignored.

The samplesheet and config refer to these files as `data/<file name>`, so the
same samplesheet works on every cluster. On each cluster, link that cluster's
copies in once with symbolic links, named exactly as in the samplesheet. A
link costs no disk space, so the large files stay where they are.

From the repo root:

```bash
SRC=/path/to/this/clusters/data
ln -s "$SRC"/STRS_HEK293_eGFP-YAP_CTCF_cont_0h_peaks.narrowPeak    data/
ln -s "$SRC"/STRS_HEK293_eGFP-YAP_CTCF_sorbitol_1h_peaks.narrowPeak data/
ln -s "$SRC"/YAPP_HEK293_eGFP-YAP_Cai_control_megaMap_inter_30.hic  data/
ln -s "$SRC"/YAPP_HEK293_eGFP-YAP_Cai_sorbitol_megaMap_inter_30.hic data/
ln -s /path/to/GRCh38.primary_assembly.genome.fa                    data/
ls -lL data/
```

If the files are in several folders, link each one from its own folder. If
the genome already has a `.fai` index, link it too
(`ln -s /path/to/GRCh38.primary_assembly.genome.fa.fai data/`); otherwise the
workflow indexes the genome on first use and writes the `.fai` here.

`ls -lL` follows the links, so a broken link shows up as an error. A dry run
(`bash ChromSimPipe.sh -n`) also stops with `MissingInputException`, naming
any file it cannot find.

To use files with other names, edit the samplesheet (and `genome` in the
config) to match. Absolute paths also work there.
