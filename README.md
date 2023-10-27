# PaintSHOP DNA QC
R script for assessing the uniqueness of DNA FISH oligo probes' mapping to dm6 genome. Probes are generated in [PaintSHOP](https://paintshop.io/) and loaded to the 'data' folder. Run the script 'probeQC.R' which will produce these plots in 'plots' directory: [karyoplot](plots/all_probes_karyo.pdf) and [genomic probe counts](plots/all_probes_counts_barps.pdf)


## Requirements
* [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) and [mamba](https://mamba.readthedocs.io/en/latest/mamba-installation.html)
* conda environment 'align' with bowtie2 and samtools installed
  ```
  conda config --add channels defaults
  conda config --add channels bioconda
  conda config --add channels conda-forge
  conda config --set channel_priority strict
  mamba create -n align -c bioconda bowtie2 samtools
  ```
