#!/bin/sh
#SBATCH --partition=kamiak
#SBATCH --account=ficklin
#SBATCH --job-name=smash_gem
#SBATCH --output=smash_gem.out

# Load needed modules.
module add anaconda3/5.1.0 java nextflow

# Run the workflow.
nextflow run main.nf \
--gem="/scidas/oryza_sativa/rice_PRJNA301554_heat_drought/GEMmaker/output/GEM/rice_heat_drought.GEM.FPKM.txt" \
--sample_size=1000
