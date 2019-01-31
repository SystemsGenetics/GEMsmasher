# GEM Smasher

Pipeline for extracting sub-networks from Gene Expression
Matrix (GEM) files.

***In Development***

*This project may not work, or give meaningful results.*

#### Author(s)
+ Tyler D. Biggs

### Testing Instructions

*A 'standard' run.*

```
nextflow run main.nf --gem="gem_file_path.GEM"
```

*Running from github.*
```
nextflow run SystemsGenetics/GEMsmasher --gem="gem_file_path.GEM"
# or
nextflow run SystemsGenetics/GEMsmasher --gem="gem_file_path.GEM" --sample_size 1000
```

*Run using a subsample of the GEM file.*

```
nextflow run main.nf --gem="gem_file_path.GEM" --sample_size 1000
```

*With Kamiak's `/scidas` mounted*

```
nextflow run main.nf \
--gem="/scidas/oryza_sativa/rice_PRJNA301554_heat_drought/GEMmaker/output/GEM/rice_heat_drought.GEM.FPKM.txt" \
--sample_size 1000
```

*In technicolor.*

```
nextflow run main.nf \
--gem="/scidas/oryza_sativa/rice_PRJNA301554_heat_drought/GEMmaker/output/GEM/rice_heat_drought.GEM.FPKM.txt" \
--sample_size=1000 | lolcat
```

*Run on Kamiak*

Submit an sbatch job that runs the workflow.

```
sbatch /scidas/tyler/GEMsmasher/scripts/smashgem.sh
```

*Cleaning up a previous run.*

This will force your `conda` image to rebuild next time you run. To
prevent this you can define the `--cacheDir` as a command line argument
to specify a different location for `conda` environments to be stored.

```
bash bin/cleanup.sh
```

### TODO

+ [x] Enable remote execution with GitHub.
+ [x] Replace basic `nextflow.config` options with CLI arguments.
+ [x] Move project to the `SystemsGenetics` GitHub repository.
+ [ ] Implement cluster scoring.
+ [ ] Implement a visualization / tree graph of the created clusters.
+ [ ] Implement variational autoencoder.



I have a process with multiple optional output files. It is part of a recursive
workflow. I am having trouble finding an elegant way to have the workflow
terminate if this process does not produce any output. Any advice?
