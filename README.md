# GEM Smasher

Pipeline for extracting sub-networks from Gene Expression
Matrix (GEM) files.

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

*Cleaning up a previous run.*

This will force your `conda` image to rebuild next time you run. To
prevent this you can define the `--cacheDir` as a command line argument
to specify a different location for `conda` environments to be stored.

```
bash bin/cleanup.sh
```

### TODO

+ [x] Enable remote execution with github.
+ [x] Replace basic `nextflow.config` options with CLI arguments.
+ [x] Move project to the SystemsGenetics Github repository.
+ [ ] Implement cluster scoring.
+ [ ] Implement a visualization / tree graph of the created clusters.
