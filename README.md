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

*Run using a subsample of the GEM file.*

```
nextflow run main.nf --gem="gem_file_path.GEM" \
--sample_size 1000
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

### TODO

+ [ ] Enable remote execution with github.
+ [x] Replace basic `nextflow.config` options with CLI arguments.
+ [ ] Move project to the SystemsGenetics Github repository.
+ [ ] Implement cluster scoring.
+ [ ] Implement a visualization / tree graph of the created clusters.
