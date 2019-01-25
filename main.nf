#!/usr/bin/env nextflow

/*
===========
GEM Smasher
===========

Pipeline for extracting sub-networks from Gene Expression
Matricies (GEM) files.
*/


// Define a help message to be displayed to the user.
def show_help() {
  log.info """

  Usage:

    nextflow run main.nf --gem 'gem_file.csv'

  Required Arguments:

    --gem             Path to the input gene expression matrix, this file
                      should be `csv` formatted, and the first column should
                      be gene indexes, and the first row should be sample
                      indexes.

  Optional Arguments:

    --tag             Optional name to be used in the construction of output
                      filenames, this defaults to 'GEM'.

    --sample_size     Number of features to sample from the gem file.

    --output_dir      Directory to save results to, by default this will be a
                      folder named 'output' created in the directory
                      GEM Smasher was launched from.

  """
}


// Define a helper function that will construct (if needed)
// arguments to be passed on to Python scripts.
def prep_kwargs(param_args) { param_args ? "--key_pairs ${param_args}" : "" }


// Begin displaying stuff.
log.info"""
===============================================================================
                ┏━━━┓                             ╭□
        ◿ □ ◺   ┃   ╠══════════════════╗       ╭□ ┤
        □ ◪ ■   ┃   ║  G.E.M. Smasher  ║     ◪ ┤  ╭■
        ◹ ■ ◤   ┃   ╠═════╤════╤═══════╝       ╰■ ┤
                ┗━━━┛     ┆    ┆                  ╰◤

===============================================================================

GEM to smash:       ${params.gem}
Use subset:         ${params.sample_size.toString().isInteger()}
Subset size:        ${params.sample_size}

"""



// Check if the help message should be displayed.
if (params.help || !params.gem){ show_help() exit 0 }

// Create a channel for the input file (the original GEM).
GEM_ORIGINAL = Channel
  .fromPath(params.gem)
  .map { file -> tuple(params.tag, file) }

// Manage the optional skipping of the subset process.
// See official nextflow examples for a discussion of this code.
(FULL_GEM_A, FULL_GEM_B) = ( params.sample_size.toString().isInteger()
                         ? [GEM_ORIGINAL, Channel.empty()]
                         : [Channel.empty(), GEM_ORIGINAL])


// Prepare the logging configuration for all processes.
// This may disallow the 'task.ext.logging' value from being null,
// perhaps it can only be an empty string.
// TODO: Consider refactoring this function, it is very ad-hoc and
//       unlikely to be very robust to configuration changes.
def prep_smasherPy_log_args(task) {
  ( task.ext.log_smasherPy
  ? "${task.ext.logging} --py_log_name ${task.ext.python_log_file}"
  : "${task.ext.logging}")
}



/*******************************************************************************
Subsample GEM.

An optional process to speed up testing.
*******************************************************************************/
process Subsample_GEM {

  tag { cluster_ID }
  publishDir "${params.output_dir}/gem_subsamples/"

  input:
    set val(cluster_ID), file(gem) from FULL_GEM_A

  output:
    set val(cluster_ID), file("sample_${cluster_ID}.csv") into SAMPLE_GEM

  // Only run this process when a user provides an integer.
  when: params.sample_size

  script:
  """
  smasher.py ${prep_smasherPy_log_args(task)} \
    -i $gem -t $cluster_ID \
  subsample-gem --size ${params.sample_size}
  """
}



/*******************************************************************************
Normalize the input GEM.

This process only operates on the original GEM.
*******************************************************************************/
process Normalize_GEM {

  tag { cluster_ID }
  publishDir "${params.output_dir}/normalized_gems/"

  input:
    set val(cluster_ID), file(gem) from FULL_GEM_B.mix(SAMPLE_GEM)

  output:
    set val(cluster_ID), file("normal_${cluster_ID}.csv") into NORMALIZED_GEM

  script:
  """
  smasher.py ${prep_smasherPy_log_args(task)} \
    -i $gem -t $cluster_ID \
  normalize-gem --method ${params.normalization_method} \
    ${prep_kwargs("${task.ext.normal_kwargs}")}
  """
}



/*******************************************************************************
File and hyper-parameter channel creation.

For now only basic file input is implemented.
*******************************************************************************/
// Create a channel for the original input file.
// The Python script will emit "None" when it can no longer create
// subsets, the `until` below captures this and ends this channel.
NORMALIZED_GEM_SUBSETS = Channel.create()
NORMAL_GEMS = NORMALIZED_GEM
  .mix(NORMALIZED_GEM_SUBSETS)
  .until { it=="None" }



/*******************************************************************************
Dimensinoal Reduction.

Future:
+ More methods, eg. tSNE, PCA.
+ Hyper-parameter optimization.
*******************************************************************************/
process Dimensional_Reduction {

  tag { cluster_id }
  publishDir "${params.output_dir}/reduction/"

  input:
    set val(cluster_id), file(normal_gem) from NORMAL_GEMS

  output:
    set val(cluster_id), file("reduction_*.csv") into REDUCED_GEM

  script:
  """
  smasher.py ${prep_smasherPy_log_args(task)} \
    -i ${normal_gem} -t ${cluster_id} \
  run-umap \
    ${prep_kwargs("${task.ext.reduction}")}
  """
}



/*******************************************************************************
Clustering of reduced dimension embeddings.

For now HDBSCAN is used.

Future:
+ Explore other clustering methods.
+ Use tree-cluster feature of HDBSCAN.
*******************************************************************************/
process Cluster {

  tag { cluster_id }

  input:
    set val(cluster_id), file(reduced_gem) from REDUCED_GEM

  output:
    set val(cluster_id), file("clusters_*.csv") into CLUSTER_LABELS_A

  script:
  """
  smasher.py ${prep_smasherPy_log_args(task)} \
    -i $reduced_gem -t $cluster_id \
  cluster \
    ${prep_kwargs("${task.ext.cluster}")}
  """
}



/*******************************************************************************
Cluster Subset Creation.

Consider the output shape of these subsets with respect to re-using
embedding models.
*******************************************************************************/
process Subset_GEM {

  tag { cluster_id }

  input:
    set val(cluster_id), file(cluster_labels) from CLUSTER_LABELS_A
    // Use a ternary operator to select either the complete GEM,
    // or the subset created earlier.
    file(original_gem) from file(
      params.sample_size.toString().isInteger()
      ? "${params.output_dir}/gem_subsamples/sample_${params.tag}.csv"
      : "${params.gem}"
    )

  output:
    file("subset_*.csv") into ALL_GEM_SUBSETS mode flatten

  script:
  """
  smasher.py ${prep_smasherPy_log_args(task)} \
    -i ${cluster_labels} -t ${params.tag} \
  subset-gem \
    --originalGEM ${original_gem} \
    --minsize ${params.min_subset_size} \
    --max_depth ${params.max_depth} \
    ${prep_kwargs("${task.ext.subset}")}
  """
}


GEM_SUBSETS = ALL_GEM_SUBSETS.map {
  file -> tuple("${params.tag}_${file.baseName}".replace("subset_", ""), file)
}



// /*******************************************************************************
//
// *******************************************************************************/
// process Score_Cluster {
//
//   tag { cluster_id }
//
//   publishDir "${prams.output_dir}/p_values/"
//
//   input:
//   set val(cluster_id), file(cluster_labels) from CLUSTER_LABELS_B
//
//   output:
//   file("*_scores.csv")
//
//   """
//   python3 ${PWD}/scripts/main.py score-cluster $cluster_labels $cluster_id
//   """
//
// }



/*******************************************************************************
Normalize GEM subsets.

Re-use the original normalization process? It may be easier to simply keep
them separate.
*******************************************************************************/
process Normalize_GEM_Subset {

  publishDir "${params.output_dir}/normal_gem_subset/"
  tag { cluster_id }

  input:
    set val(cluster_id), file(gem_subset) from GEM_SUBSETS

  output:
    set val(cluster_id), file("normal_*.csv") into NORMALIZED_GEM_SUBSETS

  script:
  """
  smasher.py ${prep_smasherPy_log_args(task)} \
    -i ${gem_subset} -t ${cluster_id} \
  normalize-gem-subset \
    --method ${params.normalization_method} \
    ${prep_kwargs("${task.ext.normal_kwargs}")}
  """
}



/*******************************************************************************

*******************************************************************************/
// process ML_Model_Subset {
//
//   publishDir "${PWD}/model_ranks/"
//
//   input:
//   file(cluster_labels) from GEM_SUBSET_B
//
//   output:
//   file(gene_ranks)
//
//   """
//   python3 ${PWD}/scripts/rfe_model.py $cluster_labels
//   """
//
// }
