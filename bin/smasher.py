#!/usr/bin/env python3

"""Smasher.py is a collection of functions used in gem-smasher.

"""


# ----------------------------------------------------------------------------
# General imports.
# ----------------------------------------------------------------------------
import os
import sys
import logging
import argparse


# ----------------------------------------------------------------------------
# Data science imports.
# ----------------------------------------------------------------------------
import umap
import hdbscan
import pandas as pd
from sklearn.preprocessing.data import QuantileTransformer
from sklearn.preprocessing import StandardScaler, MinMaxScaler, MaxAbsScaler
from scipy.sparse import csc_matrix


# ----------------------------------------------------------------------------
# Define normalization methods.
# ----------------------------------------------------------------------------
METHODS = {'standard_scaling': StandardScaler,
           'min_max_scaling': MinMaxScaler,
           'min_abs_scaling': MaxAbsScaler,
           'guassian_quantile': QuantileTransformer}


# ----------------------------------------------------------------------------
# Argparse to kwargs implementation.
# https://stackoverflow.com/a/42355279
# ----------------------------------------------------------------------------
class StoreDictKeyPair(argparse.Action):
    """Allows for passing undefined kwargs to functions.

    This prevents the need to construct argparse arguments for functions
    to be called within from the command line, or by Nextflow.
    """

    def __init__(self, option_strings, dest, nargs=None, **kwargs):
        self._nargs = nargs
        super().__init__(option_strings, dest, nargs=nargs, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        my_dict = {}
        for kv in values:
            key, value = kv.split('=')
            my_dict[key] = value
        setattr(namespace, self.dest, my_dict)


# ----------------------------------------------------------------------------
# Parser Configuration.
# ----------------------------------------------------------------------------
parser = argparse.ArgumentParser(description=__doc__)
# Create a subparser for the individual steps / functions.
subparsers = parser.add_subparsers(help='Subfunctions to call.')
# All of these functions will take an input file and a tag / label.
parser.add_argument('-i', '--infile', type=str)
parser.add_argument('-t', '--tag', type=str)

# Setup parsing for logging level.
parser.add_argument('-d', '--debug', help="Print lots of debugging statements",
                    action="store_const", dest="loglevel", const=logging.DEBUG,
                    default=logging.WARNING)
parser.add_argument('-v', '--verbose', help="Be verbose",
                    action="store_const", dest="loglevel", const=logging.INFO)


# ----------------------------------------------------------------------------
# Interpret filename.
# ----------------------------------------------------------------------------

def get_filename(path):
    """Return the filename without its path or extension."""
    return os.path.splitext(os.path.basename(path))[0]


def interp_cluster_name(path, tag):
    """(parent.GEM)_(Top Cluster ID)[_Next cluster ID]"""
    # Split by the .GEM.
    filename = get_filename(path)
    _parent, cluster_string = filename.split(tag)
    # Split by underscore.
    clusters = tuple(cluster_string.split("_"))
    # Return the rank(depth) of this cluster, a tuple of parent IDs,
    # as well as the name of its parent.
    return len(clusters), clusters


def load_original_gem(gem_path):
    """Load the original GEM, as created by GEMmaker."""
    return pd.read_table(gem_path, index_col=0).T


# ----------------------------------------------------------------------------
# Main Functions.
# ----------------------------------------------------------------------------
def normalize_gem(args: dict):
    """Normalize a .GEM file with a Quantile transformer."""
    # Read the input file.
    df = load_original_gem(args['infile'])

    # Convert the data set to a sparse matrix.
    # csc_matrix is used as it is converted to that format in use.
    sparse_df = csc_matrix(df)

    # Create the normalization model.
    model = METHODS[args['method']]
    model_kwargs = args.get('key_pairs')

    if model_kwargs is not None:
        model = model(**model_kwargs)
    else:
        model = model()

    logging.debug(f"Running the normalization: {model}"
                  f"with the arguments: {model_kwargs}")

    # Run the transform, and convert the data back into a dataframe.
    transform = model.fit_transform(sparse_df)
    output_df = pd.DataFrame(transform.toarray(), index=df.index)

    # Ensure the index has the same name for consistency.
    output_df.index.name = "Sample_ID"

    # Create the output file name based on the provided tag.
    outfile = f"normal_{args['tag']}.csv"
    # Save the file and log the results.
    output_df.to_csv(outfile)
    logging.debug(f"Transform saved to: {outfile}")


normalize_gem_parser = subparsers.add_parser("normalize-gem")
normalize_gem_parser.add_argument('--method', action='store', type=str)
normalize_gem_parser.add_argument("--key_pairs", dest="my_dict",
                                  action=StoreDictKeyPair, nargs="+")
normalize_gem_parser.set_defaults(func=normalize_gem)


# ----------------------------------------------------------------------------


def subsample_gem(args: dict):
    """Subsample the original GEM matrix."""
    logging.debug(f"Taking subsample of size: {args['size']}")
    df = load_original_gem(args['infile'])

    sub_df = df.sample(int(args['size']), axis=1)
    outfile = f"sample_{args['tag']}.csv"
    sub_df.to_csv(outfile, sep='\t')


subsample_parser = subparsers.add_parser("subsample-gem")
subsample_parser.add_argument("--size", type=int)
subsample_parser.set_defaults(func=subsample_gem)


# ----------------------------------------------------------------------------


def run_umap(args: dict):
    """Reduce dimensions of a given file."""
    # Read the input file.
    df = pd.read_csv(args['infile'], index_col=0)

    # Create the reduction model.
    model = umap.UMAP
    model_kwargs = args.get('key_pairs')

    if model_kwargs is not None:
        model = model(**model_kwargs)
    else:
        model = model()

    # Create and save the embedding.
    reduction = model.fit_transform(df)
    output_df = pd.DataFrame(reduction, index=df.index)
    outfile = f"reduction_{args['tag']}.csv"
    logging.debug(f"Transform object type: {type(reduction)}.")
    logging.debug(f"Saving reduction to: {outfile}")
    output_df.to_csv(outfile)


reduce_dimensions_parser = subparsers.add_parser("run-umap")
reduce_dimensions_parser.add_argument("--key_pairs", dest="my_dict",
                                      action=StoreDictKeyPair, nargs="+")
reduce_dimensions_parser.set_defaults(func=run_umap)


# ----------------------------------------------------------------------------


def cluster(args: dict):
    """Cluster the given embedding."""
    # Read the input file.
    df = pd.read_csv(args['infile'], index_col=0)

    # Create the clustring model, fit the data and save the results.
    model_kwargs = args.get('key_pairs')

    if model_kwargs is not None:
        model = hdbscan.HDBSCAN(**model_kwargs)
    else:
        model = hdbscan.HDBSCAN()

    model.fit(df)
    outfile = f"clusters_{args['tag']}.csv"
    output_df = pd.DataFrame({"cluster_ID": model.labels_}, index=df.index)
    output_df.to_csv(outfile)
    logging.debug(f"Saved cluster labels to: {outfile}")


cluster_parser = subparsers.add_parser("cluster")
cluster_parser.add_argument("--key_pairs", dest="my_dict",
                            action=StoreDictKeyPair, nargs="+")
cluster_parser.set_defaults(func=cluster)


# ----------------------------------------------------------------------------


def subset_gem(args):
    """Create a subset of a GEM file based on the given cluster labels."""
    # Read in the orginal GEM and the given label file.
    df = load_original_gem(args['originalGEM'])
    labels = pd.read_csv(args['infile'], index_col=0)

    # Determine the rank and ID of the current of the subset.
    rank, ids = interp_cluster_name(args['infile'], args['tag'])
    new_rank = rank + 1

    # Ensure the tree is not too deep.
    if rank >= args['max_depth']:
        logging.info(f"Rank limit of {args['max_depth']} reached.")
        sys.stdout.write("None")
        exit(0)

    # Add the clustring group IDs to the dataframe.
    df[f"cluster_ID_{new_rank}"] = labels["cluster_ID"]

    # Group the dataframe by cluster ID.
    cluster_groups = df.groupby(f"cluster_ID_{new_rank}")
    logging.debug(f"Examining {len(cluster_groups)} groups.")

    # Groupby the new cluster ID and save each group as a separate .csv.
    for key, group in cluster_groups:
        # Pull out the subset, and ensure the subset is not too small.
        # We also ignore un-clustered items (labeled -1).
        subset = df.loc[group.index.values]
        if subset.shape[0] < args['minsize'] or key == -1:
            continue

        # Build the new cluster ID and output file name, then save.
        outfile = f"subset_{new_rank}.csv"
        logging.debug(f"New subset: {new_rank} saved to: {outfile}")
        subset.to_csv(outfile)


subset_gem_parser = subparsers.add_parser("subset-gem")
subset_gem_parser.set_defaults(func=subset_gem)
subset_gem_parser.add_argument('--originalGEM', action='store', type=str)
subset_gem_parser.add_argument('--minsize', action='store', type=int)
subset_gem_parser.add_argument('--max_depth', action='store', type=int)
subset_gem_parser.add_argument("--key_pairs", dest="my_dict",
                               action=StoreDictKeyPair, nargs="+")

# ----------------------------------------------------------------------------


def score_cluster(args):
    """Place holder function for now..."""
    rank, ids = interp_cluster_name(args['infile'], args['tag'])
    outfile = f"{rank}_{'_'.join(ids)}_scores.csv"
    empty_data = pd.DataFrame(
        {f"subset_{args['tag']}_{'_'.join(ids)}": [0, 0, 0]})
    # TODO: Some code to rank clusters...
    empty_data.to_csv(outfile)


score_cluster_parser = subparsers.add_parser("score-cluster")
score_cluster_parser.set_defaults(func=score_cluster)
score_cluster_parser.add_argument("--key_pairs", dest="my_dict",
                                  action=StoreDictKeyPair, nargs="+")

# ----------------------------------------------------------------------------


def normalize_gem_subset(args):
    """Normalize a GEM subset."""
    df = pd.read_csv(args['infile'], index_col=0)

    cluster_id = df.iloc[0, -1]

    new_tag = f"{args['tag']}_{cluster_id}"
    outfile = f"normal_{new_tag}.csv"
    logging.debug(f"Normalizing: {args['infile']} with tag: {new_tag}.")

    sparse_df = csc_matrix(df)
    model = METHODS[args['method']]()
    transform = model.fit_transform(sparse_df)
    output_df = pd.DataFrame(transform.toarray(), index=df.index)
    output_df.index.name = "Sample_ID"
    output_df.to_csv(outfile)
    sys.stdout.write(new_tag)
    logging.debug(f"Transform saved to: {outfile}")


normalize_parser = subparsers.add_parser("normalize-gem-subset")
normalize_parser.set_defaults(func=normalize_gem_subset)
normalize_parser.add_argument('--method', action='store', type=str,
                              choices=list(METHODS.keys()))
normalize_parser.add_argument("--key_pairs", dest="my_dict",
                              action=StoreDictKeyPair, nargs="+")


# ----------------------------------------------------------------------------
# Code to run if this file is called as a script.
# ----------------------------------------------------------------------------
if __name__ == '__main__':
    # Parse and log the arguments from the command line.
    given_arguments = vars(parser.parse_args())

    # Set up the logger.
    logging.basicConfig(stream=sys.stderr, level=given_arguments['loglevel'])

    # Pass a dictionary of provided command line arguments as a dictionary
    # to the user-selected function.
    given_arguments['func'](given_arguments)
