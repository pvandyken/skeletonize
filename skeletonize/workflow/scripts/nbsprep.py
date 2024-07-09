from math import log10, floor
import numpy as np
import pandas as pd
import shutil as sh
from scipy.io import savemat
from snakeboost import snakemake_args


def get_metadata(file, hem = None):
    metadata = pd.read_csv(file, sep="\t")
    if hem in ("L", "R"):
        return metadata[metadata["hemisphere"] == hem].reset_index()
    return metadata.reset_index()

    
def get_hemisphere_ix_(metadata_file, hem):
    index = get_metadata(metadata_file, hem)["Label ID"]
    return np.ix_(index, index)

def get_connectomes(files):
    return np.dstack([np.loadtxt(file, delimiter=",") for file in files])

def main():
    args = snakemake_args(
        input={"connectomes": "--connectomes", "metadata": "--metadata"},
        output=["out"],
        params={
            "hemi": "--hem",
        }
    )
    if (
        not isinstance(args.input, dict) or
        'connectomes' not in args.input or
        'metadata' not in args.input
    ):
        raise TypeError("Inputs must be specified as a dict")
    if isinstance(args.output, dict):
        raise TypeError("Outputs must be specified as a single path")
    if (
        not isinstance(args.params, dict) or
        'hemi' not in args.params
    ):
        raise TypeError("Params must be specified as a dict with one key 'hemi'")


    metadata = args.input["metadata"]
    ix_ = get_hemisphere_ix_(metadata, args.params['hemi'])
    connectomes = get_connectomes(args.input['connectomes'])[(*ix_, ...)]

    savemat(args.output[0], {"data": connectomes})


if __name__ == "__main__":
    main()