import itertools as it
import numpy as np
import scipy
import h5py

from snakeboost import snakemake_args

def unpack_mat(mat, scheme):

    result = {}
    for key, val in scheme.items():
        if isinstance(val, dict):
            result[key] = unpack_mat(mat[key][0,0], val)
        elif val == "literal":
            try:
                art = mat[key].flatten()
                if len(art) > 1:
                    raise IndexError("Multiple values found")
                elif len(art):
                    result[key] = art[0]
                else:
                    result[key] = ''
            except Exception:
                print(key)
                print(mat[key])
                raise
        elif val == "arr":
            result[key] = mat[key]
        else:
            raise TypeError()
    return result

def convert_to_h5py():
    snakemake = snakemake_args()
    file = snakemake.input["matrix"]
    results = scipy.io.loadmat(file)

    UI_STRUC = {"ui": "literal", "ok": "literal"}   
    nbs = unpack_mat(results["nbs"][0,0], {
        "NBS": {
            "n": "literal",
            "con_mat": "arr",
            "pval": "arr",
            "test_stat": "arr",
        },
        "GLM": {
            "y": "arr",
            "X": "arr",
            "contrast": "arr",
            "test": "literal",
            "perms": "literal",
        },
        "STATS": {
            "thresh": "literal",
            "alpha": "literal",
            "size": "literal",
            "N": "literal",
            "test_stat": "arr"
        },
        "UI": dict(zip([
            "method",
            "test",
            "size",
            "thresh",
            "perms",
            "alpha",
            "contrast",
            "design",
            "matrices",
            "node_coor",
            "node_label",
            "exchange",
        ], it.repeat(UI_STRUC)))
        # "UI": v["UI"],
        # "STATS": v["STATS"],
    })

    # np.savetxt(file.with_suffix(".tsv"), nbs["NBS"]["test_stat"], delimiter="\t")
    with h5py.File(snakemake.output["matrix"], 'w') as f:
        # for i in range(len(nbs["NBS"]['con_mat'])):
        #     print(nbs["NBS"]["con_mat"][i])
        conmats = [
            nbs["NBS"]['con_mat'][0, i].A for i in range(nbs["NBS"]['con_mat'].shape[1])
        ]
        conmats = np.dstack(conmats) if conmats else np.ndarray((0,))
        NBS = f.create_group("nbs")
        NBS["con_mat"] = conmats
        NBS["test_stat"] = nbs["NBS"]["test_stat"]
        NBS.attrs['n'] = nbs["NBS"]['n']
        NBS.attrs['pval'] = nbs["NBS"]['pval']

if __name__ == "__main__":
    convert_to_h5py()
