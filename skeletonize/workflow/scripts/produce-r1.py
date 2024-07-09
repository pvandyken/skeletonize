from pathlib import Path
import numpy as np
import nibabel as nib
from snakeboost import snakemake_args


if __name__ == "__main__":
    args = snakemake_args()
    assert isinstance(args.input, dict)
    assert isinstance(args.output, list)

    data_path = args.input["data"]
    mask_path = args.input["mask"]
    output = Path(args.output[0])

    t1map = nib.load(data_path)
    data = t1map.get_fdata()
    mask_data = nib.load(mask_path).get_fdata().astype(bool)
    masked = np.ma.masked_where(~mask_data, data)
    r1_data = (1/masked).data
    r1_data[~mask_data] = 0

    r1 = nib.Nifti1Image(r1_data, t1map.affine)
    nib.save(r1, output)
