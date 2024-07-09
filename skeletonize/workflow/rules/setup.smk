import os
import re
import tempfile
import functools as ft
import itertools as it

from bids.layout import parse_file_entities
import more_itertools as itx
from snakebids import (
    bids, generate_inputs, filter_list, BidsPartialComponent, set_bids_spec
)
from snakebids.exceptions import ConfigError
import pandas as pd
import numpy as np

from pathlib import Path
from snakeboost import Tar, Pyscript, XvfbRun, PipEnv, Boost, Datalad, Env
import snakeboost.bash as sh


set_bids_spec("v0_11_0")

tmpdir = eval(workflow.default_resources._args.get("tmpdir"), {"system_tmpdir": tempfile.gettempdir()})

if workflow.run_local:
    workflow.shadow_prefix = os.environ.get("SLURM_TMPDIR")

###
# Input Globals
###
inputs = generate_inputs(
    bids_dir=config['bids_dir'],
    pybids_inputs=config['pybids_inputs_run'],
    derivatives=config.get("derivatives", False),
    participant_label=config.get('participant_label'),
    exclude_participant_label=config.get('exclude_participant_label'),
    pybidsdb_dir=config.get("pybidsdb_dir"),
)

if "param_map" in inputs:
    comp = "param_map"
else:
    comp = "fa"


###
# Output Globals
###

work = Path(tmpdir) / "snakemodel"
source = Path(config['output_dir']) / "sourcedata"
output = Path(config['output_dir'])
qc = Path(output)/"qc"

# Unique ID for easy naming in temporary files

def _uid(entities = None):
    return Path(
        *sorted(
            it.chain.from_iterable(
                comp.wildcards.values() 
                if isinstance(comp, BidsPartialComponent) else [comp]
                for comp in entities
            )
        )
    )

def tempout(rulename: str, *entities, extension = ""):
    uid = _uid(entities)
    return temp(work.joinpath(rulename, uid, rulename).with_suffix(extension))

def log(rulename: str, *entities):
    uid = _uid(entities)
    return Path("code", "logs", rulename, uid, rulename).with_suffix(".log")

def benchmark(rulename: str, *entities):
    uid = _uid(entities)
    return Path("code", "benchmarks", rulename, uid,  rulename).with_suffix(".tsv")

def resource(path):
    return os.path.join(workflow.basedir, "..", "resources", path)
