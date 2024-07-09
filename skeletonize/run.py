#!/usr/bin/env python3
import os
from snakebids.app import SnakeBidsApp
from snakebids.exceptions import ConfigError



def pick_component(app: SnakeBidsApp):
    param_map = app.config["pybids_inputs"]["param_map"]
    app.config["pybids_inputs_run"] = {
        "fa": app.config["pybids_inputs"]["fa"]
    }
    if "filters" in app.config["pybids_inputs"]["param_map"]:
        app.config["pybids_inputs_run"]["param_map"] = param_map

def main():
    app = SnakeBidsApp(
        os.path.abspath(os.path.dirname(__file__)),
        plugins=[
            pick_component,
        ],
    )
    app.run_snakemake()


if __name__ == "__main__":
    main()