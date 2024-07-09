from collections import defaultdict
from os import PathLike
from pathlib import Path
import pandas as pd
import more_itertools as itx
import functools as ft

def _set_subsess_index(df):
    index = [df["participant_id"].map(lambda s: s[4:])]
    drop = ["participant_id"]
    if "session_id" in df:
        index.append(df["session_id"].map(lambda s: s[4:]))
        drop.append("session_id")
    return df.set_index(index, drop=True).drop(columns=drop)

def load_metadata(
    participant_path: PathLike[str | bytes],
    phenotypes_path: PathLike[str | bytes] | None = None,
    phenotypes_spec: object | None = None,
):
    sub_df = load_participants(participant_path)
    if phenotypes_path is None or phenotypes_spec is None:
        return sub_df
    phenotypes = load_phenotypes(phenotypes_path, phenotypes_spec)
    return ft.reduce(
        lambda df1, df2: df1.join(df2, how="outer"), phenotypes, sub_df
    )

def load_participants(participant_path):
    df = pd.read_csv(participant_path, sep="\t")
    if "session_id" in df:
        raise ValueError("participants.tsv should not have a session_id column")
    return _set_subsess_index(df)

# def filter_participants(participant_path, **filters):
#     df = pd.read_csv(participant_path, sep="\t")
#     for col, values in filters.items():
#         if col not in df:
#             continue
#         df = df[df[col].isin(map(df[col].dtype.type, itx.always_iterable(values)))]
#     return df.assign(
#         participant_id=df["participant_id"].map(lambda s: s[4:]),
#         session_id=df["session_id"].map(lambda s: s[4:])
#     )

def load_phenotypes(phen_path, spec):
    by_file = defaultdict(dict)
    for field, loc in spec.items():
        by_file[loc["file"]][field] = loc.get("field", field)
    dfs = []
    for filename, field_map in by_file.items():
        file = phen_path / Path(filename).with_suffix(".tsv")
        rename_map = dict(zip(field_map.values(), field_map))
        try:
            df = (
                pd.read_csv(file, sep="\t")
                .pipe(_set_subsess_index)
                [list(rename_map)]
                .rename(columns=rename_map)
            )
        except KeyError as err:
            raise KeyError(
                f"Key error reading data from {filename}: {err.args[0]}"
            ) from err
        for field in field_map:
            attrs = spec[field]
            if "map" in attrs:
                df[field] = df[field].map(attrs["map"])
            if "valid" in attrs:
                df[field] = (
                    df[field].where(df[field].isin(attrs["valid"]))
                )
            if "invalid" in attrs:
                df[field] = (
                    df[field].where(~df[field].isin(attrs["invalid"]))
                )
        dfs.append(df)
    return dfs
        
    # result = copy.copy(self)
    # metadata = ft.reduce(
    #     lambda df1, df2: df1.join(df2, how="outer"), dfs, self.metadata
    # )
        
