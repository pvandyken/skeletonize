def get_entity_filters(filters):
    return {
        key: val for key, val in filters.items()
        if key not in {"scope", "extension"} and "_" not in key
    }

def log_transform(*args):
    def _log_transform(x):
        return np.log10(x)

    return preproc.FunctionTransformer(_log_transform)


def impute_transform(df, field, strategy, fill=None):
    return impute.SimpleImputer(
        strategy=strategy,
        fill_value=(
            df[field].dtype.type(fill) if fill is not None else None
        )
    )

def _squeeze_text(s):
    u = s.unique()
    if len(u) > 1:
        raise ValueError(f"Differing text values when grouping:\n{s}")
    return u.sum()

def group_sessions(df):
    aggs = {col: _squeeze_text if df[col].dtype == "O" else "mean" for col in df}
    return df.groupby("participant_id").agg(aggs)


@ft.lru_cache(None)
def _get_participants():
    if config.get("metadata_mapping"):
        with open(config["metadata_mapping"]) as f:
            phenotypes = yaml.safe_load(f)
    else:
        phenotypes = None
    df = load_metadata(
        Path(config['bids_dir'], 'participants.tsv'),
        Path(config['bids_dir'], 'phenotypes'),
        phenotypes
    )
    def squeeze(item):
        if len(item) == 1:
            return item[0]
        return item

    df = df.loc[
        [
            squeeze(e) for e in zip(
                *inputs[comp][tuple(inputs.subj_wildcards)].zip_lists.values()
            )
            if squeeze(e) in df.index
        ]
    ]

    if len(inputs.subj_wildcards) == 1 and len(df.index.names) == 2:
        df = group_sessions(df)

    # if "session" in inputs[comp].wildcards and set(inputs[comp].entities["session"]) == {"1", "2"}:
    #     df = df.loc[[
    #         s for s in inputs[comp].entities["subject"]
    #         if len(
    #             set(inputs[comp].filter(subject=s)["session"]) & {"1", "2"}
    #         ) == 2
    #     ]]

    filters = {
        filt[0]: filt[1:] for filt in config.get("filter_participants", []) or []
    }
    for col, values in filters.items():
        if col not in df:
            continue
        df = df[df[col].isin(map(df[col].dtype.type, itx.always_iterable(values)))]

    if config["transform"]:
        transforms = {
            "log": log_transform,
            "impute": impute_transform
        }
        for col, *txfs in config["transform"]:
            txfs_parsed = [shlex.split(txf) for txf in txfs]
            pp = pipeline.make_pipeline(
                *(
                    transforms[txf](df, col, *txf_args)
                    for txf, *txf_args in txfs_parsed
                )
            )
            transformed = compose.make_column_transformer(
                (pp, [col])
            ).fit_transform(df)
            df[col] = transformed[:, 0]

    df = _handle_nan(df)


    # for term in terms:
    #     if np.any(pd.isnull(df.get(term, []))):
    #         raise ConfigError(
    #             f"Found nan in column '{term}'. NaNs must either be explicitly "
    #             "filtered using --skip-nan, or an imputation strategy must be used"
    #         )
    return df.sort_index()

def _handle_nan(df):
    indices = df.index.names
    na_action = 'drop' if config['skip_nan'] else 'raise'
    return df.loc[
        dmatrix(
            config["design"],
            df.reset_index().set_index(indices, drop=False),
            NA_action=na_action,
            return_type="dataframe"
        ).index
    ]

@ft.lru_cache(None)
def _unpack_contrasts():
    contrasts, *_params = it.zip_longest(*config["ttest"], fillvalue=None)
    params = [
        dict(
            param.split("=", 1) for param in row if param is not None
        )
        for row in zip(*_params)
    ]
    return contrasts, params or [{} for _ in range(len(contrasts))]

def _get_terms(df):
    return dmatrix(config["design"], df).design_info.term_names

def expand_subjects_from_group(path, allow_missing=False, **kwargs):
    participants = _get_participants()
    def inner(wcards):
        filtered = (
            pd.DataFrame(inputs[comp].zip_lists)
            .set_index(list(inputs.subj_wildcards.keys()))
            .loc[participants.index]
            .sort_index()
            .reset_index()
            .rename(columns={"participant_id": "subject"})
        )
        # _, counts = np.unique(
        #     filtered.index.get_level_values("subject"), return_counts=True
        # )
        # if len(np.unique(counts)) > 1:
        #     raise ConfigError(
        #         "Not all subjects have the same number of data files. Use filters to "
        #         "ensure each subject appears the same number of times."
        #     )
        if _get_design_matrix().shape[0] != filtered.shape[0]:
            raise ConfigError(
                "Number of data files different than number of rows in design matrix. "
                "You may need to add additional filters"
            )
        return expand(
            expand(
                path,
                zip,
                allow_missing=True,
                **filtered,
            ),
            allow_missing=allow_missing,
            **kwargs,
            **wcards,
        )
    return inner


def _get_design_matrix():
    df = _get_participants()
    na_action = 'drop' if config['skip_nan'] else 'raise'
    return dmatrix(config["design"], df.reset_index(), NA_action=na_action)


def _get_contrast_matrix(dm):
    contrasts, _ = _unpack_contrasts()
    contrast = ", ".join(contrasts)
    return dm.design_info.linear_constraint(contrast).coefs

def _get_exchange_blocks():
    if (exchange_field := config.get("exchange_blocks")) is None:
        return None
    participants = _get_participants()
    _, indices = np.unique(
        participants.reset_index()[exchange_field],
        return_inverse=True
    )
    _, validate_counts = np.unique(indices, return_counts=True)
    if len(np.unique(validate_counts)) > 1:
        msg = ["Each exchange block must have the same number of entries"]
        for size, group in (
            participants
            .reset_index()
            .groupby(exchange_field)
            .count()
            .max(axis=1)
            .to_frame(0).groupby(0)[0]
        ):
            msg.append(f"Values with {size}: " + ", ".join(group.index))
        # raise ConfigError("\n".join(msg))
    return indices


class Counter:
    def __init__(self):
        self.data = {
            "": 0
        }
    
    def next(self, item):
        if item not in self.data:
            self.data[item] = 0
            return item
        self.data[item] += 1
        return f"{item}{self.data[item]}"

def _get_contrast_labels(dm):
    labels = []
    counter = Counter()
    for contrast, params in zip(*_unpack_contrasts()):
        try:
            num_contrasts = dm.design_info.linear_constraint(contrast).coefs.shape[0]
        except:
            print(dm.design_info.column_names)
            raise
        name = params.get("name", "")
        for _ in range(num_contrasts):
            labels.append(counter.next(name))
    return labels


rule test_patsy_dm:
    input:
        Path(config['bids_dir'], 'participants.tsv'),
    run:
        dm = _get_design_matrix()
        part = _get_participants()
        dm_df = pd.DataFrame(dm, columns=dm.design_info.column_names)
        dm_df.index = part.index

        print(dm_df.to_csv(sep="\t"), end='')
        # print(dm.design_info.column_names)
        # print(dm)
        # con = _get_contrast_matrix(dm)
        # print("Contrast")
        # print(repr(con))
        # if config.get("exchange_blocks") is not None:
        #     exch = _get_exchange_blocks()
        #     print("Exchange Blocks")
        #     print(exch)

rule make_patsy_dm:
    input:
        Path(config['bids_dir'], 'participants.tsv'),
    output:
        mat=tempout("make_patsy_design_matrix", extension=".mat.txt"),
        con=tempout("make_patsy_design_matrix", extension=".con.txt"),
    run:
        dm = _get_design_matrix()
        print(f"saving to {output.mat}")
        np.savetxt(output.mat, dm, fmt="%i")
        np.savetxt(output.con, _get_contrast_matrix(dm), fmt="%i")


rule make_exchange_spec:
    input:
        Path(config['bids_dir'], 'participants.tsv'),
    output:
        tempout("make_exchange_spec", extension=".xch.txt"),
    run:
        exch = _get_exchange_blocks()
        np.savetxt(output[0], exch, fmt="%i")


rule split_contrast_matrix:
    input:
        rules.make_patsy_dm.output.con,
    output:
        tempout("split_contrast_matrix", "{contrast_label}", extension=".con.txt"),
    params:
        contrast_idx=lambda wcards: (
            _get_contrast_labels(_get_design_matrix()).index(wcards["contrast_label"])
            + 1
        )
    shell:
        "sed '{params.contrast_idx}q;d' {input} > {output}"


rule convert_design_matrices_to_fsl:
    input:
        mat=rules.make_patsy_dm.output["mat"],
        con=rules.split_contrast_matrix.output[0]
    output:
        mat=tempout("convert_design_matrices_to_fsl", "{contrast_label}", extension=".mat"),
        con=tempout("convert_design_matrices_to_fsl", "{contrast_label}", extension=".con"),
    envmodules:
        "StdEnv/2020",
        "gcc/9.3.0",
        "fsl/6.0.4"
        
    threads: 1
    shell:
        """
        Text2Vest {input.mat} {output.mat}
        Text2Vest {input.con} {output.con}
        """

rule convert_exchange_matrix_to_fsl:
    input:
        rules.make_exchange_spec.output[0],
    output:
        tempout("convert_exchange_matrix_to_fsl", "{contrast_label}", extension=".grp"),
    envmodules:
        "StdEnv/2020",
        "gcc/9.3.0",
        "fsl/6.0.4"
    shell:
        """
        Text2Vest {input} {output}
        """
