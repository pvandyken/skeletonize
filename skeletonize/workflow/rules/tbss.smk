rule create_average_image:
    input:
        inputs["fa"].expand()
    output:
        mean=bids(source, suffix="mean.nii.gz"),
        mask=bids(source, suffix="mask.nii.gz"),
    log: log("create_average_image")
    benchmark: benchmark("create_average_image")
    group: 'tbss'
    threads: 1
    envmodules:
        "StdEnv/2020",
        "gcc/9.3.0",
        "fsl/6.0.4",
    resources:
        mem_mb=6000,
        runtime=5,
    shadow: 'minimal'
    shell:
        """
        fslmerge -t stack.nii.gz {input}
        fslmaths stack.nii.gz -max 0 -Tmin -bin mask.nii.gz -odt char
        fslmaths stack.nii.gz -mas mask.nii.gz stack.nii.gz
        fslmaths stack.nii.gz -Tmean {output.mean}
        fslmaths mask.nii.gz -Tmean {output.mask}
        """

rule skeletonize_average_image:
    input:
        rules.create_average_image.output.mean
    output:
        tempout("skeletonize_average_image", extension=".nii.gz")
    log: log("skeletonize_average_image")
    benchmark: benchmark("skeletonize_average_image")
    group: 'tbss'
    envmodules:
        "StdEnv/2020",
        "gcc/9.3.0",
        "fsl/6.0.4",
    threads: 1
    resources:
        mem_mb=1000,
        runtime=1,
    shell:
        "tbss_skeleton -i {input} -o {output}"


rule threshold_skeletonized_image:
    input:
        rules.skeletonize_average_image.output
    output:
        tempout("threshold_skeletonized_image", extension=".nii.gz")
    log: log("threshold_skeletonized_image")
    benchmark: benchmark("threshold_skeletonized_image")
    group: "tbss"
    threads: 1
    envmodules:
        "StdEnv/2020",
        "gcc/9.3.0",
        "fsl/6.0.4",
    resources:
        mem_mb=1000,
        runtime=1,
    params:
        threshold=0.2
    shell:
        "fslmaths {input} -thr {params.threshold} -bin {output}"

rule get_distance_map:
    input:
        skeleton=rules.threshold_skeletonized_image.output,
        mask=rules.create_average_image.output.mask,
    output:
        tempout("get_distance_map", extension=".nii.gz")
    log: log("get_distance_map")
    benchmark: benchmark("get_distance_map")
    group: 'tbss'
    threads: 1
    resources:
        mem_mb=1000,
        runtime=1,
    shadow: 'minimal'
    envmodules:
        "StdEnv/2020",
        "gcc/9.3.0",
        "fsl/6.0.4",
    shell:
        """
        fslmaths {input.mask} -mul -1 -add 1 -add {input.skeleton} \
            dist_mask.nii.gz
        distancemap -i dist_mask.nii.gz -o {output}
        """



rule project_onto_skeleton:
    input:
        fa=inputs["fa"].path,
        distance=rules.get_distance_map.output,
        mean=rules.create_average_image.output.mean,
        mask=rules.create_average_image.output.mask,
        **(
            {"data": inputs[comp].path} if comp != "fa" else {}
        )
    output:
        bids(
            output,
            suffix="skeletonized.nii.gz",
            **inputs[comp].wildcards,
        )
    log: log("project_onto_skeleton", inputs[comp])
    benchmark: benchmark("project_onto_skeleton", inputs[comp])
    group: 'tbss'
    threads: 1
    resources:
        mem_mb=500,
        runtime=1,
    envmodules:
        "StdEnv/2020",
        "gcc/9.3.0",
        "fsl/6.0.4",
    params:
        threshold=0.2,
        param_map=lambda wcards, input: str(input.data) if comp != "fa" else ""
    shadow: 'minimal'
    shell:
        """
        fslmaths "{input.fa}" -mas "{input.mask}" fa_masked.nii.gz
        if [[ -n "{params.param_map}" ]]; then
            fslmaths "{params.param_map}" -mas "{input.mask}" map_masked.nii.gz
            param_map="-a map_masked.nii.gz"
        else
            param_map=
        fi

        tbss_skeleton -i {input.mean} -p {params.threshold} {input.distance} \\
            ${{FSLDIR}}/data/standard/LowerCingulum_1mm fa_masked.nii.gz \\
            {output} $param_map
        """
