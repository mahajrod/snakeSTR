localrules: add_flanks

rule add_flanks:
    priority: 1000
    input:
        str_loci_bed_path
    output:
        out_dir_path / "str/str_loci_with_flanks.bed"
    params:
        flank_len=config["flank_len"]
    log:
        std=output_dict["log"] / "add_flanks.log",
        cut=output_dict["log"] / "add_flanks.cut.log",
        cluster_log=output_dict["cluster_log"] / "add_flanks.cluster.log",
        cluster_err=output_dict["cluster_error"] / "add_flanks.cluster.err",
    benchmark:
        output_dict["benchmark"] / "add_flanks.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["add_flanks"],
        time=parameters["time"]["add_flanks"],
        mem=parameters["memory_mb"]["add_flanks"],
    threads:
        parameters["threads"]["add_flanks"]
    shell:
         " add_flanks_to_bed.py -i {input} -l {params.flank_len} -r {params.flank_len} 2>{log.std} | "
         " cut -f 1-3 > {output} 2>{log.cut} "
