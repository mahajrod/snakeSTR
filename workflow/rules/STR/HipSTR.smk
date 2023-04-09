
rule hipSTR:
    priority: 1000
    input:
        bams=expand(out_dir_path / "realigned_bam/{sample_id}/{sample_id}.realigned.bam",
                    sample_id=sample_id_list),
        reference=reference_path,
        regions=str_repeat_bed_path
    output:
        vcf=out_dir_path / "str/hipSTR.vcf.gz",
        viz=out_dir_path / "str/hipSTR.viz.gz",
    params:
        bams=lambda wildcards: ",".join(sample_id_list),
        min_reads=parameters["tool_options"]["hipSTR"][config["parameter_set"]]["min_reads"],
    log:
        std=output_dict["log"] / "hipSTR.log",
        err=output_dict["log"] / "hipSTR.err",
        cluster_log=output_dict["cluster_log"] / "hipSTR.cluster.log",
        cluster_err=output_dict["cluster_error"] / "hipSTR.cluster.err",
    benchmark:
        output_dict["benchmark"] / "hipSTR.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["hipSTR"],
        time=parameters["time"]["hipSTR"],
        mem=parameters["memory_mb"]["hipSTR"],
    threads:
        parameters["threads"]["hipSTR"]
    shell:
         " HipSTR --bams {params.bams} --fasta {input.reference} --regions {input.regions} "
         " --min-reads {params.min_reads} --str-vcf {output.vcf} --log {log.std} --viz-out {output.viz} 2>{log.err}"
