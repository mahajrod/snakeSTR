
rule pretextmap: #
    input:
        bam=out_dir_path  / ("{assembly_stage}/{assembler}/{haplotype}/alignment/%s.{assembly_stage}.{assembler}.{haplotype}.bwa.filtered.rmdup.bam"  % config["genome_name"]),
    output:
        map=out_dir_path  / ("{assembly_stage}/{assembler}/{haplotype}/alignment/%s.{assembly_stage}.{assembler}.{haplotype}.bwa.filtered.rmdup.map.pretext"  % config["genome_name"]),
    params:
        min_mapq=parameters["tool_options"]["pretextmap"]["mapq"],
        sortby=parameters["tool_options"]["pretextmap"]["sortby"],
        sortorder=parameters["tool_options"]["pretextmap"]["sortorder"]
    log:
        view=output_dict["log"]  / "pretextmap.{assembler}.{assembly_stage}.{haplotype}.view.log",
        map=output_dict["log"]  / "pretextmap.{assembler}.{assembly_stage}.{haplotype}.map.log",
        cluster_log=output_dict["cluster_log"] / "pretextmap.{assembler}.{assembly_stage}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "pretextmap.{assembler}.{assembly_stage}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "pretextmap.{assembler}.{assembly_stage}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["pretextmap"] ,
        time=parameters["time"]["pretextmap"],
        mem=parameters["memory_mb"]["pretextmap"]
    threads: parameters["threads"]["pretextmap"]

    shell:
        " samtools view -h {input.bam} 2>{log.view} | "
        " PretextMap -o {output.map} --sortby {params.sortby} --sortorder {params.sortorder} "
        " --mapq {params.min_mapq} > {log.map} 2>&1"

