
rule bwa_map: #
    input:
        reference=reference_path,
        forward_read=rules.extract_properly_mapped_reads.output.forward_read,
        reverse_read=rules.extract_properly_mapped_reads.output.reverse_read,
    output:
        bam=out_dir_path  / "intermediate_bam/{sample_id}/{sample_id}.intermediate.bam"
    params:
        samtools_sort_threads=parameters["threads"]["samtools_sort"],
        samtools_fixmate_threads=parameters["threads"]["samtools_fixmate"],
        samtools_sort_memory=parameters["memory_mb"]["samtools_per_thread"]
    log:
        map=output_dict["log"]  / "bwa_map.{sample_id}.map.log",
        fixmate=output_dict["log"]  / "bwa_map.{sample_id}.fixmate.log",
        sort=output_dict["log"]  / "bwa_map.{sample_id}.sort.log",
        index=output_dict["log"] / "bwa_map.{sample_id}.index.log",
        cluster_log=output_dict["cluster_log"] / "bwa_map.{sample_id}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "bwa_map.{sample_id}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "bwa_map.{sample_id}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["bwa_map"] ,
        time=parameters["time"]["bwa_map"],
        mem=parameters["memory_mb"]["bwa_map"] + parameters["memory_mb"]["samtools_per_thread"] * parameters["threads"]["samtools_sort"]
    threads: parameters["threads"]["bwa_map"] + parameters["threads"]["samtools_sort"] + parameters["threads"]["samtools_fixmate"]

    shell:
        " bwa mem -t {threads} -R  \'@RG\\tID:{wildcards.sample_id}\\tPU:x\\tSM:{wildcards.sample_id}\\tPL:Illumina\\tLB:{wildcards.sample_id}\' "
        " {input.reference} {input.forward_read} {input.reverse_read} 2>{log.map} | "
        " samtools fixmate -@ {params.samtools_fixmate_threads} -m - -  2> {log.fixmate} | "
        " samtools sort -@ {params.samtools_sort_threads} -m {params.samtools_sort_memory}m -o {output.bam} 1>{log.sort} 2>&1;"
        " samtools index {output.bam} >{log.index} 2>&1"

rule realign_bam: #
    input:
        reference=reference_path,
        bam=out_dir_path / "intermediate_bam/{sample_id}/{sample_id}.intermediate.bam",
        bed=str_loci_bed_path
    output:
        bam=out_dir_path  / "realigned_bam/{sample_id}/{sample_id}.realigned.bam"
    log:
        std=output_dict["log"]  / "realign_bam.{sample_id}.map.log",
        cluster_log=output_dict["cluster_log"] / "realign_bam.{sample_id}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "realign_bam.{sample_id}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "realign_bam.{sample_id}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["realign_bam"] ,
        time=parameters["time"]["realign_bam"],
        mem=parameters["memory_mb"]["realign_bam"]
    threads: parameters["threads"]["realign_bamp"]
    shell:
        " gatk -T IndelRealigner -targetIntervals {input.bed} -I {input.bam} -o {output.bam} "
        " -R {input.reference} > {log.std} 2>&1"
