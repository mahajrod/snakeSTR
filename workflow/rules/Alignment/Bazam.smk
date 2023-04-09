
rule extract_properly_mapped_reads:
    input:
        bam=rules.extract_alignments.output.bam,
        bed=str_loci_bed_path
    output:
        forward=temp(out_dir_path / "extracted_fastq/{sample_id}/{sample_id}.extracted_1.fastq"),
        reverse=temp(out_dir_path / "extracted_fastq/{sample_id}/{sample_id}.extracted_2.fastq")
    log:
        std=output_dict["log"]  / "extract_properly_mapped_reads.{sample_id}.log",
        cluster_log=output_dict["cluster_log"] / "extract_properly_mapped_reads.{sample_id}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "extract_properly_mapped_reads.{sample_id}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "extract_properly_mapped_reads.{sample_id}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["extract_properly_mapped_reads"] ,
        time=parameters["time"]["extract_properly_mapped_reads"],
        mem=parameters["memory_mb"]["extract_properly_mapped_reads"]
    threads: parameters["threads"]["extract_properly_mapped_reads"]

    shell:
        " bazam -bam {input.bam} -L {input.bed} -r1 {output.forward} -r2 {output.reverse} > {log.std} 2>&1"
