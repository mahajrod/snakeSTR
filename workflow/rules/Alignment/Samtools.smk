
rule extract_alignments: #
    input:
        bam=input_dir_path / ("alignment/{sample_id}/{sample_id}%s.bam" % config["bam_suffix"]),
        bed=rules.add_flanks.output
    output:
        bam=out_dir_path  / "extracted_bam/{sample_id}/{sample_id}.extracted.bam"
    log:
        extract=output_dict["log"]  / "extract_alignments.{sample_id}.extract.log",
        index=output_dict["log"] / "extract_alignments.{sample_id}.index.log",
        cluster_log=output_dict["cluster_log"] / "extract_alignments.{sample_id}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "extract_alignments.{sample_id}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "extract_alignments.{sample_id}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["extract_alignments"] ,
        time=parameters["time"]["extract_alignments"],
        mem=parameters["memory_mb"]["extract_alignments"]
    threads: parameters["threads"]["extract_alignments"]

    shell:
        " samtools view -bhM -L {input.bed} -o {output.bam} {input.bam} > {log.extract} 2>&1;"
        " samtools index {output.bam} > {log.index} 2>&1"
