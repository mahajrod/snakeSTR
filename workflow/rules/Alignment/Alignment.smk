
rule bwa_map: #
    input:
        index=rules.bwa_index.output.index,
        reference=out_dir_path  / ("{assembly_stage}/{assembler}/%s.{assembly_stage}.{assembler}.{haplotype}.fasta" % config["genome_name"]),
        fastq=output_dict["data"] / ("fastq/hic/raw/{fileprefix}%s" % config["fastq_extension"])
    output:
        bam=out_dir_path  / ("{assembly_stage}/{assembler}/{haplotype}/alignment/%s.{assembly_stage}.{assembler}.{haplotype}.bwa.filtered.{fileprefix}.bam" % config["genome_name"]) # TODO: make_temp
    params:
        id="{0}_hic".format(config["genome_name"])
    log:
        map=output_dict["log"]  / "bwa_map.{assembler}.{assembly_stage}.{haplotype}.{fileprefix}.map.log",
        sort=output_dict["log"]  / "bwa_map.{assembler}.{assembly_stage}.{haplotype}.{fileprefix}.sort.log",
        filter=output_dict["log"]  / "bwa_map.{assembler}.{assembly_stage}.{haplotype}.{fileprefix}.filter.log",
        cluster_log=output_dict["cluster_log"] / "bwa_map.{assembler}.{assembly_stage}.{haplotype}.{fileprefix}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "bwa_map.{assembler}.{assembly_stage}.{haplotype}.{fileprefix}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "bwa_map.{assembler}.{assembly_stage}.{haplotype}.{fileprefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["bwa_map"] ,
        time=parameters["time"]["bwa_map"],
        mem=parameters["memory_mb"]["bwa_map"]
    threads: parameters["threads"]["bwa_map"]

    shell:
        " bwa mem -t {threads} -R  \'@RG\\tID:{params.id}\\tPU:x\\tSM:{params.id}\\tPL:Illumina\\tLB:x\' "
        " {input.reference} {input.fastq} 2>{log.map} | filter_five_end.pl 2>{log.filter} | samtools view -Sb - > {output.bam} 2>{log.sort} "


rule bam_merge_pairs:
    input:
        forward_bam=out_dir_path  / ("{assembly_stage}/{assembler}/{haplotype}/alignment/%s.{assembly_stage}.{assembler}.{haplotype}.bwa.filtered.{pairprefix}%s.bam" % (config["genome_name"], input_forward_suffix_dict["hic"])),
        reverse_bam=out_dir_path  / ("{assembly_stage}/{assembler}/{haplotype}/alignment/%s.{assembly_stage}.{assembler}.{haplotype}.bwa.filtered.{pairprefix}%s.bam" % (config["genome_name"], input_reverse_suffix_dict["hic"])),
        reference_fai=rules.ref_faidx.output.fai
    output:
        bam=out_dir_path  / ("{assembly_stage}/{assembler}/{haplotype}/alignment/%s.{assembly_stage}.{assembler}.{haplotype}.bwa.filtered.{pairprefix}.bam" % config["genome_name"]) # TODO: make_temp
    params:
        min_mapq=parameters["tool_options"]["two_read_bam_combiner"]["mapq"],
        sort_threads=parameters["threads"]["samtools_sort"],
        sort_memory=parameters["memory_mb"]["samtools_sort"],
        tmp_dir=lambda wildcards: out_dir_path  / "hic_scaffolding/{0}/{1}/alignment/{2}".format(wildcards.assembler, wildcards.haplotype, wildcards.pairprefix)
    log:
        merge=output_dict["log"] / "bam_merge_pairs.{assembler}.{assembly_stage}.{haplotype}.{pairprefix}.merge.log",
        view=output_dict["log"] / "bam_merge_pairs.{assembler}.{assembly_stage}.{haplotype}.{pairprefix}.view.log",
        sort=output_dict["log"] / "bam_merge_pairs.{assembler}.{assembly_stage}.{haplotype}.{pairprefix}.sort.log",
        cluster_log=output_dict["cluster_log"] / "bam_merge_pairs.{assembler}.{assembly_stage}.{haplotype}.{pairprefix}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "bam_merge_pairs.{assembler}.{assembly_stage}.{haplotype}.{pairprefix}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "bam_merge_pairs.{assembler}.{assembly_stage}.{haplotype}.{pairprefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["two_read_bam_combiner"] ,
        time=parameters["time"]["two_read_bam_combiner"],
        mem=parameters["memory_mb"]["two_read_bam_combiner"] + parameters["memory_mb"]["samtools_sort"]
    threads: parameters["threads"]["two_read_bam_combiner"] + parameters["threads"]["samtools_sort"]

    shell:
        " two_read_bam_combiner.pl {input.forward_bam} {input.reverse_bam} samtools {params.min_mapq} 2>{log.merge} | "
        " samtools view -bS -t {input.reference_fai} - 2>{log.view} | "
        " samtools sort -m {params.sort_memory}M -@ {params.sort_threads} -o {output.bam} 2>{log.sort}"

rule bam_merge_files:
    input:
        bams=expand(out_dir_path  / ("{assembly_stage}/{assembler}/{haplotype}/alignment/%s.{assembly_stage}.{assembler}.{haplotype}.bwa.filtered.{pairprefix}.bam" % config["genome_name"]),
                    allow_missing=True,
                    pairprefix=input_pairprefix_dict["hic"]),
        reference_fai=rules.ref_faidx.output.fai,
        reference=out_dir_path  / ("{assembly_stage}/{assembler}/%s.{assembly_stage}.{assembler}.{haplotype}.fasta" % config["genome_name"])
    output:
        bam=out_dir_path  / ("{assembly_stage}/{assembler}/{haplotype}/alignment/%s.{assembly_stage}.{assembler}.{haplotype}.bwa.filtered.bam"  % config["genome_name"]) # TODO: make temp
    params:
        sort_threads=parameters["threads"]["samtools_sort"]
    log:
        std=output_dict["log"] / "bam_merge_files.{assembler}.{assembly_stage}.{haplotype}.log",
        cluster_log=output_dict["cluster_log"] / "bam_merge_files.{assembler}.{assembly_stage}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "bam_merge_files.{assembler}.{assembly_stage}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "bam_merge_files.{assembler}.{assembly_stage}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["samtools_sort"] ,
        time=parameters["time"]["samtools_sort"],
        mem=parameters["memory_mb"]["samtools_sort"]
    threads: parameters["threads"]["samtools_sort"]

    shell:
        " samtools merge -@ {params.sort_threads} --reference  -o {output.bam} {input.bams} 1>{log.std} 2>&1"

rule rmdup:
    input:
        bam=rules.bam_merge_files.output.bam
    output:
        bam=out_dir_path  / ("{assembly_stage}/{assembler}/{haplotype}/alignment/%s.{assembly_stage}.{assembler}.{haplotype}.bwa.filtered.rmdup.bam"  % config["genome_name"]),
        bai=out_dir_path  / ("{assembly_stage}/{assembler}/{haplotype}/alignment/%s.{assembly_stage}.{assembler}.{haplotype}.bwa.filtered.rmdup.bam.bai"  % config["genome_name"]),
        dup_stats=out_dir_path  / ("{assembly_stage}/{assembler}/{haplotype}/alignment/%s.{assembly_stage}.{assembler}.{haplotype}.bwa.filtered.rmdup.stats"  % config["genome_name"]),
        bam_stats=out_dir_path  / ("{assembly_stage}/{assembler}/{haplotype}/alignment/%s.{assembly_stage}.{assembler}.{haplotype}.bwa.filtered.rmdup.bam.stats"  % config["genome_name"])
    params:
        tmp_dir=lambda wildcards: out_dir_path  / "hic_scaffolding/{0}/{1}/alignment/temp".format(wildcards.assembler, wildcards.haplotype)
    log:
        std=output_dict["log"] / "rmdup.{assembler}.{assembly_stage}.{haplotype}.log",
        cluster_log=output_dict["cluster_log"] / "rmdup.{assembler}.{assembly_stage}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "rmdup.{assembler}.{assembly_stage}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "rmdup.{assembler}.{assembly_stage}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["rmdup"] ,
        time=parameters["time"]["rmdup"],
        mem=parameters["memory_mb"]["rmdup"]
    threads: parameters["threads"]["rmdup"]

    shell:
        " picard -Xmx{resources.mem}m MarkDuplicates -I {input} -O {output.bam} " 
        " --REMOVE_DUPLICATES true -M {output.dup_stats}  --TMP_DIR {params.tmp_dir}1>{log.std} 2>&1;"
        " samtools index {output.bam}; "
        " get_stats.pl {output.bam} > {output.bam_stats}"

rule bam2bed:
    input:
        bam=rules.rmdup.output.bam
    output:
        bed=out_dir_path  / ("{assembly_stage}/{assembler}/{haplotype}/alignment/%s.{assembly_stage}.{assembler}.{haplotype}.bwa.filtered.rmdup.bed"  % config["genome_name"])
    log:
        convert=output_dict["log"] / "bam2bed.{assembler}.{assembly_stage}.{haplotype}.convert.log",
        sort=output_dict["log"] / "bam2bed.{assembler}.{assembly_stage}.{haplotype}.sort.log",
        cluster_log=output_dict["cluster_log"] / "bam2bed.{assembler}.{assembly_stage}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "bam2bed.{assembler}.{assembly_stage}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "bam2bed.{assembler}.{assembly_stage}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["bam2bed"] ,
        time=parameters["time"]["bam2bed"],
        mem=parameters["memory_mb"]["bam2bed"]
    threads: parameters["threads"]["bam2bed"]

    shell:
        " bamToBed -i {input.bam} 2>{log.convert} | sort -k 4 > {output.bed} 2>{log.sort}"
