
rule hipSTR:
    priority: 1000
    input:
        bams=expand(out_dir_path / "realigned_bam/{sample_id}/{sample_id}.realigned.bam",
                    sample_id=sample_id_list),
        reference=reference_path,
        regions=str_repeat_bed_path
    output:
        vcf=out_dir_path / "str/hipSTR.raw.vcf.gz",
        viz=out_dir_path / "str/hipSTR.raw.viz.gz",
    params:
        bams=lambda wildcards: ",".join(expand(out_dir_path / "realigned_bam/{sample_id}/{sample_id}.realigned.bam",
                                               sample_id=sample_id_list)),
        min_reads=parameters["tool_options"]["hipSTR"]["min_reads"],
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

rule filter_hipSTR:
    priority: 1000
    input:
        vcf=rules.hipSTR.output.vcf
    output:
        vcf=out_dir_path / "str/hipSTR.filtered.vcf.gz",
    params:
        min_call_qual=parameters["tool_options"]["hipSTR"]["min_call_qual"],
        max_call_flank_indel=parameters["tool_options"]["hipSTR"]["max_call_flank_indel"],
        max_call_stutter=parameters["tool_options"]["hipSTR"]["max_call_stutter"],
        min_call_allele_bias=parameters["tool_options"]["hipSTR"]["min_call_allele_bias"],
        min_call_strand_bias=parameters["tool_options"]["hipSTR"]["min_call_strand_bias"],
    log:
        std=output_dict["log"] / "filter_hipSTR.log",
        cluster_log=output_dict["cluster_log"] / "filter_hipSTR.cluster.log",
        cluster_err=output_dict["cluster_error"] / "filter_hipSTR.vluster.err",
    benchmark:
        output_dict["benchmark"] / "filter_hipSTR.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["filter_hipSTR"],
        time=parameters["time"]["filter_hipSTR"],
        mem=parameters["memory_mb"]["filter_hipSTR"],
    threads:
        parameters["threads"]["filter_hipSTR"]
    shell:
         " filter_vcf.py  --vcf {input.vcf} "
         " --min-call-qual {params.min_call_qual} "
         " --max-call-flank-indel {params.max_call_flank_indel} "
         " --max-call-stutter {params.max_call_stutter}  "
         " --min-call-allele-bias {params.min_call_allele_bias} "
         " --min-call-strand-bias {params.min_call_strand_bias} | gzip -c > {output.vcf} 2>{log.std}"

def postprocess_allel_file_func(wildcards):
    if config["stage_coretools"]["admixture"] == "structure":
        if parameters["tool_options"]["structure"]["config_file_parameters"]["POPDATA"] == 1:
            return " | sed '1s/^[^\t]*\t[^\t]*\t//' " # remove labels of first and second columns from header for compatibility with STRUCTURE if POPDATA is 1
        elif parameters["tool_options"]["structure"]["config_file_parameters"]["POPDATA"] == 0:
            return " | sed '1s/^[^\t]*\t//' " # remove first column label from header for compatibility with STRUCTURE if POPDATA is 0
        else:
            raise ValueError("ERROR!!! Unrecognized value of POPDATA parameter (STRUCTURE config)!")
    else:
        return ""


rule convert_hipSTR_vcf:
    priority: 1000
    input:
        vcf=out_dir_path / "str/hipSTR.{stage}.vcf.gz"
    output:
        str_tab=out_dir_path / "str/hipSTR.allels.{stage}.str.tab",
        loci_tab=out_dir_path / "str/hipSTR.allels.{stage}.loci.tab",
        pop_tab=out_dir_path / "str/hipSTR.allels.{stage}.pop.tab"
    params:
        absent_allel_alias=config["absent_allel_alias"],
        len_file=" -l {0} --amplicon_id_column_name {1} --amplicon_len_column_name {2} ".format(config["str_loci_len_file"],
                                                                                                parameters["tool_options"]["convert_hipSTR_vcf"]["amplicon_id_column_name"],
                                                                                                parameters["tool_options"]["convert_hipSTR_vcf"]["amplicon_len_column_name"],) if config["str_loci_len_file"] else "",
        postprocessing=postprocess_allel_file_func,
        add_population_column=" --add_population_column " if parameters["tool_options"]["structure"]["config_file_parameters"]["POPDATA"] == 1 else "",
        pop_file=" --pop_file {0} ".format(config["pop_file"]) if config["pop_file"] else ""
    log:
        str=output_dict["log"] / "convert_hipSTR_vcf.{stage}.str.log",
        loci=output_dict["log"] / "convert_hipSTR_vcf.{stage}.loci.log",
        cluster_log=output_dict["cluster_log"] / "convert_hipSTR_vcf.{stage}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "convert_hipSTR_vcf.{stage}.cluster.err",
    benchmark:
        output_dict["benchmark"] / "convert_hipSTR_vcf.{stage}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["convert_hipSTR_vcf"],
        time=parameters["time"]["convert_hipSTR_vcf"],
        mem=parameters["memory_mb"]["convert_hipSTR_vcf"],
    threads:
        parameters["threads"]["convert_hipSTR_vcf"]
    shell:
         " convert_STR_vcf_to_allel_length.py {params.add_population_column} {params.pop_file} --encode_ids "
         " --pop_df_file {output.pop_tab} "
         " -i {input.vcf} {params.len_file} 2>{log.loci} {params.postprocessing} > {output.loci_tab};"
         " convert_STR_vcf_to_allel_length.py {params.add_population_column} {params.pop_file} --encode_ids "
         " -i {input.vcf}  2>{log.str} {params.postprocessing} > {output.str_tab};"
