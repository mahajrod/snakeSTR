localrules: generate_config_for_structure

rule generate_config_for_structure:
    priority: 1000
    input:
        loci_tab=out_dir_path / "str/hipSTR.allels.{stage}.loci.tab"
    output:
        config=out_dir_path / "admixture/structure/{stage}/structure.K{K}.R{run}.config",
        extra_config=out_dir_path / "admixture/structure/{stage}/structure.K{K}.R{run}.config.extra"
    params:
        output_prefix=lambda wildcards: out_dir_path / "admixture/structure/{0}/structure.K{1}.R{2}".format(wildcards.stage,
                                                                                                           wildcards.K,
                                                                                                           wildcards.run)
    log:
        cluster_log=output_dict["cluster_log"] / "generate_config_for_structure:.cluster.{stage}.{K}.{run}.log",
        cluster_err=output_dict["cluster_error"] / "generate_config_for_structure:.cluster.{stage}.{K}.{run}.err",

    benchmark:
        output_dict["benchmark"] / "generate_config_for_structure:.benchmark.{stage}.{K}.{run}.txt"
    #conda:
    #    config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["generate_config_for_structure"],
        time=parameters["time"]["generate_config_for_structure"],
        mem=parameters["memory_mb"]["generate_config_for_structure"],
    threads:
        parameters["threads"]["generate_config_for_structure"]
    run:
        with open(input.loci_tab, "r") as in_fd:
            number_of_loci = len(in_fd.readline().split("\t"))
            number_of_individuals = 0
            for line in in_fd:
                number_of_individuals += 1
            number_of_individuals //= 2

        Path(output.config).parent.mkdir(parents=True, exist_ok=True)

        with open(output.config, "w") as out_fd:
            out_fd.write("#define OUTFILE {0}\n".format(params.output_prefix))
            out_fd.write("#define INFILE {0}\n".format(input.loci_tab))
            out_fd.write("#define MISSING {0}\n".format(config["absent_allel_alias"]))
            out_fd.write("#define PLOIDY {0}\n".format(config["ploidy"]))
            out_fd.write("#define MAXPOPS {0}\n".format(wildcards.K))
            out_fd.write("#define NUMINDS {0}\n".format(number_of_individuals))
            out_fd.write("#define NUMLOCI {0}\n".format(number_of_loci))

            for parameter in parameters["tool_options"]["structure"]["config_file_parameters"]:
                out_fd.write("#define {0} {1}\n".format(parameter,
                                                        parameters["tool_options"]["structure"]["config_file_parameters"][parameter]))
        with open(output.extra_config,"w") as out_fd:
            for extra_parameter in parameters["tool_options"]["structure"]["extra_config_file_parameters"]:
                out_fd.write("#define {0} {1}\n".format(extra_parameter,
                                                        parameters["tool_options"]["structure"]["config_file_parameters"][extra_parameter]))

rule structure:
    priority: 1000
    input:
        loci_tab=out_dir_path / "str/hipSTR.allels.{stage}.loci.tab",
        config=rules.generate_config_for_structure.output.config,
        extra_config=rules.generate_config_for_structure.output.extra_config
    output:
        res=out_dir_path / "admixture/structure/{stage}/structure.K{K}.R{run}_f",

    params:
        output_prefix=lambda wildcards: out_dir_path / "admixture/structure/{0}/structure.K{1}.R{2}".format(wildcards.stage,
                                                                                                           wildcards.K,
                                                                                                           wildcards.run),
        #extraparams_file=parameters["tool_options"]["structure"]["extraparams_file"]
    log:
        std=output_dict["log"] / "structure.{stage}.{K}.{run}.log",
        err=output_dict["log"] / "structure.{stage}.{K}.{run}.err",
        cluster_log=output_dict["cluster_log"] / "structure.cluster.{stage}.{K}.{run}.log",
        cluster_err=output_dict["cluster_error"] / "structure.cluster.{stage}.{K}.{run}.err",
    benchmark:
        output_dict["benchmark"] / "structure.benchmark.{stage}.{K}.{run}.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["structure"],
        time=parameters["time"]["structure"],
        mem=parameters["memory_mb"]["structure"],
    threads:
        parameters["threads"]["structure"]
    shell:
         #" NUM_LINES=`wc -l {input.loci_tab} | cut -f 1 -d ' '`; "
         #" let NUM_INDIV=(${{NUM_LINES}}-1)/2;"
         #" NUM_LOCI=`head -n 1 {input.loci_tab} | awk -F'\t' '{{print NF}}'`; "
         " structure -m {input.config} -e {input.extra_config} > {log.std} 2>{log.err}" # -K {wildcards.K} -L ${{NUM_LOCI}} -N ${{NUM_INDIV}}  -i {input.loci_tab} -o {params.output_prefix}
