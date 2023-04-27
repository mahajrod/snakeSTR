import pandas as pd

localrules: create_files_for_pong, generate_pong_script

rule create_files_for_pong:
    priority: 1000
    input:
        q_table_list=expand(out_dir_path / "admixture/structure/{stage}/{loci_subset}/structure.{loci_subset}.K{K}.clumpp.permutted.R_{run}",
                         run=structure_run_id_list,
                         K=parameters["tool_options"]["structure"]["K_list"],
                         allow_missing=True),
        pop_tab=rules.convert_hipSTR_vcf.output.pop_tab
    output:
        filemap=out_dir_path / "admixture/structure/{stage}/{loci_subset}/pong.filemap",
        ind2pop=out_dir_path / "admixture/structure/{stage}/{loci_subset}/pong.ind2pop",
        pop_names=out_dir_path / "admixture/structure/{stage}/{loci_subset}/pong.pop_names"
    log:
        std=output_dict["log"] / "create_files_for_pong.{stage}.{loci_subset}.log",
        cluster_log=output_dict["cluster_log"] / "create_files_for_pong.cluster.{stage}.{loci_subset}.log",
        cluster_err=output_dict["cluster_error"] / "create_files_for_pong.cluster.{stage}.{loci_subset}.err",
    benchmark:
        output_dict["benchmark"] / "create_files_for_pong.benchmark.{stage}.{loci_subset}.txt"
    #conda:
    #    config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["create_files_for_pong"],
        time=parameters["time"]["create_files_for_pong"],
        mem=parameters["memory_mb"]["create_files_for_pong"],
    threads:
        parameters["threads"]["create_files_for_pong"]
    run:
        #create filemap
        with open(output.filemap, "w") as out_fd:
            for K in parameters["tool_options"]["structure"]["K_list"]:
                for run in structure_run_id_list:
                    #out_fd.write("{0}_K{1}_R{2}\t{1}\t{3}/admixture/structure/{0}/structure.K{1}.clumpp.permutted.R_{2}\n".format(wildcards.stage,
                    #                                                                                                              K,
                    #                                                                                                              run,
                    #                                                                                                              str(out_dir_path.absolute())))
                    out_fd.write("{0}_K{1}_R{2}\t{1}\tstructure.{3}.K{1}.clumpp.permutted.R_{2}\n".format(wildcards.stage,
                                                                                                          K,
                                                                                                          run,
                                                                                                          wildcards.loci_subset
                                                                                                          ))

        #create ind2pop and pop_names files
        pop_df = pd.read_csv(input.pop_tab, sep="\t", header=0, )
        #print(pop_df)
        if config["pop_file"]:
            pop_df[["pop_id"]].to_csv(output.ind2pop, sep="\t", header=False,index=False)
            pop_df[["pop_id", "pop_id"]].to_csv(output.pop_names,sep="\t",header=False,index=False)
        else:
            pop_df[["id"]].to_csv(output.ind2pop,sep="\t",header=False,index=False)
            pop_df[["id", "id"]].to_csv(output.pop_names,sep="\t",header=False,index=False)

def detect_col_number_to_ignore(wildcards):
    if config["stage_coretools"]["admixture"] == "structure":
        return 5
    elif config["stage_coretools"]["admixture"] == "admixture":
        return 0

rule generate_pong_script:
    priority: 1000
    input:
        filemap=rules.create_files_for_pong.output.filemap,
        ind2pop=rules.create_files_for_pong.output.ind2pop,
        pop_names=rules.create_files_for_pong.output.pop_names
    output:
        #dir=out_dir_path / "admixture/structure/{stage}/pong/",
        pong_script=out_dir_path / "admixture/structure/{stage}/{loci_subset}/pong.sh"
    params:
        ignore_cols=detect_col_number_to_ignore,
        out_dir=lambda wildcards: out_dir_path / "admixture/structure/{0}/{1}pong/".format(wildcards.stage, wildcards.loci_subset)
    log:
        std=output_dict["log"] / "generate_pong_script.{stage}.{loci_subset}.log",
        cluster_log=output_dict["cluster_log"] / "generate_pong_script.cluster.{stage}.{loci_subset}.log",
        cluster_err=output_dict["cluster_error"] / "generate_pong_script.cluster.{stage}.{loci_subset}.err",
    benchmark:
        output_dict["benchmark"] / "generate_pong_script.benchmark.{stage}.{loci_subset}.txt"
    #conda:
    #    config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["generate_pong_script"],
        time=parameters["time"]["generate_pong_script"],
        mem=parameters["memory_mb"]["generate_pong_script"],
    threads:
        parameters["threads"]["generate_pong_script"]
    run:
        with open(output.pong_script, "w") as out_fd:
            out_fd.write("#!/usr/bin/env bash\n")
            out_fd.write("pong -c {0} -m {1} --ind2pop {2} --pop_names {3} -o {4} \n".format(params.ignore_cols,
                                                                                             Path(input.filemap).name,
                                                                                             Path(input.ind2pop).name,
                                                                                             Path(input.pop_names).name,
                                                                                             Path(params.out_dir).name))
        shell("chmod +x {output.pong_script}")

