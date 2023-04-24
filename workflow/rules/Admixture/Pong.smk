import pandas as pd

localrules: create_files_for_pong

rule create_files_for_pong:
    priority: 1000
    input:
        q_table_list=expand(out_dir_path / "admixture/structure/{stage}/structure.K{K}.clumpp.permutted.R_{run}",
                         run=structure_run_id_list,
                         K=parameters["tool_options"]["structure"]["K_list"],
                         allow_missing=True),
        pop_tab=rules.convert_hipSTR_vcf.output.pop_tab
    output:
        filemap=out_dir_path / "admixture/structure/{stage}/pong.filemap",
        ind2pop=out_dir_path / "admixture/structure/{stage}/pong.ind2pop",
        pop_names=out_dir_path / "admixture/structure/{stage}/pong.pop_names"
    log:
        std=output_dict["log"] / "create_files_for_pong.{stage}.log",
        cluster_log=output_dict["cluster_log"] / "create_files_for_pong.cluster.{stage}.log",
        cluster_err=output_dict["cluster_error"] / "create_files_for_pong.cluster.{stage}.err",
    benchmark:
        output_dict["benchmark"] / "create_files_for_pong.benchmark.{stage}.txt"
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
        with open(output.filemap) as out_fd:
            for K in parameters["tool_options"]["structure"]["K_list"]:
                for run in structure_run_id_list:
                    out_fd.write("{0}.K{1}.R{2}\t{1}\t{3}/admixture/structure/{0}/structure.K{1}.clumpp.permutted.R_{2}\n".format(wildcards.stage, 
                                                                                                                                  K, 
                                                                                                                                  run,
                                                                                                                                  str(out_dir_path)))
        #create ind2pop and pop_names files
        pop_df = pd.read_csv(input.pop_tab, sep="tsv", header=0, )
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

rule pong:
    priority: 1000
    input:
        q_table_list=rules.clumpp.output.permutted,
        filemap=rules.create_files_for_pong.output.filemap,
        ind2pop=rules.create_files_for_pong.output.ind2pop,
        pop_names=rules.create_files_for_pong.output.pop_names
    output:
        dir=out_dir_path / "admixture/structure/{stage}/pong/"
    params:
        ignore_cols=detect_col_number_to_ignore
    log:
        std=output_dict["log"] / "pong.{stage}.log",
        cluster_log=output_dict["cluster_log"] / "pong.cluster.{stage}.log",
        cluster_err=output_dict["cluster_error"] / "pong.cluster.{stage}.err",
    benchmark:
        output_dict["benchmark"] / "pong.benchmark.{stage}.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["pong"],
        time=parameters["time"]["pong"],
        mem=parameters["memory_mb"]["pong"],
    threads:
        parameters["threads"]["pong"]
    shell:
        " pong -c {params.ignore_cols} -m {input.filemap} --ind2pop {input.ind2pop} --pop_names {input.pop_names} "
        " -o {output.dir} > {log.std} 2>&1"

