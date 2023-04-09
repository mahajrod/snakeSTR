import os
#import logging
import shutil
from copy import deepcopy
from collections import OrderedDict
from collections.abc import Mapping
from copy import deepcopy
from pathlib import Path, PosixPath

import yaml

import pandas as pd

#logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%d/%m/%Y %I:%M:%S %p')

#---- Functions ----
def p_distance(seq_a, seq_b, seq_len):
    dist = 0
    for i in range(0, seq_len):
        if seq_a[i] != seq_b[i]:
            dist += 1
    return dist

def get_common_prefix_ans_suffixes(seq_a, seq_b):
    seq_a_len = len(seq_a)
    seq_b_len = len(seq_b)

    prefix = ""
    for i in range(0, min(seq_a_len, seq_b_len)):
        if seq_a[i] != seq_b[i]:
           return prefix, seq_a[i:], seq_b[i:]
        prefix += seq_a[i]
    return prefix, "", ""

def convert_posixpath2str_in_dict(dictionary):
    output_dictionary = deepcopy(dictionary)
    for entry in output_dictionary:
        if isinstance(output_dictionary[entry], PosixPath):
            output_dictionary[entry] = str(output_dictionary[entry])
        else:
            if not isinstance(output_dictionary[entry], Mapping): # check if existing entry is not dictionary or dictionary like
                continue # exit from recursion
            output_dictionary[entry] = convert_posixpath2str_in_dict(output_dictionary[entry])

    return output_dictionary


def copy_absent_entries(input_dictionary, output_dictionary):
    for entry in input_dictionary:
        if entry not in output_dictionary:
            output_dictionary[entry] = deepcopy(input_dictionary[entry])
        else:
            if not isinstance(output_dictionary[entry], Mapping): # check if existing entry is not dictionary or dictionary like
                continue # exit from recursion
            copy_absent_entries(input_dictionary[entry], output_dictionary[entry])

#----


#-- Initialization of path variables from config file --
#logging.info("Initialization of path variables...")
#---- Initialization of path variables for input----
input_dir_path = Path(config["input_dir"])

input_dict = {}

#---- Initialization of path variables for output ----
out_dir_path = Path(config["out_dir"])
output_dict = {}

for first_level_sub_dir in config["first_level_subdir_list"]:
    output_dict[first_level_sub_dir] = out_dir_path / config["{0}_subdir".format(first_level_sub_dir)]


#--------

#----

#---- Checking input files ----
#logging.info("Checking input files...")

input_filedict = {}
input_file_prefix_dict = {}
input_forward_suffix_dict = {}
input_reverse_suffix_dict = {}
input_pairprefix_dict = {}

for d_type in fastq_based_data_type_set:
    input_filedict[d_type] = find_fastqs(input_dict[d_type]["fastq_dir"], fastq_extension=config["fastq_extension"])
    input_file_prefix_dict[d_type] = list(map(lambda s: str(s.name)[:-len(config["fastq_extension"])],
                                                input_filedict[d_type]))
# check filenames of paired data
for d_type in set(config["paired_fastq_based_data"]) & fastq_based_data_type_set:
   if (len(input_filedict[d_type]) % 2) != 0:
        raise ValueError("ERROR!!! {0} fastq files seems to be unpaired or misrecognized".format(d_type))
   for forward, reverse in zip(input_filedict[d_type][::2], input_filedict[d_type][1::2]):
        #print(forward, reverse)
        if p_distance(str(forward), str(reverse), len(str(forward))) > 1:
            raise ValueError("ERROR!!! Forward and reverse read files differs by more than one symbol:\n\t{0}\n\t{1}".format(str(forward),
                                                                                                                             str(reverse)))
#---- Initialize tool parameters ----
parameters = config["parameters"][config["parameter_set"]] # short alias for used set of parameters


#---- Configure stages ----
config["stage_list"] = []

# Select configuration and combine stages from all mega_stages in a single list without nesting
if config["mode"] == "preprocessing":
    mega_stage_list = ["preprocessing"]
elif config["mode"] == "qc":
    mega_stage_list = ["preprocessing", "qc"]
elif config["mode"] == "assembly":
    mega_stage_list = ["preprocessing", "qc", "assembly"]
else:
    raise ValueError("ERROR!!! Unknown mode: %s" % config["mode"])

for mega_stage in mega_stage_list:
    custom_megastage_entry = "custom_" + mega_stage + "_stages"
    if (custom_megastage_entry in config) and (config[custom_megastage_entry]):
        config["stage_list"].append(config[custom_megastage_entry])
    else:
        config["stage_list"] += config["allowed_stage_list"][mega_stage][config[mega_stage + "_mode"]][config["starting_point"]]


#---- Save configuration and input files ----
final_config_yaml = output_dict["config"] / "config.final.yaml"
final_input_yaml = output_dict["config"] / "input.final.yaml"

os.makedirs(output_dict["config"], exist_ok=True)

with open(final_config_yaml, 'w') as final_config_fd, open(final_input_yaml, 'w') as final_input_fd:
    yaml.dump(convert_posixpath2str_in_dict(config), final_config_fd, default_flow_style=False, sort_keys=False)
    yaml.dump(convert_posixpath2str_in_dict(input_dict), final_input_fd, default_flow_style=False, sort_keys=False)


#-------------------------------------------
localrules: all

results_list = []

#---- Create output filelist ----


#---- Final rule ----
rule all:
    input:
        results_list
        #results_dict[config["mode"]]
#----

#---- Include section ----
include: "workflow/rules/Preprocessing/Files.smk"



#----
