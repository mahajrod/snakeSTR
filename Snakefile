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

for subdir in config["out_subdir_list"]:
    output_dict[subdir] = out_dir_path / subdir

#--------

#----

#---- Checking input files ----
#logging.info("Checking input files...")
sample_id_list = []
for path in (input_dir_path / "alignment" ).glob("*"):
    if path.is_dir():
        sample_id_list.append(path.name)

if not sample_id_list:
    raise ValueError("ERROR!!! No samples were detected...")

absent_bam_sample_list = []

for sample in sample_id_list:
    if not (input_dir_path / "alignment" / sample / (sample + config["bam_suffix"] + ".bam")).exists():
        absent_bam_sample_list.append(sample)

if absent_bam_sample_list:
    raise ValueError("ERROR! Bam files we not found for following samples: {0}".format(",".join(absent_bam_sample_list)))

# check presence of bed files
str_loci_bed_path = Path(config["str_loci_bed"])
if not str_loci_bed_path.exists(): # check if file is located in the input_folder
    str_loci_bed_path = input_dir_path / "bed"/ config["str_loci_bed"]

str_repeat_bed_path = Path(config["str_repeat_bed"])
if not str_repeat_bed_path.exists():  # check if file is located in the input_folder
    str_repeat_bed_path = input_dir_path / "bed" / config["str_repeat_bed"]

reference_path = Path(config["reference"])
if not reference_path.exists():  # check if file is located in the input_folder
    reference_path = input_dir_path / "reference" / config["reference"]


#---- Initialize tool parameters ----
parameters = config["parameters"][config["parameter_set"]] # short alias for used set of parameters



#---- Save configuration and input files ----
final_config_yaml = out_dir_path / "config/config.final.yaml"

os.makedirs(out_dir_path / "config/", exist_ok=True)

with open(final_config_yaml, 'w') as final_config_fd:
    yaml.dump(convert_posixpath2str_in_dict(config), final_config_fd, default_flow_style=False, sort_keys=False)

#-------------------------------------------
localrules: all

results_list = [out_dir_path / "str/hipSTR.filtered.vcf",
                expand(out_dir_path / "admixture/structure/{stage}/structure.K{K}.R{run}_f",
                      stage=["raw", "filtered"],
                      K=parameters["tool_options"]["structure"]["K_list"],
                      run=[0, 1, 2])]

#---- Create output filelist ----


#---- Final rule ----
rule all:
    input:
        results_list
#----

#---- Include section ----
include: "workflow/rules/Preprocessing/Files.smk"
include: "workflow/rules/Alignment/Samtools.smk"
include: "workflow/rules/Alignment/Bazam.smk"
include: "workflow/rules/Alignment/Alignment.smk"
include: "workflow/rules/STR/HipSTR.smk"
include: "workflow/rules/Admixture/STRUCTURE.smk"

#----
