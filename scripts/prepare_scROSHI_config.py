############### prepare_scROSHI_config.py ###############
# Author: Marcus Lindberg
# Date: 05/10/2020
#
# Requirements: Python 3.7.4
#
# A simple script for filtering celltype_config files based on cell types present in a GMX marker file.
# As input for scROSHI, it is critical that the celltype_config file does not contain cell types not present in the GMX file.
# The final config file is saved as a tsv file.
#
# Example usage:
# python prepare_scROSHI_config.py --gmx {marker_file.gmx} --config_groundtruth {config_file.tsv} --outfile {modified_config_file.tsv}

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--gmx', metavar='type=string',
                    help='Input path to GMX marker file.')
parser.add_argument('--config_groundtruth', help='Input path to ground-truth config file.')
parser.add_argument('--outfile', help='Output path to tsv file to which matching celltype_config will be written.')
opt = parser.parse_args()

# Read in header of GMX marker file.
gmx_file = open(opt.gmx).read()
gmx_header = gmx_file.split('\n', 1)[0].split("\t")

# Read config file and create a dictionary mapping major cell types to minor subtypes.
config_dictionary = {}
config_subtypes = []
config_major_types = []
with open(opt.config_groundtruth) as config_file:
    next(config_file)
    for line in config_file:
       (major_type, sub_type) = line.split()
       config_dictionary[major_type] = sub_type.split(",")
       config_major_types.append(major_type)
       config_subtypes.append(sub_type)

config_subtypes = ",".join(set(config_subtypes)).split(",")
if "none" in config_subtypes:
    config_subtypes.remove("none")
config_major_types = ",".join(set(config_major_types)).split(",")
all_types = config_major_types + config_subtypes

# Summarize cell types present within GMX and config file.
print("\n", "#"*20, "\n", "Summary of GMX and config file:", "\n", "#"*20, "\n")
print("GMX cell types: ", gmx_header, "\n")
print("Config file: ", config_dictionary, "\n")

print("Config major cell types: ", config_major_types, "\n")
print("Config minor cell types: ", config_subtypes, "\n")
print("Combined config cell types: ", all_types, "\n")

# Match cell types found in both GMX and config file.
print("\n", "#"*20, "\n", "Matches between GMX and config file:", "\n", "#"*20, "\n")
matched_types = [cell_type for cell_type in gmx_header if cell_type in all_types]
print("Matching cell types: ", matched_types)

types_in_config_only = [cell_type for cell_type in all_types if not cell_type in gmx_header]
types_in_gmx_only = [cell_type for cell_type in gmx_header if not cell_type in all_types]

# Validate config file integrity.
if types_in_config_only:
    print("Cell types in config file but not in GMX: ", types_in_config_only, "\n")
if types_in_gmx_only:
    print("\033[1m" + "WARNING: Config file is missing cell types found in GMX. Make sure config file is correct." + "\033[0m")
    print("Cell types in GMX but not in config file: ", [cell_type for cell_type in gmx_header if not cell_type in all_types])

# Save final config file as tsv file.
print("\n", "#"*20, "\n", "Final config file:", "\n", "#"*20, "\n")
config_outfile = open(opt.outfile, "w")
outfile_header = "Major" + "\t" + "Subtype"
config_outfile.write(outfile_header + "\n")
print(outfile_header)
for cell_type in matched_types:
    if cell_type in config_major_types:
        subtypes = [subtype for subtype in config_dictionary[cell_type] if subtype in matched_types]
        if not len(subtypes):
            subtypes = ["none"]
        outfile_line = cell_type + "\t" + ",".join(subtypes)
        config_outfile.write(outfile_line + "\n")
        print(outfile_line)

config_outfile.close()
