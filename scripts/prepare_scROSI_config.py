import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('--gmx', metavar='type=string',
                    help='Input path to GMX marker file.')
parser.add_argument('--config_groundtruth', help='Input path to ground-truth config file.')
parser.add_argument('--outfile', help='Output path to tsv file to which matching celltype_config will be written.')
opt = parser.parse_args()

# read in header of GMX file
# read in config file
# compare which of the ground truth cell types are present in the header you just extracted

gmx_file = open(opt.gmx).read()
gmx_header = gmx_file.split('\n', 1)[0].split("\t")

config_dictionary = {}
config_subtypes = []
config_major_types = []
with open(opt.config_groundtruth) as config_file:
    next(config_file)
    for line in config_file:
       (major_type, sub_type) = line.split()
       config_dictionary[major_type] = sub_type.split()
       config_major_types.append(major_type)
       config_subtypes.append(sub_type)

config_subtypes = ",".join(set(config_subtypes)).split(",")
config_subtypes.remove("none")
config_major_types = ",".join(set(config_major_types)).split(",")
all_types = config_major_types + config_subtypes
print("GMX cell types: ", gmx_header, "\n")
print("Config file: ", config_dictionary, "\n")

print("Config major cell types: ", config_major_types, "\n")
print("Config minor cell types: ", config_subtypes, "\n")
print("Combined config cell types: ", all_types, "\n")

matched_types = [cell_type for cell_type in gmx_header if cell_type in all_types]
print("Matching cell types: ", matched_types)
print("Cell types in config file but not in GMX: ", [cell_type for cell_type in all_types if not cell_type in gmx_header])
print("Cell types in GMX but not in config file: ", [cell_type for cell_type in gmx_header if not cell_type in all_types])

gmx_table = pd.read_table(opt.gmx, skip_blank_lines = False, usecols = matched_types) # need to make sure cell types in config file are in GMX...
gmx_table.fillna('', inplace=True)
print(gmx_table)
print("Saving modified config file to ", opt.outfile, "...")
gmx_table.to_csv(opt.outfile, sep="\t", index = False)
print("Done.")
