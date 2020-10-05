import argparse
import pandas as pd
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('--gmx', metavar='type=string',
                    help='Input path to GMX marker file.')
parser.add_argument('--config', help='Input path to ground-truth config file.')
parser.add_argument('--output_dir', help='Output directory in which to write matching celltype_config.')
parser.add_argument('--output_name', help='Output filename to prepend to GMX.')
opt = parser.parse_args()

# read in header of GMX file
# read in config file
# compare which of the ground truth cell types are present in the header you just extracted

gmx_file = open(opt.gmx).read()
gmx_header = gmx_file.split('\n', 1)[0].split("\t")

gmx_out_path = opt.output_dir + opt.output_name + "trimmed.gmx"

config_dictionary = {}
config_subtypes = []
config_major_types = []
with open(opt.config) as config_file:
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

index = 0
cell_types_kept_indices = []
for each in gmx_header:
    if each in all_types:
        print("Matched cell types: GMX - ", each)
        cell_types_kept_indices.append(index)
    else:
        print("Not found: ",each)
    index += 1

print(cell_types_kept_indices)
cell_types_kept_indices_string = [i+1 for i in cell_types_kept_indices]
cell_types_kept_indices_string = [str(index) for index in cell_types_kept_indices_string]
cell_types_kept_indices_string = ",".join(cell_types_kept_indices_string)
print("\n","Cell types kept: ", [gmx_header[cell_type] for cell_type in cell_types_kept_indices], "\n")
print("GMX columns kept: ", cell_types_kept_indices_string, "\n")

#import pandas
gmx_out_path = opt.output_dir + opt.output_name + "_trimmed.gmx"

gmx_table = pd.read_table(opt.gmx, skip_blank_lines = False, usecols = cell_types_kept_indices)
gmx_table.fillna('', inplace=True)
print(gmx_table)
print("Saving Pandas DF as tsv...")
gmx_table.to_csv((gmx_out_path+".pd"), sep="\t", index = False)
print("Done.")
#simplified#
all_types = [cell_type for cell_type in gmx_header if cell_type in all_types]
gmx_table = pd.read_table(opt.gmx, skip_blank_lines = False, usecols = all_types) # need to make sure cell types in config file are in GMX...
gmx_table.fillna('', inplace=True)
print(gmx_table)
#import subprocess
def cut_gmx():
    cut_command = ["cut","-d",r"	","-f",f"{cell_types_kept_indices_string}",f"{opt.gmx}"]
    gmx_out = open(gmx_out_path, "w")
    print("Saving modified GMX file to " + gmx_out_path  + "...")
    subprocess.call(cut_command, stdout = gmx_out)
    gmx_out.close()
