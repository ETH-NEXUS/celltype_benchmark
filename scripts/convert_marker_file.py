# -*- coding: utf-8 -*-
"""
Making cell typing tool-compatible markers file from Linda's PBMC marker file

Options: garnett, cellassign
##################################
Garnett marker file template:
##################################
>Monocytes
expressed: CD14, FCGR1A, CD68, S100A12
references: https://www.biolegend.com/essential_markers

>B cells
expressed: CD19, MS4A1, CD79A

>T cells
expressed: CD3D, CD3E, CD3G
##################################
##################################
CellAssign marker file template:
##################################
list("cell_type1" = c("gene1","gene2"))
##################################

"""
#%%
marker_file_path = "/cluster/project/nexus/benchmarking/celltyping/marker_files/PBMC/mel_required_celltypes_PBMC.gmx"

#%%
def read_marker_file(marker_file_path):
    marker_file = open(marker_file_path,"r").read()
    names = marker_file.split('\n')[0].split('\t')
    
    """
    Make a dictionary with each population/identity being the key, and a list of marker genes is its value.
    """
    identity_genes = {}
    index = 2
    for population in names:
        identity_genes[population] = list(filter(None,marker_file.split('\n')[index].split('\t'))) # Make sure to filter out blank spaces in gene list
        index += 1
    return identity_genes
#%%
# Simple conversion to Garnett format
def write_garnett(identity_genes):
    out_file = open("/cluster/project/nexus/benchmarking/celltyping/marker_files/PBMC/PBMC_Garnett_markers.txt","w")
    for k,v in identity_genes.items():
        if "(" in k or ")" in k:
            pop_name = str(k).replace("(","_")
            pop_name = pop_name.replace(")","_")
        else:
            pop_name = str(k)
        
        out_file.write(">"+pop_name+"\nexpressed: "+", ".join(v)+"\n\n")
    out_file.close()
#%%
def write_cellassign(identity_genes):
    out_file = open("/cluster/project/nexus/benchmarking/celltyping/marker_files/PBMC/PBMC_CellAssign_markers.txt","w")
    out_file.write(str(list(marker_dict.keys()))+"\n")
    for v in marker_dict.values():
        out_file.write(",".join(v)+"\n")
    out_file.close()
#%%
def write_generic_markerfile(identity_genes):
    out_file = open("/cluster/project/nexus/benchmarking/celltyping/marker_files/PBMC/PBMC_markers.txt","w")
    for k,v in identity_genes.items():
        out_file.write(k+","+",".join(v)+"\n")
    out_file.close()
#%%
marker_dict = read_marker_file(marker_file_path)
write_cellassign(marker_dict)
write_generic_markerfile(marker_dict)