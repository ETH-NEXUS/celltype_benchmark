## Broad PBMC list
import re
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--original_labels', metavar='type=string',
                    help='Input path to file with barcodes and labels.')
parser.add_argument('--label_dictionary', metavar='type=string',
                    help='Input path to file with matched labels dictionary.')
parser.add_argument('--subtype_dictionary', metavar='type=string',
                    help='Input path to file mapping subtypes to major cell types.')
parser.add_argument('--output', help='Output path to final ground truth label file.')
opt = parser.parse_args()

#file_path = "/cluster/work/nexus/benchmarking/celltyping/data/broad/2020-Mar-Census-Adult-Immune-10x/2020-Mar-Census-Adult-Immune-10x_annotated_v1.scp.metadata.txt"
#file_path = "/cluster/work/nexus/benchmarking/celltyping/data/broad/2020-Mar-Census-Newborn-Blood-10x/2020-Mar-Census-Newborn-10x_annotated_v1.scp.metadata.txt"

mode = ""
if("broad" in opt.original_labels.lower()):
    mode = "broad"
elif("zheng" in opt.original_labels.lower()):
    mode = "zheng"
print(mode)

in_list = []
out_list = []

label_dict = {}
with open(opt.label_dictionary) as f:
    for line in f:
       (key, val) = line.strip('\n').split('\t')
       label_dict[key] = val

subtype_dict = {}
with open(opt.subtype_dictionary) as f:
    for line in f:
       (key, val) = line.strip('\n').split('\t')
       subtype_dict[key] = val.split(',')

print(label_dict.items())
print(subtype_dict.items())

'''
label_references = {"CD8":"CD8.T.cells","cytotoxic":"CD8.T.cells","Tem":"CD8.T.cells","effector":"CD8.T.cells",
                    "CD4":"CD4.T.cells","helper":"CD4.T.cells","T.cells.CD4":"CD4.T.cells",
                    "T-helper cell including regulatory T cell":"CD4.T.cells",
                    "regulatory":"T.reg.cells","Treg":"T.reg.cells",
                    "NK":"NK.cells","natural killer":"NK.cells", "CD56":"NK.cells",
                    "CD14":"Monocytes","monocyte":"Monocytes",
                    "macrophage":"Monocytes","myeloid":"Monocytes",
                    "plasmacytoid dendritic cell":"Plasmacytoid.dendritic.cells",
                    "dendritic":"Dendritic.cells",
                    "memory B cell":"B.cells.memory","naive B cell":"B.cells.naive",
                    "precursor B cell":"B.cells.precursor","pro-B cell":"B.cells.precursor",
                    "CD19":"B.cells","B.cell":"B.cells", 
                    "B cell":"B.cells","plasma":"B.cells",
                    "gamma-delta":"Gamma-delta.T.cells",
                    "gamma":"Gamma-delta.T.cells","delta":"Gamma-delta.T.cells",
                    "gd T": "Gamma-delta.T.cells",
                    "CD34": "Stem.cells", "progen":"Stem.cells", "hematopoei":"Stem.cells",
                    "unannotated":"Unknown", "unknown":"Unknown"
                   }

zheng_label_references = {
  "-1":"b_cells_fastqs",
  "-2":"cd14_monocytes_fastqs",
  "-3":"cd34_fastqs",
  "-4":"cd4_t_helper_fastqs",
  "-5":"cd56_nk_fastqs",
  "-6":"cytotoxic_t_fastqs",
  "-7":"memory_t_fastqs",
  "-8":"naive_cytotoxic_fastqs",
  "-9":"naive_t_fastqs",
  "10":"regulatory_t_fastqs",
}
'''

#base_path = "/cluster/project/nexus/benchmarking/celltyping/"
#out_file = open(r"{}ground_truth.tsv".format(base_path),'w')
out_file = open(opt.output,'w')
header = "Barcode"+"\t"+"Original_type"+"\t"+"Coerced_major_type"+"\t"+"Coerced_minor_type"+"\n"
out_file.write(header)
#re.match(pattern="",string="")

# Broad PBMC lists:
if(mode == "broad"):
    with open(opt.original_labels) as f: 
     for line in f:
       barcode = line.split("\t")[0]
       cell_type = line.split("\t")[1]
       in_list.append(cell_type)
      # for k,v in label_references.items():
       for k,v in label_dict.items():
           if re.search(pattern=k, string=cell_type):
               out_list.append(v)
               coerced_type_major = ""
               for major,subtype in subtype_dict.items():
                   print(k,v,major,subtype)
                   if(v in subtype):
                       coerced_type_major = major
                   elif(v in major):
                       coerced_type_major = v
                       v = "none"
                       #print("Major: ",major, " | Subtype: ",v);
                       break
               out_file.write(barcode+"\t"+cell_type+"\t"+coerced_type_major+"\t"+v+"\n")
                  #print(k,"------",cell_type)
               break

    out_file.close()
#print(out_list)
  
elif(mode == "zheng"):
#Zheng:
    with open(barcode_path) as f:
        for barcode in f:
          barcode = line.split("\t")[0]
          cell_type = barcode.strip()[-2:]
          in_list.append(cell_type)
          for k,v in label_references.items():
              if re.search(pattern=k.lower(), string=zheng_label_references[cell_type].lower()):
                  #out_list.append(v)
                  out_file.write(barcode.strip()+"\t"+zheng_label_references[cell_type]+"\t"+v+"\n")
                  print(k,"------",cell_type)
                  break

    out_file.close()



# See how often a cell type occurs in output & see if it matches original list
from collections import Counter

print(Counter(in_list))
#Counter({'naive thymus-derived CD4-positive, alpha-beta T cell': 49169, 'CD14-positive monocyte': 44237, 'cytotoxic T cell': 38076, 'helper T cell': 28814, 'naive B cell': 26547, 'natural killer cell': 20554, 'CD8-positive, alpha-beta T cell': 20464, 'erythrocyte': 14198, 'memory B cell': 9139, 'precursor B cell': 8200, 'unannotated': 6018, 'bone marrow hematopoietic cell': 4603, 'conventional dendritic cell': 4494, 'CD14-positive, CD16-positive monocyte': 4065, 'pro-B cell': 3398, 'plasmacytoid dendritic cell': 2735, 'plasma cell': 2687, 'lymphocyte': 1336, 'megakaryocyte': 901, 'mesenchymal stem cell of the bone marrow': 277, 'annotated_cell_identity.ontology_label': 1, 'group': 1})
# CD8s: 38076 + 20464 = 58540; CD4s: 49169 + 28814 = 77983;
print(Counter(out_list))
#Counter({'CD4_T_cells': 77983, 'CD8_T_cells': 58540, 'B_cells': 52706, 'Monocytes': 48302, 'NK_cells': 20554, 'Unknown': 6018, '': 4494})
