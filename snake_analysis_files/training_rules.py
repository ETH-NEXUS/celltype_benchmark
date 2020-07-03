# Collection of rules necessary to run the celltype benchmark training pipeline

if not 'CV_FOLD_IN' in globals():
    CV_FOLD_IN = INPUTDIR
if not 'CV_FOLD_OUT' in globals():
    CV_FOLD_OUT = OUTDIR + 'cv_k_index/'

# performs the splitting of input files into subsets for the cross validation
rule cross_validation:
        input:
                sce_in = CV_FOLD_IN + "{sample}.genes_cells_filtered.corrected.RDS"
        output:
                barcodes_out = expand(CV_FOLD_OUT + '{{sample}}.indexFold_{k_index}.RDS', k_index = K_FOLD_INDEX),
		barcodes_all = CV_FOLD_OUT + '{sample}.indexAll.RDS'
        params:
                lsfoutfile = CV_FOLD_OUT + '{sample}.cv_fold.lsfout.log',
                lsferrfile = CV_FOLD_OUT + '{sample}.cv_fold.lsferr.log',
                scratch = config['statistics']['cross_validation']['scratch'],
                mem = config['statistics']['cross_validation']['mem'],
                time = config['statistics']['cross_validation']['time'],
                params = config['statistics']['cross_validation']['params'],
                output_dir = CV_FOLD_OUT
        threads:
                config['statistics']['cross_validation']['threads']
        benchmark:
                CV_FOLD_OUT + '{sample}.cv_fold.benchmark'
        shell:
                config['statistics']['cross_validation']['call'] + ' --SCE {input.sce_in} --outputDirec {params.output_dir}'


if not 'MARKER_GENES_IN' in globals():
    MARKER_GENES_IN = INPUTDIR
if not 'MARKER_GENES_OUT' in globals():
    MARKER_GENES_OUT = OUTDIR + 'marker_genes/'

print(MARKER_GENES_IN)
# identify marker genes
rule find_markers:
        input:
                sce_in =  MARKER_GENES_IN + '{sample}.genes_cells_filtered.corrected.RDS'
        output:
                gmx_out = MARKER_GENES_OUT + '{sample}.markers.gmx'
        params:
                lsfoutfile = MARKER_GENES_OUT + '{sample}.find_markers.lsfout.log',
                lsferrfile = MARKER_GENES_OUT + '{sample}.find_markers.lsferr.log',
                scratch = config['markers']['find_markers']['scratch'],
                mem = config['markers']['find_markers']['mem'],
                time = config['markers']['find_markers']['time'],
                params = config['markers']['find_markers']['params'],
                sample_name = '{sample}',
                output_dir = MARKER_GENES_OUT
        threads:
                config['markers']['find_markers']['threads']
        benchmark:
                MARKER_GENES_OUT + '{sample}.find_markers.benchmark'
        shell:
                config['markers']['find_markers']['call'] + ' --SCE {input.sce_in} --outputDirec {params.output_dir}  --sampleName {params.sample_name}'

# convert from gmx to tsv
rule convert_markers:
        input:
                gmx_in = MARKER_GENES_OUT + '{sample}.markers.gmx'
        output:
                tsv_out = MARKER_GENES_OUT + '{sample}.markers.garnett.tsv'
        params:
                lsfoutfile = MARKER_GENES_OUT + '{sample}.convert_markers.lsfout.log',
                lsferrfile = MARKER_GENES_OUT + '{sample}.convert_markers.lsferr.log',
                scratch = config['markers']['convert_markers']['scratch'],
                mem = config['markers']['convert_markers']['mem'],
                time = config['markers']['convert_markers']['time'],
                params = config['markers']['convert_markers']['params']
        threads:
                config['markers']['convert_markers']['threads']
        benchmark:
                MARKER_GENES_OUT + '{sample}.convert_markers.benchmark'
        shell:
                config['markers']['convert_markers']['call'] + ' --gmx_file {input.gmx_in} --tsv_file {output.tsv_out}'

if not 'METHOD_TRAIN_IN' in globals():
    METHOD_TRAIN_IN = CV_FOLD_OUT
if not 'METHOD_TRAIN_OUT' in globals():
    METHOD_TRAIN_OUT = OUTDIR + 'crossvalidation/'
# method rf; training crossvalidaton
rule rf_train:
        input:
                sce_in =  INPUTDIR + '{sample}.genes_cells_filtered.corrected.RDS',
                barcodes = METHOD_TRAIN_IN + '{sample}.indexFold_{k_index}.RDS'
        output:
                pred_labels = METHOD_TRAIN_OUT + '{sample}.indexFold_{k_index}.rf_predicted_labels.tsv'
        params:
                lsfoutfile = METHOD_TRAIN_OUT + '{sample}.indexFold_{k_index}.rf_train.lsfout.log',
                lsferrfile = METHOD_TRAIN_OUT + '{sample}.indexFold_{k_index}.rf_train.lsferr.log',
                scratch = config['methods']['rf']['train']['scratch'],
                mem = config['methods']['rf']['train']['mem'],
                time = config['methods']['rf']['train']['time'],
                params = config['methods']['rf']['train']['params'],
                output_dir = METHOD_TRAIN_OUT,
                sample_name = "{sample}"
        threads:
                config['methods']['rf']['train']['threads']
        benchmark:
                METHOD_TRAIN_OUT + '{sample}.indexFold_{k_index}.rf_train.benchmark'
        shell:
                config['methods']['rf']['train']['call'] + ' --SCE {input.sce_in} --CVindex {input.barcodes} --outputDirec {params.output_dir} --sampleName {params.sample_name}'

if not 'METHOD_ALL_TRAIN_OUT' in globals():
    METHOD_ALL_TRAIN_OUT = OUTDIR + "trained_model/"
# method rf, training complete model
rule rf_train_all:
        input:
                sce_in = INPUTDIR + '{sample}.genes_cells_filtered.corrected.RDS',
                barcodes = METHOD_TRAIN_IN + '{sample}.indexAll.RDS'
        output:
                pred_labels = METHOD_ALL_TRAIN_OUT + '{sample}.rf_predicted_labels_all.tsv'
        params:
                lsfoutfile = METHOD_ALL_TRAIN_OUT + '{sample}.rf_train_all.lsfout.log',
                lsferrfile = METHOD_ALL_TRAIN_OUT + '{sample}.rf_train_all.lsferr.log',
                scratch = config['methods']['rf']['train']['scratch'],
                mem = config['methods']['rf']['train']['mem'],
                time = config['methods']['rf']['train']['time'],
                params = config['methods']['rf']['train']['params'],
                output_dir = METHOD_ALL_TRAIN_OUT,
                sample_name = "{sample}"
        threads:
                config['methods']['rf']['train']['threads']
        benchmark:
                METHOD_ALL_TRAIN_OUT + '{sample}.rf_train_all.benchmark'
        shell:
                config['methods']['rf']['train']['call'] + ' --SCE {input.sce_in} --CVindex {input.barcodes} --outputDirec {params.output_dir} --sampleName {params.sample_name}'
# method svm, training crossvalidation
rule svm_train:
        input:
                sce_in =  INPUTDIR + '{sample}.genes_cells_filtered.corrected.RDS',
                barcodes = METHOD_TRAIN_IN + '{sample}.indexFold_{k_index}.RDS'
        output:
                pred_labels = METHOD_TRAIN_OUT + '{sample}.indexFold_{k_index}.svm_predicted_labels.tsv'
        params:
                lsfoutfile = METHOD_TRAIN_OUT + '{sample}.indexFold_{k_index}.svm_train.lsfout.log',
                lsferrfile = METHOD_TRAIN_OUT + '{sample}.indexFold_{k_index}.svm_train.lsferr.log',
                scratch = config['methods']['svm']['train']['scratch'],
                mem = config['methods']['svm']['train']['mem'],
                time = config['methods']['svm']['train']['time'],
                params = config['methods']['svm']['train']['params'],
                output_dir = METHOD_TRAIN_OUT,
                sample_name = "{sample}"
        threads:
                config['methods']['svm']['train']['threads']
        benchmark:
                METHOD_TRAIN_OUT + '{sample}.indexFold_{k_index}.svm_train.benchmark'
        shell:
                config['methods']['svm']['train']['call'] + ' --SCE {input.sce_in} --CVindex {input.barcodes} --outputDirec {params.output_dir} --sampleName {params.sample_name}'

# method svm, training, complete model
rule svm_train_all:
        input:
                sce_in = INPUTDIR + '{sample}.genes_cells_filtered.corrected.RDS',
                barcodes = METHOD_TRAIN_IN + '{sample}.indexAll.RDS'
        output:
                pred_labels = METHOD_ALL_TRAIN_OUT + '{sample}.svm_predicted_labels_all.tsv'
        params:
                lsfoutfile = METHOD_ALL_TRAIN_OUT + '{sample}.svm_train_all.lsfout.log',
                lsferrfile = METHOD_ALL_TRAIN_OUT + '{sample}.svm_train_all.lsferr.log',
                scratch = config['methods']['svm']['train']['scratch'],
                mem = config['methods']['svm']['train']['mem'],
                time = config['methods']['svm']['train']['time'],
                params = config['methods']['svm']['train']['params'],
                output_dir = METHOD_ALL_TRAIN_OUT,
                sample_name = "{sample}"
        threads:
                config['methods']['svm']['train']['threads']
        benchmark:
                METHOD_ALL_TRAIN_OUT + '{sample}.svm_train_all.benchmark'
        shell:
                config['methods']['svm']['train']['call'] + ' --SCE {input.sce_in} --CVindex {input.barcodes} --outputDirec {params.output_dir} --sampleName {params.sample_name}'

# helper script - converts SCE to CDS (input for Garnett)
rule sce_to_cds:
        input:
                sce_in = INPUTDIR + '{sample}.genes_cells_filtered.corrected.RDS'
        output:
                cds_out = CDS_OUT + '{sample}.monocle3_cds.RDS'
        params:
                lsfoutfile = METHOD_TRAIN_OUT + '{sample}.sce_to_cds.lsfout.log',
                lsferrfile = METHOD_TRAIN_OUT + '{sample}.sce_to_cds.lsferr.log',
                scratch = config['preprocessing']['sce_to_cds']['scratch'],
                mem = config['preprocessing']['sce_to_cds']['mem'],
                time = config['preprocessing']['sce_to_cds']['time'],
                params = config['preprocessing']['sce_to_cds']['params'],
                output_dir = CDS_OUT
                sample_name = "{sample}"
        threads:
                config['preprocessing']['sce_to_cds']['threads']
        conda:
                config['methods']['garnett']['conda']
        benchmark:
                METHOD_TRAIN_OUT + '{sample}.indexFold_{k_index}.garnett_train.benchmark'
        shell:
                config['preprocessing']['sce_to_cds']['call'] + ' --sce_in {input.sce_in} --sample_name {params.sample_name} --cds_out {output.cds_out}'

# method garnett, training, crossvalidation
rule garnett_train:
        input:
                cds_in =  CDS_OUT + '{sample}.monocle3_cds.RDS',
                barcodes = METHOD_TRAIN_IN + '{sample}.indexFold_{k_index}.RDS'
        output:
                pred_labels = METHOD_TRAIN_OUT + '{sample}.indexFold_{k_index}.garnett_predicted_labels.tsv'
        params:
                lsfoutfile = METHOD_TRAIN_OUT + '{sample}.indexFold_{k_index}.garnett_train.lsfout.log',
                lsferrfile = METHOD_TRAIN_OUT + '{sample}.indexFold_{k_index}.garnett_train.lsferr.log',
                scratch = config['methods']['garnett']['train']['scratch'],
                mem = config['methods']['garnett']['train']['mem'],
                time = config['methods']['garnett']['train']['time'],
                params = config['methods']['garnett']['train']['params'],
                output_dir = METHOD_TRAIN_OUT,
                sample_name = "{sample}"
        threads:
                config['methods']['garnett']['train']['threads']
        conda:
                config['methods']['garnett']['conda']
        benchmark:
                METHOD_TRAIN_OUT + '{sample}.indexFold_{k_index}.garnett_train.benchmark'
        shell:
                config['methods']['garnett']['train']['call'] + ' --cds {input.cds_in} --barcodes_index {input.barcodes} --output_dir {params.output_dir} --pred_labels_out {output.pred_labels} --sample_name {params.sample_name}'

# method garnett, training, complete model
rule garnett_train_all:
        input:
                cds_in =  CDS_OUT + '{sample}.monocle3_cds.RDS'
        output:
                pred_labels = METHOD_ALL_TRAIN_OUT + '{sample}.garnett_predicted_labels_all.tsv'
        params:
                lsfoutfile = METHOD_ALL_TRAIN_OUT + '{sample}.garnett_train_all.lsfout.log',
                lsferrfile = METHOD_ALL_TRAIN_OUT + '{sample}.garnett_train_all.lsferr.log',
                scratch = config['methods']['garnett']['train']['scratch'],
                mem = config['methods']['garnett']['train']['mem'],
                time = config['methods']['garnett']['train']['time'],
                params = config['methods']['garnett']['train']['params'],
                output_dir = METHOD_ALL_TRAIN_OUT,
                sample_name = "{sample}"
        threads:
                config['methods']['garnett']['train']['threads']
        conda:
                config['methods']['garnett']['conda']
        benchmark:
                METHOD_ALL_TRAIN_OUT + '{sample}.garnett_train_all.benchmark'
        shell:
                config['methods']['garnett']['train']['call'] + ' --cds {input.cds_in} --output_dir {params.output_dir} --pred_labels_out {output.pred_labels} --sample_name {params.sample_name}'

# helper function to retrieve files necessary for the crossvalidation evaluation
def getTrainedSubsets(wildcards):
        return expand(METHOD_TRAIN_OUT + wildcards.sample + '.indexFold_{k_index}.{method}_predicted_labels.tsv', k_index = K_FOLD_INDEX, method = SELECTED_METHODS)

# get the crossvalidation summary
rule cv_summary_train:
        input:
                cv_files = getTrainedSubsets,
                sce_in = METHOD_TRAIN_IN + "{sample}.indexAll.RDS"
        output:
                METHOD_TRAIN_OUT + '{sample}.cv_summary.tsv'
        params:
                lsfoutfile = METHOD_TRAIN_OUT + '{sample}.cv_summary.lsfout.log',
                lsferrfile = METHOD_TRAIN_OUT + '{sample}.cv_summary.lsferr.log',
                scratch = config['statistics']['cv_summary']['scratch'],
                mem = config['statistics']['cv_summary']['mem'],
                time = config['statistics']['cv_summary']['time'],
                params = config['statistics']['cv_summary']['params'],
                output_dir = METHOD_TRAIN_OUT,
                sample_name = "{sample}"
        threads:
                config['statistics']['cv_summary']['threads']
        benchmark:
                METHOD_TRAIN_OUT + '{sample}.cv_summary.benchmark'
        shell:
                config['statistics']['cv_summary']['call'] + ' --CVfiles {input.cv_files} --SCE {input.sce_in} --outputDirec {params.output_dir} --sampleName {params.sample_name}'
