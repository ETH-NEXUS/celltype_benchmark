# Collection of rules necessary to run the celltype benchmark training pipeline

if not 'CV_FOLD_IN' in globals():
    CV_FOLD_IN = INPUTDIR
if not 'CV_FOLD_OUT' in globals():
    CV_FOLD_OUT = OUTDIR + 'cv_k_index/'

# performs the splitting of input files into subsets for the cross validation
rule cross_validation:
        input:
                sce_in = CV_FOLD_IN + "{sample}.RDS"
        output:
                barcodes_out = expand(CV_FOLD_OUT + '{{sample}}.indexFold_{k_index}.RDS', k_index = K_FOLD_INDEX),
		barcodes_all = CV_FOLD_OUT + '{sample}.index_All.RDS'
        params:
                lsfoutfile = CV_FOLD_OUT + '{sample}.cv_fold.lsfout.log',
                lsferrfile = CV_FOLD_OUT + '{sample}.cv_fold.lsferr.log',
                scratch = config['statistics']['cross_validation']['scratch'],
                mem = config['statistics']['cross_validation']['mem'],
                time = config['statistics']['cross_validation']['time'],
                params = config['statistics']['cross_validation']['params'],
                k_folds = config['statistics']['cross_validation']['k_folds'],
                sample_name = "{sample}",
                output_dir = CV_FOLD_OUT
        threads:
                config['statistics']['cross_validation']['threads']
        benchmark:
                CV_FOLD_OUT + '{sample}.cv_fold.benchmark'
        shell:
                config['statistics']['cross_validation']['call'] + ' --SCE {input.sce_in} --sample_name {params.sample_name} --outputDirec {params.output_dir} --kfold {params.k_folds}'


if not 'MARKER_GENES_IN' in globals():
    MARKER_GENES_IN = INPUTDIR
if not 'MARKER_GENES_OUT' in globals():
    MARKER_GENES_OUT = OUTDIR + 'marker_genes/'

#print(MARKER_GENES_IN)
# identify marker genes
rule find_markers:
        input:
                sce_in =  MARKER_GENES_IN + '{sample}.RDS'
        output:
                gmx_out = MARKER_GENES_OUT + '{sample}.gmx',
		tsv_garnett_out = MARKER_GENES_OUT + '{sample}.garnett.tsv'
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

if not 'METHOD_TRAIN_IN' in globals():
    METHOD_TRAIN_IN = CV_FOLD_OUT
if not 'METHOD_TRAIN_OUT' in globals():
    METHOD_TRAIN_OUT = OUTDIR + 'crossvalidation/'
# method rf; training crossvalidaton
rule rf_train:
        input:
                sce_in =  INPUTDIR + '{sample}.RDS',
                barcodes = METHOD_TRAIN_IN + '{sample}.indexFold_{k_index}.RDS'
        output:
                pred_labels = METHOD_TRAIN_OUT + '{sample}.indexFold_{k_index}.rf_predicted_labels.csv'
        params:
                lsfoutfile = METHOD_TRAIN_OUT + '{sample}.indexFold_{k_index}.rf_train.lsfout.log',
                lsferrfile = METHOD_TRAIN_OUT + '{sample}.indexFold_{k_index}.rf_train.lsferr.log',
                scratch = config['methods']['rf']['train']['scratch'],
                mem = config['methods']['rf']['train']['mem'],
                time = config['methods']['rf']['train']['time'],
                params = config['methods']['rf']['train']['params'],
                output_dir = METHOD_TRAIN_OUT,
                sample_name = "{sample}.indexFold_{k_index}",
                method_name = "rf"
        threads:
                config['methods']['rf']['train']['threads']
        benchmark:
                METHOD_TRAIN_OUT + '{sample}.indexFold_{k_index}.rf_train.benchmark'
        shell:
                config['methods']['rf']['train']['call'] + ' --SCE {input.sce_in} --CVindex {input.barcodes} --outputDirec {params.output_dir} --sampleName {params.sample_name} --method {params.method_name}'

if not 'METHOD_ALL_TRAIN_OUT' in globals():
    METHOD_ALL_TRAIN_OUT = OUTDIR + "trained_model/"
# method rf, training complete model
rule rf_train_all:
        input:
                sce_in = INPUTDIR + '{sample}.RDS',
                barcodes = METHOD_TRAIN_IN + '{sample}.index_All.RDS'
        output:
                model_out = METHOD_ALL_TRAIN_OUT + '{sample}.rf_model.RDS'
        params:
                lsfoutfile = METHOD_ALL_TRAIN_OUT + '{sample}.rf_train_all.lsfout.log',
                lsferrfile = METHOD_ALL_TRAIN_OUT + '{sample}.rf_train_all.lsferr.log',
                scratch = config['methods']['rf']['train']['scratch'],
                mem = config['methods']['rf']['train']['mem'],
                time = config['methods']['rf']['train']['time'],
                params = config['methods']['rf']['train']['params'],
                output_dir = METHOD_ALL_TRAIN_OUT,
                sample_name = "{sample}",
                method_name = "rf"
        threads:
                config['methods']['rf']['train']['threads']
        benchmark:
                METHOD_ALL_TRAIN_OUT + '{sample}.rf_train_all.benchmark'
        shell:
                config['methods']['rf']['train']['call'] + ' --SCE {input.sce_in} --CVindex {input.barcodes} --outputDirec {params.output_dir} --sampleName {params.sample_name} --method {params.method_name}'
# method svm, training crossvalidation
rule svm_train:
        input:
                sce_in =  INPUTDIR + '{sample}.RDS',
                barcodes = METHOD_TRAIN_IN + '{sample}.indexFold_{k_index}.RDS'
        output:
                pred_labels = METHOD_TRAIN_OUT + '{sample}.indexFold_{k_index}.svm_predicted_labels.csv'
        params:
                lsfoutfile = METHOD_TRAIN_OUT + '{sample}.indexFold_{k_index}.svm_train.lsfout.log',
                lsferrfile = METHOD_TRAIN_OUT + '{sample}.indexFold_{k_index}.svm_train.lsferr.log',
                scratch = config['methods']['svm']['train']['scratch'],
                mem = config['methods']['svm']['train']['mem'],
                time = config['methods']['svm']['train']['time'],
                params = config['methods']['svm']['train']['params'],
                output_dir = METHOD_TRAIN_OUT,
                sample_name = "{sample}.indexFold_{k_index}",
                method_name = "svm"
        threads:
                config['methods']['svm']['train']['threads']
        benchmark:
                METHOD_TRAIN_OUT + '{sample}.indexFold_{k_index}.svm_train.benchmark'
        shell:
                config['methods']['svm']['train']['call'] + ' --SCE {input.sce_in} --CVindex {input.barcodes} --outputDirec {params.output_dir} --sampleName {params.sample_name} --method {params.method_name}'

# method svm, training, complete model
rule svm_train_all:
        input:
                sce_in = INPUTDIR + '{sample}.RDS',
                barcodes = METHOD_TRAIN_IN + '{sample}.index_All.RDS'
        output:
                model_out = METHOD_ALL_TRAIN_OUT + '{sample}.svm_model.RDS'
        params:
                lsfoutfile = METHOD_ALL_TRAIN_OUT + '{sample}.svm_train_all.lsfout.log',
                lsferrfile = METHOD_ALL_TRAIN_OUT + '{sample}.svm_train_all.lsferr.log',
                scratch = config['methods']['svm']['train']['scratch'],
                mem = config['methods']['svm']['train']['mem'],
                time = config['methods']['svm']['train']['time'],
                params = config['methods']['svm']['train']['params'],
                output_dir = METHOD_ALL_TRAIN_OUT,
                sample_name = "{sample}",
                method_name = "svm"
        threads:
                config['methods']['svm']['train']['threads']
        benchmark:
                METHOD_ALL_TRAIN_OUT + '{sample}.svm_train_all.benchmark'
        shell:
                config['methods']['svm']['train']['call'] + ' --SCE {input.sce_in} --CVindex {input.barcodes} --outputDirec {params.output_dir} --sampleName {params.sample_name} --method {params.method_name}'

# method fcNN, training crossvalidation
rule fcnn_train:
        input:
                sce_in =  INPUTDIR + '{sample}.RDS',
                barcodes = METHOD_TRAIN_IN + '{sample}.indexFold_{k_index}.RDS'
        output:
                pred_labels = METHOD_TRAIN_OUT + '{sample}.indexFold_{k_index}.fcnn_predicted_labels.csv'
        params:
                lsfoutfile = METHOD_TRAIN_OUT + '{sample}.indexFold_{k_index}.fcnn_train.lsfout.log',
                lsferrfile = METHOD_TRAIN_OUT + '{sample}.indexFold_{k_index}.fcnn_train.lsferr.log',
                scratch = config['methods']['fcnn']['train']['scratch'],
                mem = config['methods']['fcnn']['train']['mem'],
                time = config['methods']['fcnn']['train']['time'],
                params = config['methods']['fcnn']['train']['params'],
                output_dir = METHOD_TRAIN_OUT,
                sample_name = "{sample}.indexFold_{k_index}",
                method_name = "fcnn"
        threads:
                config['methods']['fcnn']['train']['threads']
        benchmark:
                METHOD_TRAIN_OUT + '{sample}.indexFold_{k_index}.fcnn_train.benchmark'
        shell:
                config['methods']['fcnn']['train']['call'] + ' --SCE {input.sce_in} --CVindex {input.barcodes} --outputDirec {params.output_dir} --sampleName {params.sample_name} --method {params.method_name}'

# method fcNN, training, complete model
rule fcnn_train_all:
        input:
                sce_in = INPUTDIR + '{sample}.RDS',
                barcodes = METHOD_TRAIN_IN + '{sample}.index_All.RDS'
        output:
                model_out = METHOD_ALL_TRAIN_OUT + '{sample}.fcnn_model.h5'
        params:
                lsfoutfile = METHOD_ALL_TRAIN_OUT + '{sample}.fcnn_train_all.lsfout.log',
                lsferrfile = METHOD_ALL_TRAIN_OUT + '{sample}.fcnn_train_all.lsferr.log',
                scratch = config['methods']['fcnn']['train']['scratch'],
                mem = config['methods']['fcnn']['train']['mem'],
                time = config['methods']['fcnn']['train']['time'],
                params = config['methods']['fcnn']['train']['params'],
                output_dir = METHOD_ALL_TRAIN_OUT,
                sample_name = "{sample}",
                method_name = "fcnn"
        threads:
                config['methods']['fcnn']['train']['threads']
        benchmark:
                METHOD_ALL_TRAIN_OUT + '{sample}.fcnn_train_all.benchmark'
        shell:
                config['methods']['fcnn']['train']['call'] + ' --SCE {input.sce_in} --CVindex {input.barcodes} --outputDirec {params.output_dir} --sampleName {params.sample_name} --method {params.method_name}'

# helper script - converts SCE to CDS (input for Garnett)
rule sce_to_cds:
        input:
                sce_in = INPUTDIR + '{sample}.RDS'
        output:
                cds_out = CDS_OUT + '{sample}.monocle3_cds.RDS'
        params:
                lsfoutfile = CDS_OUT + '{sample}.sce_to_cds.lsfout.log',
                lsferrfile = CDS_OUT + '{sample}.sce_to_cds.lsferr.log',
                scratch = config['preprocessing']['sce_to_cds']['scratch'],
                mem = config['preprocessing']['sce_to_cds']['mem'],
                time = config['preprocessing']['sce_to_cds']['time'],
                params = config['preprocessing']['sce_to_cds']['params'],
                output_dir = CDS_OUT,
                sample_name = "{sample}",
                features_filepath = config['preprocessing']['sce_to_cds']['params']['features_filepath']
        threads:
                config['preprocessing']['sce_to_cds']['threads']
        conda:
                config['methods']['garnett']['conda']
        benchmark:
                CDS_OUT + '{sample}.sce_to_cds.benchmark'
        shell:
                config['preprocessing']['sce_to_cds']['call'] + ' --sce_in {input.sce_in} --sample_name {params.sample_name} --output_dir {params.output_dir} --features_file {params.features_filepath}'

# method garnett, training, crossvalidation
rule garnett_train:
        input:
                cds_in =  CDS_OUT + '{sample}.monocle3_cds.RDS',
                barcodes = METHOD_TRAIN_IN + '{sample}.indexFold_{k_index}.RDS',
                marker_file = MARKER_GENES_OUT + '{sample}.garnett.tsv'
        output:
                pred_labels = METHOD_TRAIN_OUT + '{sample}.indexFold_{k_index}.garnett_predicted_labels.csv'
        params:
                lsfoutfile = METHOD_TRAIN_OUT + '{sample}.indexFold_{k_index}.garnett_train.lsfout.log',
                lsferrfile = METHOD_TRAIN_OUT + '{sample}.indexFold_{k_index}.garnett_train.lsferr.log',
                scratch = config['methods']['garnett']['train']['scratch'],
                mem = config['methods']['garnett']['train']['mem'],
                time = config['methods']['garnett']['train']['time'],
                params = config['methods']['garnett']['train']['params'],
                output_dir = METHOD_TRAIN_OUT,
                sample_name = "{sample}.indexFold_{k_index}"
        threads:
                config['methods']['garnett']['train']['threads']
        conda:
                config['methods']['garnett']['conda']
        benchmark:
                METHOD_TRAIN_OUT + '{sample}.indexFold_{k_index}.garnett_train.benchmark'
        shell:
                config['methods']['garnett']['train']['call'] + ' --cds {input.cds_in} --barcodes_index {input.barcodes} --output_dir {params.output_dir} --marker_file {input.marker_file} --sample_name {params.sample_name}'

# method garnett, training, complete model
rule garnett_train_all:
        input:
                cds_in =  CDS_OUT + '{sample}.monocle3_cds.RDS',
                barcodes = METHOD_TRAIN_IN + '{sample}.index_All.RDS',
                marker_file = MARKER_GENES_OUT + '{sample}.garnett.tsv'
        output:
                model_out = METHOD_ALL_TRAIN_OUT + '{sample}.garnett_model.RDS'
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
                config['methods']['garnett']['train']['call'] + ' --cds {input.cds_in} --barcodes_index {input.barcodes} --output_dir {params.output_dir} --marker_file {input.marker_file} --sample_name {params.sample_name}'

# helper function to retrieve files necessary for the crossvalidation evaluation
def getTrainedSubsets(wildcards):
        return expand(METHOD_TRAIN_OUT + wildcards.sample + '.indexFold_{k_index}.{method}_predicted_labels.csv', k_index = K_FOLD_INDEX, method = SELECTED_METHODS)

if not 'CV_SUMMARY_IN' in globals():
    CV_SUMMARY_IN = INPUTDIR
if not 'CV_SUMMARY_OUT' in globals():
    CV_SUMMARY_OUT = OUTDIR + "cv_summary/

# get the crossvalidation summary
rule cv_summary_train:
        input:
                cv_files = METHOD_TRAIN_OUT,
                sce_in = CV_SUMMARY_IN + "{sample}.RDS",
        output:
                accuracy = CV_SUMMARY_OUT + '{sample}_byClass.csv',
                f1_score = CV_SUMMARY_OUT + '{sample}_F1.csv',
                overall = CV_SUMMARY_OUT + '{sample}_overall.csv'
        params:
                lsfoutfile = CV_SUMMARY_OUT + '{sample}.cv_summary.lsfout.log',
                lsferrfile = CV_SUMMARY_OUT + '{sample}.cv_summary.lsferr.log',
                scratch = config['statistics']['cv_summary']['scratch'],
                mem = config['statistics']['cv_summary']['mem'],
                time = config['statistics']['cv_summary']['time'],
                params = config['statistics']['cv_summary']['params'],
                output_dir = CV_SUMMARY_OUT,
                barcodes_dir = CV_FOLD_OUT,
                sample_name = "{sample}",
                tools = ",".join(config['resources']['celltype_methods'])
        threads:
                config['statistics']['cv_summary']['threads']
        benchmark:
                CV_SUMMARY_OUT + '{sample}.cv_summary.benchmark'
        shell:
                config['statistics']['cv_summary']['call'] + ' --CVfiles {input.cv_files} --CVIndex {params.barcodes_dir} --SCE {input.sce_in} --outputDirec {params.output_dir} --sampleName {params.sample_name} --method {params.tools}'
