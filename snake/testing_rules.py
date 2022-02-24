# Collection of rules necessary to run the celltype benchmark testing pipeline

if not 'CELLTYPING_IN' in globals():
    CELLTYPING_IN = OUTDIR + "trained_model/"
if not 'CELLTYPING_OUT' in globals():
    CELLTYPING_OUT = OUTDIR + 'celltyping/'


# method rf, testing
rule rf_test:
        input:
                sce_in =  INPUTDIR + '{sample}.RDS',
                trained_model = CELLTYPING_IN + '{model}.rf_model.RDS'
        output:
                pred_labels = CELLTYPING_OUT + '{sample}/rf/predicted_labels.{model}.csv'
        params:
                lsfoutfile = CELLTYPING_OUT + '{sample}/rf/predicted_labels.{model}.lsfout.log',
                lsferrfile = CELLTYPING_OUT + '{sample}/rf/predicted_labels.{model}.lsferr.log',
                scratch = config['methods']['rf']['test']['scratch'],
                mem = config['methods']['rf']['test']['mem'],
                time = config['methods']['rf']['test']['time']
        threads:
                config['methods']['rf']['test']['threads']
        benchmark:
                CELLTYPING_OUT + '{sample}/rf/predicted_labels.{model}.benchmark'
        shell:
                config['methods']['rf']['test']['call'] + ' --SCE {input.sce_in} --rf_model {input.trained_model} --outputFile {output.pred_labels}'

# method svm, testing
rule svm_test:
        input:
                sce_in =  INPUTDIR + '{sample}.RDS',
                model = CELLTYPING_IN + '{model}.svm_model.RDS'
        output:
                pred_labels = CELLTYPING_OUT + '{sample}/svm/predicted_labels.{model}.csv'
        params:
                lsfoutfile = CELLTYPING_OUT + '{sample}/svm/predicted_labels.{model}.lsfout.log',
                lsferrfile = CELLTYPING_OUT + '{sample}/svm/predicted_labels.{model}.lsferr.log',
                scratch = config['methods']['svm']['test']['scratch'],
                mem = config['methods']['svm']['test']['mem'],
                time = config['methods']['svm']['test']['time']
        threads:
                config['methods']['svm']['test']['threads']
        benchmark:
                CELLTYPING_OUT + '{sample}/svm/predicted_labels.{model}.benchmark'
        shell:
                config['methods']['svm']['test']['call'] + ' --SCE {input.sce_in} --svm_model {input.model} --outputFile {output.pred_labels}'

# method fcNN, testing
rule fcnn_test:
        input:
                sce_in =  INPUTDIR + '{sample}.RDS',
                model = CELLTYPING_IN + '{model}.fcnn_model.h5'
        output:
                pred_labels = CELLTYPING_OUT + '{sample}/fcnn/predicted_labels.{model}.csv'
        params:
                lsfoutfile = CELLTYPING_OUT + '{sample}/fcnn/predicted_labels.{model}.lsfout.log',
                lsferrfile = CELLTYPING_OUT + '{sample}/fcnn/predicted_labels.{model}.lsferr.log',
                scratch = config['methods']['fcnn']['test']['scratch'],
                mem = config['methods']['fcnn']['test']['mem'],
                time = config['methods']['fcnn']['test']['time']
        threads: 
                config['methods']['fcnn']['test']['threads']
        benchmark:
                CELLTYPING_OUT + '{sample}/fcnn/predicted_labels.{model}.benchmark'
        shell:
                config['methods']['fcnn']['test']['call'] + ' --SCE {input.sce_in} --nn_model {input.model} --outputFile {output.pred_labels}'

if not 'CDS_OUT' in globals():
    CDS_OUT = OUTDIR + "marker_genes/"

# method garnett, training, crossvalidation
rule garnett_test:
        input:
                cds_in =  CDS_OUT + '{sample}.monocle3_cds.RDS',
                model = CELLTYPING_IN + '{model}.garnett_model.RDS'
        output:
                pred_labels = CELLTYPING_OUT + '{sample}/garnett/predicted_labels.{model}.csv'
        params:
                lsfoutfile = CELLTYPING_OUT + '{sample}/garnett/predicted_labels.{model}.lsfout.log',
                lsferrfile = CELLTYPING_OUT + '{sample}/garnett/predicted_labels.{model}.lsferr.log',
                scratch = config['methods']['garnett']['test']['scratch'],
                mem = config['methods']['garnett']['test']['mem'],
                time = config['methods']['garnett']['test']['time']
        threads:
                config['methods']['garnett']['test']['threads']
        benchmark:
                CELLTYPING_OUT + '{sample}/garnett/predicted_labels.{model}.benchmark'
        shell:
                config['methods']['garnett']['test']['call'] + ' --cds {input.cds_in} --classifier {input.model} --outputFile {output.pred_labels}'

if not 'SCROSHI_IN' in globals():
    SCROSHI_IN = OUTDIR + "marker_genes/"

# prep config file based on present cell types, for scROSHI
rule prepare_scROSHI_config_file:
        input:
                celltype_gmx = SCROSHI_IN + '{sample}.gmx'
        output:
                celltype_config = SCROSHI_IN + '{sample}.celltype_config.tsv'
        params:
                lsfoutfile = SCROSHI_IN + '{sample}.prep_celltype_config.lsfout.log',
                lsferrfile = SCROSHI_IN + '{sample}.prep_celltype_config.lsferr.log',
                scratch = config['methods']['prep_scROSHI_config']['scratch'],
                mem = config['methods']['prep_scROSHI_config']['mem'],
                time = config['methods']['prep_scROSHI_config']['time'],
                celltype_config_gt = config['resources']['celltype_config_gt']
        #threads:
        #        config['methods']['prep_scROSHI_config']['threads']
        benchmark:
                SCROSHI_IN + '{sample}.prep_celltype_config.benchmark'
        shell:
                config['methods']['prep_scROSHI_config']['call'] + ' --gmx {input.celltype_gmx} --config_groundtruth {params.celltype_config_gt} --outfile {output.celltype_config}'

# scROSHI, testing
rule scROSHI:
        input:
                sce =  INPUTDIR + '{sample}.RDS', 
                celltype_gmx = SCROSHI_IN + '{sample}.gmx',
		celltype_config = SCROSHI_IN + '{sample}.celltype_config.tsv'
        output:
                pred_labels = CELLTYPING_OUT + '{sample}/scROSHI/predicted_labels.csv'
        params:
                lsfoutfile = CELLTYPING_OUT + '{sample}/scROSHI/predicted_labels.lsfout.log',
                lsferrfile = CELLTYPING_OUT + '{sample}/scROSHI/predicted_labels.lsferr.log',
                scratch = config['methods']['scROSHI']['scratch'],
                mem = config['methods']['scROSHI']['mem'],
                time = config['methods']['scROSHI']['time'],
                params = config['methods']['scROSHI']['params'],
		outputDirec = CELLTYPING_OUT,
		sampleName = '{sample}'
        threads:
                config['methods']['scROSHI']['threads']
        benchmark:
                CELLTYPING_OUT + '{sample}/scROSHI/predicted_labels.benchmark'
        shell:
                config['methods']['scROSHI']['call'] + ' --SCE {input.sce} --celltype_lists {input.celltype_gmx} --celltype_config {input.celltype_config} --outputDirec {params.outputDirec}{params.sampleName}/scROSHI/ --sampleName {params.sampleName}'

rule evaluate_classifier:
        input:
                sce_in = INPUTDIR + "{sample}.RDS",
                predicted_labels = CELLTYPING_OUT + "{sample}/{method}/predicted_labels.{model}.csv"
        output:
                CELLTYPING_OUT + '{sample}/{method}/{sample}.{model}.confusion_{method}.csv',
                CELLTYPING_OUT + '{sample}/{method}/{sample}.{model}.f1-score_{method}.csv',
                CELLTYPING_OUT + '{sample}/{method}/{sample}.{model}.population-size_{method}.csv',
                CELLTYPING_OUT + '{sample}/{method}/{sample}.{model}.summary_{method}.csv'
        params:
                lsfoutfile = CELLTYPING_OUT + '{sample}/{method}/{sample}.{model}.cv_summary.lsfout.log',
                lsferrfile = CELLTYPING_OUT + '{sample}/{method}/{sample}.{model}.cv_summary.lsferr.log',
                scratch = config['statistics']['cv_summary']['scratch'],
                mem = config['statistics']['cv_summary']['mem'],
                time = config['statistics']['cv_summary']['time'],
                params = config['statistics']['cv_summary']['params'],
                output_dir = CELLTYPING_OUT + "{sample}/{method}/",
                sample_name = "{sample}.{model}",
                method_name = "{method}"
        threads:
                config['statistics']['cv_summary']['threads']
        benchmark:
                CELLTYPING_OUT + '{sample}/{method}/{sample}.{model}.{method}.cv_summary.benchmark'
        shell:
                config['statistics']['calc_accuracy']['call'] + ' --sce {input.sce_in} --predicted_labels_path {input.predicted_labels} --output_dir {params.output_dir} --sample_name {params.sample_name} --method_name {params.method_name}'
