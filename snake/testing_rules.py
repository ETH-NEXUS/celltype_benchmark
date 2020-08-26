# Collection of rules necessary to run the celltype benchmark testing pipeline

if not 'CELLTYPING_IN' in globals():
    CELLTYPING_IN = OUTDIR + "trained_model/"
if not 'CELLTYPING_OUT' in globals():
    CELLTYPING_OUT = OUTDIR + 'test_results/'


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
                time = config['methods']['rf']['test']['time'],
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
                time = config['methods']['svm']['test']['time'],
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
                time = config['methods']['fcnn']['test']['time'],
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
                pred_labels = CELLTYPING_OUT + '{sample}/garnett/predicted_labels_{model}.csv'
        params:
                lsfoutfile = CELLTYPING_OUT + '{sample}/garnett/predicted_labels_{model}.lsfout.log',
                lsferrfile = CELLTYPING_OUT + '{sample}/garnett/predicted_labels_{model}.lsferr.log',
                scratch = config['methods']['garnett']['test']['scratch'],
                mem = config['methods']['garnett']['test']['mem'],
                time = config['methods']['garnett']['test']['time'],
        threads:
                config['methods']['garnett']['test']['threads']
        conda:
                config['methods']['garnett']['conda']
        benchmark:
                CELLTYPING_OUT + '{sample}/garnett/predicted_labels_{model}.benchmark'
        shell:
                config['methods']['garnett']['test']['call'] + ' --cds {input.cds_in} --classifier {input.model} --outputFile {output.pred_labels}'

if not 'SCROSI_IN' in globals():
    SCROSI_IN = OUTDIR + "marker_genes/"
# scROSI, testing
rule scROSI:
        input:
                sce_in =  INPUTDIR + '{sample}.RDS', 
                gmx_in = SCROSI_IN + '{model}.gmx'
        output:
                pred_labels = CELLTYPING_OUT + '{sample}/scROSI/predicted_labels_{model}.csv'
        params:
                lsfoutfile = CELLTYPING_OUT + '{sample}.nexus_test.lsfout.log',
                lsferrfile = CELLTYPING_OUT + '{sample}.nexus_test.lsferr.log',
                scratch = config['methods']['nexus']['test']['scratch'],
                mem = config['methods']['nexus']['test']['mem'],
                time = config['methods']['nexus']['test']['time'],
                params = config['methods']['nexus']['test']['params'],
                output_dir = CELLTYPING_OUT,
                sample_name = "{sample}",
                method_name = "nexus"
        threads:
                config['methods']['nexus']['test']['threads']
        benchmark:
                CELLTYPING_OUT + '{sample}.nexus_test.benchmark'
        shell:
                config['methods']['nexus']['test']['call'] + ' --SCE {input.sce_in} --marker_list {input.gmx_in} --outputDirec {params.output_dir} --sampleName {params.sample_name} -- method {params.method_name}'

# evaluate models & calculate accuracy statistics (F1-score)
rule evaluate_classifier:
        input:
                sce_in = INPUTDIR + "{sample}.genes_cells_filtered.corrected.RDS",
                predicted_labels = CELLTYPING_OUT + "{sample}.{method}.predicted_labels.csv"
        output:
                CELLTYPING_OUT + '{sample}.confusion_{method}.csv',
                CELLTYPING_OUT + '{sample}.f1-score_{method}.csv',
                CELLTYPING_OUT + '{sample}.population-size_{method}.csv',
                CELLTYPING_OUT + '{sample}.summary_{method}.csv'
        params:
                lsfoutfile = CELLTYPING_OUT + '{sample}.cv_summary.lsfout.log',
                lsferrfile = CELLTYPING_OUT + '{sample}.cv_summary.lsferr.log',
                scratch = config['statistics']['cv_summary']['scratch'],
                mem = config['statistics']['cv_summary']['mem'],
                time = config['statistics']['cv_summary']['time'],
                params = config['statistics']['cv_summary']['params'],
                output_dir = CELLTYPING_OUT,
                sample_name = "{sample}",
                method_name = "{method}"
        threads:
                config['statistics']['cv_summary']['threads']
        benchmark:
                CELLTYPING_OUT + '{sample}.cv_summary.benchmark'
        shell:
                config['statistics']['cv_summary']['call'] + ' --sce {input.sce_in} --predicted_labels_path {input.predicted_labels} --output_dir {params.output_dir} --sample_name {params.sample_name} --method_name {params.method_name}'
