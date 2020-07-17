# Collection of rules necessary to run the celltype benchmark testing pipeline

if not 'METHOD_TEST_IN' in globals():
    METHOD_TEST_IN = OUTDIR + "trained_model/"
if not 'METHOD_TEST_OUT' in globals():
    METHOD_TEST_OUT = OUTDIR + 'test_results/'

# helper function to retrieve all trained models (saved as .RDS)
def get_trained_models(wildcards):
        return expand(METHOD_TRAIN_OUT + wildcards.sample + '_{method}_model.RDS', method = SELECTED_METHODS) #needs to differentiate model trained on ALL cells

# method rf, testing
rule rf_test:
        input:
                sce_in =  INPUTDIR + '{sample}.genes_cells_filtered.corrected.RDS',
                model = METHOD_TEST_IN + '{sample}.rf.model.RDS'
        output:
                pred_labels = METHOD_TEST_OUT + '{sample}.rf.predicted_labels.csv'
        params:
                lsfoutfile = METHOD_TEST_OUT + '{sample}.rf_test.lsfout.log',
                lsferrfile = METHOD_TEST_OUT + '{sample}.rf_test.lsferr.log',
                scratch = config['methods']['rf']['test']['scratch'],
                mem = config['methods']['rf']['test']['mem'],
                time = config['methods']['rf']['test']['time'],
                params = config['methods']['rf']['test']['params'],
                output_dir = METHOD_TEST_OUT,
                sample_name = "{sample}",
                method_name = "rf"
        threads:
                config['methods']['rf']['test']['threads']
        benchmark:
                METHOD_TEST_OUT + '{sample}.rf_test.benchmark'
        shell:
                config['methods']['rf']['test']['call'] + ' --SCE {input.sce_in} --rf_model {input.model} --outputDirec {params.output_dir} --sampleName {params.sample_name} -- method {params.method_name}'

# method svm, testing
rule svm_test:
        input:
                sce_in =  INPUTDIR + '{sample}.genes_cells_filtered.corrected.RDS',
                model = METHOD_TEST_IN + '{sample}.svm.model.RDS'
        output:
                pred_labels = METHOD_TEST_OUT + '{sample}.svm.predicted_labels.csv'
        params:
                lsfoutfile = METHOD_TEST_OUT + '{sample}.svm_test.lsfout.log',
                lsferrfile = METHOD_TEST_OUT + '{sample}.svm_test.lsferr.log',
                scratch = config['methods']['svm']['test']['scratch'],
                mem = config['methods']['svm']['test']['mem'],
                time = config['methods']['svm']['test']['time'],
                params = config['methods']['svm']['test']['params'],
                output_dir = METHOD_TEST_OUT,
                sample_name = "{sample}",
                method_name = "svm"
        threads:
                config['methods']['svm']['test']['threads']
        benchmark:
                METHOD_TEST_OUT + '{sample}.svm_test.benchmark'
        shell:
                config['methods']['svm']['test']['call'] + ' --SCE {input.sce_in} --rf_model {input.model} --outputDirec {params.output_dir} --sampleName {params.sample_name} --method {params.method_name}'

# method fcNN, testing
rule fcnn_test:
        input:
                sce_in =  INPUTDIR + '{sample}.genes_cells_filtered.corrected.RDS',
                model = METHOD_TEST_IN + '{sample}.fcnn.model.RDS'
        output:
                pred_labels = METHOD_TEST_OUT + '{sample}.fcnn.predicted_labels.csv'
        params:
                lsfoutfile = METHOD_TEST_OUT + '{sample}.fcnn_test.lsfout.log',
                lsferrfile = METHOD_TEST_OUT + '{sample}.fcnn_test.lsferr.log',
                scratch = config['methods']['fcnn']['test']['scratch'],
                mem = config['methods']['fcnn']['test']['mem'],
                time = config['methods']['fcnn']['test']['time'],
                params = config['methods']['fcnn']['test']['params'],
                output_dir = METHOD_TEST_OUT,
                sample_name = "{sample}",
                method_name = "fcnn"
        threads: 
                config['methods']['fcnn']['test']['threads']
        benchmark:
                METHOD_TEST_OUT + '{sample}.fcnn_test.benchmark'
        shell:
                config['methods']['fcnn']['test']['call'] + ' --SCE {input.sce_in} --rf_model {input.model} --outputDirec {params.output_dir} --sampleName {params.sample_name} -- method {params.method_name}'

if not 'CDS_OUT' in globals():
    CDS_OUT = OUTDIR + "marker_genes/"

# method garnett, training, crossvalidation
rule garnett_test:
        input:
                cds_in =  CDS_OUT + '{sample}.monocle3_cds.RDS',
                model = METHOD_TEST_IN + '{sample}.garnett.model.RDS'
        output:
                pred_labels = METHOD_TEST_OUT + '{sample}.garnett.predicted_labels.csv'
        params:
                lsfoutfile = METHOD_TEST_OUT + '{sample}.garnett_test.lsfout.log',
                lsferrfile = METHOD_TEST_OUT + '{sample}garnett_test.lsferr.log',
                scratch = config['methods']['garnett']['test']['scratch'],
                mem = config['methods']['garnett']['test']['mem'],
                time = config['methods']['garnett']['test']['time'],
                params = config['methods']['garnett']['test']['params'],
                output_dir = METHOD_TEST_OUT,
                sample_name = "{sample}"
        threads:
                config['methods']['garnett']['test']['threads']
        conda:
                config['methods']['garnett']['conda']
        benchmark:
                METHOD_TEST_OUT + '{sample}.garnett_test.benchmark'
        shell:
                config['methods']['garnett']['test']['call'] + ' --cds {input.cds_in} --model {input.model} --output_dir {params.output_dir} --pred_labels_out {output.pred_labels} --sample_name {params.sample_name}'

# Michael's method, testing
rule nexus_test:
        input:
                sce_in =  INPUTDIR + '{sample}.genes_cells_filtered.corrected.RDS', 
                gmx_in = METHOD_TEST_IN + '{sample}.gmx'
        output:
                pred_labels = METHOD_TEST_OUT + '{sample}.nexus.predicted_labels.csv'
        params:
                lsfoutfile = METHOD_TEST_OUT + '{sample}.nexus_test.lsfout.log',
                lsferrfile = METHOD_TEST_OUT + '{sample}.nexus_test.lsferr.log',
                scratch = config['methods']['nexus']['test']['scratch'],
                mem = config['methods']['nexus']['test']['mem'],
                time = config['methods']['nexus']['test']['time'],
                params = config['methods']['nexus']['test']['params'],
                output_dir = METHOD_TEST_OUT,
                sample_name = "{sample}",
                method_name = "nexus"
        threads:
                config['methods']['nexus']['test']['threads']
        benchmark:
                METHOD_TEST_OUT + '{sample}.nexus_test.benchmark'
        shell:
                config['methods']['nexus']['test']['call'] + ' --SCE {input.sce_in} --marker_list {input.gmx_in} --outputDirec {params.output_dir} --sampleName {params.sample_name} -- method {params.method_name}'

# evaluate models & calculate accuracy statistics (F1-score)
rule evaluate_classifier:
        input:
                sce_in = INPUTDIR + "{sample}.genes_cells_filtered.corrected.RDS",
                predicted_labels = METHOD_TEST_OUT + "{sample}.{method}.predicted_labels.csv"
        output:
                METHOD_TEST_OUT + '{sample}.confusion_{method}.csv',
                METHOD_TEST_OUT + '{sample}.f1-score_{method}.csv',
                METHOD_TEST_OUT + '{sample}.population-size_{method}.csv',
                METHOD_TEST_OUT + '{sample}.summary_{method}.csv'
        params:
                lsfoutfile = METHOD_TEST_OUT + '{sample}.cv_summary.lsfout.log',
                lsferrfile = METHOD_TEST_OUT + '{sample}.cv_summary.lsferr.log',
                scratch = config['statistics']['cv_summary']['scratch'],
                mem = config['statistics']['cv_summary']['mem'],
                time = config['statistics']['cv_summary']['time'],
                params = config['statistics']['cv_summary']['params'],
                output_dir = METHOD_TEST_OUT,
                sample_name = "{sample}",
                method_name = "{method}"
        threads:
                config['statistics']['cv_summary']['threads']
        benchmark:
                METHOD_TEST_OUT + '{sample}.cv_summary.benchmark'
        shell:
                config['statistics']['cv_summary']['call'] + ' --sce {input.sce_in} --predicted_labels_path {input.predicted_labels} --output_dir {params.output_dir} --sample_name {params.sample_name} --method_name {params.method_name}'
