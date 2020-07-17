def get_sample_names():
    output = [] 
    if output == []:
        if not 'SAMPLEMAPPING' in globals():
            return ['NOMAPPINGFILE']
        try:
            open(SAMPLEMAPPING, "r")
        except IOError:
            return ['NOMAPPINGFILE']
        sampleMap = dict()
        with open(SAMPLEMAPPING, "r") as f:
            for line in f:
                if line.strip() != "":
                    lineSplit = line.strip().split()
                    sample = lineSplit[1]
                    if not (sample in output):
                        output.append(sample)
    return output

def get_model_names():
    output = []
    if output == []:
        if not 'METHOD_TEST_IN' in globals():
            METHOD_TEST_IN = OUTDIR + "trained_model/"
        try:
            dir_files = os.listdir(METHOD_TEST_IN)
        except IOError:
            raise("Trained model dir not detected. Please make sure that training pipeline was run first.")
        else:
            for each_file in dir_files:
                if "model" in each_file:
                    output.append(each_file)
    return output
