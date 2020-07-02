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
