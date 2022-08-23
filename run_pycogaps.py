''' 
this script reads parameters from the command line to run CoGAPS
supports integration with genepattern notebook
'''

import sys
import csv
import ast

if __name__ == '__main__':
    from PyCoGAPS.config import *
    from PyCoGAPS.parameters import *
    from PyCoGAPS.pycogaps_main import CoGAPS
    import pickle
    import argparse
    
    print("This vignette was built using pycogaps version", getVersion())

    ''' 
    command line args which are all parameters to CoGAPS
        - only --path arg is required
        - all other args are optional, have default values
    '''
    parser = argparse.ArgumentParser()

    parser.add_argument('--path', type=str, required=True)
    parser.add_argument('--resultFile', type=str, default='result.pkl')
    
    # standard params
    parser.add_argument('--nPatterns', type=int, default=3)
    parser.add_argument('--nIterations', type=int, default=1000)
    parser.add_argument('--seed', type=int, default=0)
    parser.add_argument('--useSparseOptimization', type=str, default="False", choices=["False","True"])
    
    # run params
    parser.add_argument('--nThreads', type=int, default=1)
    parser.add_argument('--messages', type=str, default="True", choices=["False","True"])
    parser.add_argument('--outputFrequency', type=int, default=500)
    parser.add_argument('--uncertainty', type=str, default=None) # read in as file, matrix
    parser.add_argument('--checkpointOutFile', type=str, default='gaps_checkpoint.out')
    parser.add_argument('--checkpointInFile', type=str, default="")
    parser.add_argument('--transposeData', type=str, default="False", choices=["False","True"])
    parser.add_argument('--workerID', type=int, default=1)
    parser.add_argument('--asynchronousUpdates', type=str, default=False)
    parser.add_argument('--nSnapshots', type=int, default=0)
    parser.add_argument('--snapshotPhase', type=str, default='sampling', choices=['sampling', 'equilibration', 'all'])
    
    # sparsity params
    parser.add_argument('--alphaA', type=float, default=0.01)
    parser.add_argument('--alphaP', type=float, default=0.01)
    parser.add_argument('--maxGibbsMassA', type=int, default=100)
    parser.add_argument('--maxGibbsMassP', type=int, default=100)
    
    # distributed params
    parser.add_argument('--distributed', type=str, default=None)
    parser.add_argument('--nSets', type=int, default=4)
    parser.add_argument('--cut', type=int, default=None)
    parser.add_argument('--minNS', type=int, default=None)
    parser.add_argument('--maxNS', type=int, default=None)
    parser.add_argument('--explicitSets', type=str, default=None) # read in as file, list
    parser.add_argument('--samplingAnnotation', type=str, default=None) # read in as file, list
    parser.add_argument('--samplingWeight', type=str, default=None) # read in as file, dictionary
    
    # additional params
    parser.add_argument('--subsetIndices', type=str, default=None) # read in as file, list
    parser.add_argument('--subsetDim', type=int, default=0, choices=[0,1])
    parser.add_argument('--geneNames', type=str, default=None) # read in as file, list
    parser.add_argument('--sampleNames', type=str, default=None) # read in as file, list
    parser.add_argument('--fixedPatterns', type=str, default=None) # read in as file, matrix
    parser.add_argument('--whichMatrixFixed', type=str, default=None, choices=['A', 'P'])
    parser.add_argument('--takePumpSamples', type=str, default="False", choices=["False","True"])
    parser.add_argument('--hdfKey', type=str, default=None)
    parser.add_argument('--hdfRowKey', type=str, default=None)
    parser.add_argument('--hdfColKey', type=str, default=None)
    
    initial_params = ["path", "resultFile"]
    
    standard_params = ["nPatterns", "nIterations", "seed", "useSparseOptimization"]

    run_params = ["nThreads", "messages", "outputFrequency", "uncertainty", "checkpointOutFile", "checkpointInterval", 
                "checkpointInFile", "transposeData", "workerID", "asynchronousUpdates", 
                "nSnapshots", "snapshotPhase"]

    sparsity_params = ["alphaA", "alphaP", "maxGibbsMassA", "maxGibbsMassP"]
    
    distributed_params = ["distributed", "nSets", "cut", "minNS", "maxNS", 
                "explicitSets", "samplingAnnotation", "samplingWeight"]
   
    additional_params = ["subsetIndices", "subsetDim", "geneNames", "sampleNames",   
                "fixedPatterns", "whichMatrixFixed", "takePumpSamples", 
                "hdfKey", "hdfRowKey", "hdfColKey"]

    bool_params = ["messages", "useSparseOptimization", "asynchronousUpdates", "takePumpSamples"]
    
    '''
    parse all args and set as parameters for CoGAPS
    '''
    args = parser.parse_args()

    data_path = args.path

    def to_bool(val):
        if val == "True":
            return True
        else:
            return False
    
    args.transposeData = to_bool(args.transposeData)
    
    params = CoParams(path=data_path, transposeData=args.transposeData, 
                      hdfKey=args.hdfKey, hdfRowKey=args.hdfRowKey,
                      hdfColKey=args.hdfColKey)


    prm_dict = vars(args)

    ## read in file arguments as vectors/lists/matrices ##
    list_params = ["explicitSets", "samplingAnnotation", "samplingWeight", "subsetIndices", "geneNames", "sampleNames"]
    def file_to_type(k, file_path):
        if file_path is None:
            return file_path
        if k == "samplingWeight": # read this as txt -> dictionary
            with open(file_path) as f:
                data = f.read()
                return ast.literal_eval(data)
        with open(file_path, newline='') as f: # read all others as csv -> list
            reader = csv.reader(f)
            return list(reader)[0]

    for k,v in prm_dict.items():
        if ((k not in initial_params) and (k not in distributed_params) and (k not in ("fixedPatterns", "uncertainty"))):
            if k in list_params:
                v = file_to_type(k, v)
            if k in bool_params:
                v = to_bool(v)
                print(type(v))
            print(k)
            setParam(params, k, v)
            print('done')
    
     # set fixed patterns from additional params
    if args.fixedPatterns is not None:
        params.setFixedPatterns(fixedPatterns=args.fixedPatterns, whichMatrixFixed=args.whichMatrixFixed)

    # set distributed parameters
    setParam(params, 'distributed', args.distributed)
    if args.distributed is not None:
        params.setAnnotationWeights(annotation=args.samplingAnnotation, weight=args.samplingWeight)
        params.setDistributedParams(nSets=args.nSets, cut=args.cut, minNS=args.minNS, maxNS=args.maxNS)

    '''
    run CoGAPS, save result
    '''
    result = CoGAPS(data_path, params, transposeData=args.transposeData, uncertainty=args.uncertainty)

    # save CoGAPS result
    print("Pickling...", end='\r')
    pickle.dump(result, open(args.resultFile, "wb"))
    print("Pickling complete!")
