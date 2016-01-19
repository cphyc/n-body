#!/bin/python
from pathlib import Path
import pandas as pd

workPath = Path('.') #Path('/mnt/lnec/travail/bpagani/')
callerPath = Path('.') # Path('/obs/bpagani/n-body')
def parseParams(f, runData):
    line = f.readlines()[0].split()
    runData['nodes'] = int(line[0])
    runData['mpi_per_node'] = int(line[2])
    runData['mpi_procs'] = runData['nodes']*runData['mpi_per_node']
    runData['omp_per_mpi'] = int(line[6].replace('(', ''))

def parseSimulLog(f, runData):
    line = f.readlines()[-1].split()
    runData['total_time'] = float(line[1].replace('s', ''))
    runData['user_time'] = float(line[3].replace('s', ''))
    runData['system_time'] = float(line[5].replace('s', ''))
    runData['cpu'] = int(line[7].replace('%', ''))

def parseOutput(outputFile, runData):
    line = outputFile.readline()
    while '# Simulation parameters' not in line:
        line = outputFile.readline()

    nl = lambda : outputFile.readline().replace('\n', '')

    # parse the file
    runData['dt'] = float(nl().split(':')[-1])
    runData['npoints'] = int(nl().split(':')[-1])
    runData['maxtime'] = float(nl().split(':')[-1])
    runData['maxiter'] = int(nl().split(':')[-1])
    runData['nprocs'] = int(nl().split(':')[-1])
    runData['N'] = int(nl().split(':')[-1])
    diagRaw = nl().split(':')[-1]
    if 'T' in diagRaw:
        runData['diag'] = True
    else:
        runData['diag'] = False

    memoryRaw = nl().split(':')[-1]
    if 'T' in memoryRaw:
        runData['memory'] = True
    else:
        runData['memory'] = False

    runData['memoryFactor'] = int(nl().split(':')[-1])

def analyseRun(path):
    run = str(path).split('.')[-1]
    simulLogFile = (path / 'simul.log').open()
    paramsFile = (path / 'params.log').open()

    outputFile = (callerPath / ('slurm-'+str(run)+'.out')).open()

    runData = {}
    parseParams(paramsFile, runData)
    parseSimulLog(simulLogFile, runData)
    parseOutput(outputFile, runData)

    runData['run'] = run
    return runData

def keepRun(run):
    return run > 1

if __name__ == '__main__':
    data = pd.DataFrame()
    counter = 0

    # iterate all paths matching 'run.*'
    # only keep directories and runs > MIN_RUN_NUMBER
    def keepDir(_dir):
        b = _dir.is_dir()

        dirStr = str(_dir)
        run = int(dirStr.split('.')[-1])

        b = b and keepRun(run)
        return b

    dirList = [d for d in Path('.').glob('run.*') if keepDir(d)]
    for _dir in dirList:
        # read data
        dataDict = analyseRun(_dir)

        # copy data in panda data frame
        for key in dataDict.keys():
            data.loc[counter,key] = dataDict[key]

        counter += 1

    # do something useful with it, or not
    print(data)
