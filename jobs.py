from itertools import product
import re
import os

nodes = [1, 2, 4, 8, 16]
mpi_per_node = [1, 2, 4, 8]
flag_diag_arr = ['.true.', '.false.']
flag_mpi_arr = [0, 1, 2, 3, 4]

def replace_fun(line, **kwa):
    newLine = re.sub(r'<OMP_NUM_THREADS>', str(kwa['omp_num_threads']), line)
    newLine = re.sub(r'<MPI_PROCS>', str(kwa['mpi_procs']), newLine)
    newLine = re.sub(r'<PROCS_PER_NODE>', str(kwa['procs_per_node']), newLine)
    newLine = re.sub(r'<FLAG_MPI>', str(kwa['flag_mpi']), newLine)
    newLine = re.sub(r'<FLAG_DIAG>', str(kwa['flag_diag']), newLine)
    newLine = re.sub(r'<JOB_NAME>', str(kwa['job_name']), newLine)
    newLine = re.sub(r'<NUM_NODES>', str(kwa['num_nodes']), newLine)
    return newLine

if not os.path.isdir('jobs'):
    os.mkdir('jobs')
counter = 0
for i, j, diag, mpi in product(nodes, mpi_per_node, flag_diag_arr, flag_mpi_arr):
    for k in range(1, 12//j+1):
        counter += 1
        string = ('{i} nodes, {j} mpi/node, ({ij} total), {k} OpenMP/mpi ({jk}/node, {ijk} total)' +
                  ', flag_diag={diag}, flag_mpi={mpi} ({counter})')
        print(string.format(i= i, j= j, k= k, ij=i*j, jk=j*k, ijk=i*j*k, diag=diag,
                            mpi=mpi, counter=counter))

        with open("job-hybrid.slurm", "r") as f:
            lines = f.readlines()

            if diag == '.true.':
                diagstr = 'diag'
            else:
                diagstr = ''

            jobName = 'i{}j{}k{}{}mpi{}'.format(i,j,k,diagstr,mpi)

        with open("jobs/"+jobName, 'w') as f:
            for line in lines:
                s = replace_fun(line,
                                num_nodes=i,
                                omp_num_threads=k,
                                mpi_procs=i*j,
                                procs_per_node=j,
                                flag_mpi=mpi,
                                flag_diag=diag,
                                job_name=jobName)
                f.write(s)

