#!/usr/bin/env python3
""" Task manager for running the dask pipeline for running structure alignment, 
    as defined in the input list file.
    USAGE: 
        python3 dask_taskmgr.py [-h] --timings-file TIMINGS_FILE.csv \
                                     --tskmgr-log-file TSKMGR.log \
                                     --alignments-pdb-list-file ALN_PATH \
                                     --scheduler-file SCHEDULER_FILE \
                                     --nIterations 1
                                     
    INPUT: 
        -h, --help      show this help message and exit
        --timings-file TIMINGS_FILE.csv, -ts TIMINGS_FILE.csv 
                        CSV file for task processing timings
        --tskmgr-log-file TSKMGR.log, -log TSKMGR.log 
                        path for a file within which logging output for the 
                        workflow will be written
        --alignments-pdb-list-file ALN_PATH, -inp ALN_PATH
                        path to a list file that contains the paths to the pair
                        of structure models to be aligned; third element in each
                        line is the alignment method to be used for the calc,
                        mirroring the options for USalign -mm parameter. 5 and
                        6 for fNS and sNS, respectively.
        --scheduler-file SCHEDULER_FILE, -s SCHEDULER_FILE
                        dask scheduler file; optional, default = None
        --nIterations, 1, -niters 1 
                        integer, number of iterations of each alignment to 
                        perform, default = 1
"""

import sys
import time
import io
import re
import argparse
import platform
import logging
import itertools
import numpy as np
import csv
from uuid import uuid4
from pathlib import Path

import subprocess
from subprocess import CalledProcessError

import dask
import dask.config
from distributed import Client, Worker, as_completed, get_worker

#######################################
### LOGGING FUNCTIONS
#######################################

def append_timings(csv_writer, file_object, hostname, worker_id, start_time, 
                   stop_time, AlnID):
    """ append the task timings to the CSV timings file
    :param csv_writer: CSV to which to append timings
    :param hostname: on which the processing took place
    :param worker_id: of the dask worker that did the processing
    :param start_time: start time in *NIX epoch seconds
    :param stop_time: stop time in same units
    :param nSuccesses: str or int identifying the alignment
    """
    csv_writer.writerow({'hostname'   : hostname,
                         'worker_id'  : worker_id,
                         'start_time' : start_time,
                         'stop_time'  : stop_time,
                         'AlnID' : AlnID})
    file_object.flush()


def setup_logger(name, log_file, level=logging.INFO):
    """To setup as many loggers as you want"""
    formatter = logging.Formatter('%(asctime)s    %(levelname)s       %(message)s')
    handler = logging.FileHandler(log_file)
    handler.setFormatter(formatter)

    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.addHandler(handler)

    return logger


def clean_logger(logger):
    """To cleanup the logger instances once we are done with them"""
    for handle in logger.handlers:
        handle.flush()
        handle.close()
        logger.removeHandler(handle)


#######################################
### PARSING FUNCTION
#######################################

def parse_usalign_file(file_object):
    """ Parse a USalign alignment file

    #BUG: regex to get struct1 and struct2 expect local or global paths that 
          include some slashes... so this propagates to the alignment workflow 
          where the list files need to include the paths... grumble

    INPUT:
    :param file_object: file object that can be read. many methods of creating 
                        such an object; ex: io.StringIO or "with open()" 
    RETURNS:
    :return: results, start_time, stop_time, return_value
    :return results: dictionary of quantitative results associated with the
                     alignment. 
            Keys:
                'struct1': string, path to the first structure in the alignment;
                           the mobile or query structure
                'struct2': string, path to the second structure in the 
                           alignment; the target structure
                'struct1_chainID': string, chain ID string associated with the 
                                   query structure
                'struct2_chainID': string, chain ID string associated with the 
                                   target structure
                'TMscore1': float, TMscore value normalized by the query 
                            structure's length
                'TMscore2': float, TMscore value normalized by the target 
                            structure's length
                'RMSD': float, root mean square deviation between the aligned 
                        atoms
                'SeqIDAli': float, sequence identity of aligned residues
                'Len1': float, the query structure's length
                'Len2': float, the target structure's length
                'LenAligned': float, the number of residues used in the 
                'dt': float, reported time for alignment; units: seconds
    """
    results = {}
    # file_object is a TextIOBase text stream created by open() or io.StringIO 
    # or some other way; here we just collect all lines
    lines = file_object.readlines()

    ### gather relevant lines common across all USalign outputs
    # query structure lines
    struct1_lines = [line.strip() for line in lines if 'Structure_1:' in line]
    # finds file name; finds the first instance of a match; takes the stem of
    # that path. 
    results['struct1'] = Path(re.search(r'[\w-]+\.', struct1_lines[0])[0]).stem
    #results['struct1'] = re.search(r'((?<!\w)(\.{1,2})?(?<!\/)(\/((\\\b)|[^ \b%\|:\n\"\\\/])+)+\/?)', struct1_lines[0])[0]    # finds absolute or relative paths; finds the first instance of a match. 
    # finds the chainID that was used in the alignment
    results['struct1_chainID'] = re.search(r'(?<=[:])\w+',struct1_lines[0])[0]
    # gathering quantitative metrics associated with TM-score, Length of query 
    # structure, and d0 value for the alignment
    results['TMscore1'], results['Len1'], results['d0_1'] = [float(elem) for elem in re.findall(r'(?<=[=\s](?=\d))\d*[\.]?\d*',struct1_lines[2])]
    
    # target structure lines
    struct2_lines = [line.strip() for line in lines if 'Structure_2:' in line]
    # finds file name; finds the first instance of a match; takes the stem of
    # that path. 
    results['struct2'] = Path(re.search(r'[\w-]+\.', struct2_lines[0])[0]).stem 
    #results['struct2'] = re.search(r'((?<!\w)(\.{1,2})?(?<!\/)(\/((\\\b)|[^ \b%\|:\n\"\\\/])+)+\/?)', struct2_lines[0])[0]    # finds absolute or relative paths; finds the first instance of a match. 
    # find the chainID that was used in the alignment
    results['struct2_chainID'] = re.search(r'(?<=[:])\w+',struct2_lines[0])[0] 
    # gathering quantitative metrics associated with TM-score, Length of target
    # structure, and d0 value for the alignment
    results['TMscore2'], results['Len2'], results['d0_2'] = [float(elem) for elem in re.findall(r'(?<=[=\s](?=\d))\d*[\.]?\d*',struct2_lines[2])] 
    
    # Alignment overview line
    aln_lines = [line.strip() for line in lines if 'Aligned length=' in line]
    # gathering quantitative metrics associated with Length of Alignment, RMSD,
    # and sequence ID for the aligned residues
    results['LenAligned'], results['RMSD'], results['SeqIDAli'] = [float(elem) for elem in re.findall(r'(?<=[=\s](?=\d))\d*[\.]?\d*',aln_lines[0])]

    # collect USalign reported timing
    time_lines = [line.strip() for line in lines if 'Total CPU time' in line]
    results['dt'] = float(re.findall(r'\d*[\.]?\d+',time_lines[0])[0])

    # collect translation and rotation lines if present
    trans_rot_array = np.array([line.split()[1:] for line in lines if line[0] in ['0','1','2']],dtype=float)
    if trans_rot_array.size > 0:
        results['trans_vector'] = trans_rot_array[:,0]
        results['rot_matrix'] = trans_rot_array[:,1:]

    return results


#######################################
### DASK RELATED FUNCTIONS
#######################################

def submit_subprocess_pipeline(alignments, nIterations):
    """
    alignment task function for the dask workflow
    INPUT:
        :param alignments: tuple of two strings, where the strings are paths to
                  the two proteins to be aligned nIterations times.
                  elem[0] is the path string pointing at the query protein
                  elem[1] is the path string associated with target protein,
                          being aligned to.
                  elem[2] is the USalign -mm parameter value to be used. 5 for 
                          fNS and 6 for sNS. 
        :param nIterations: int, number of iterations to perform.
    RETURNS:
        :return: platform information i.e. "hostname."
        :return: worker identification string.
        :return: start_time, end_time: floats, start and stop epoch times for 
                 total time of the iterations, including parsing time.
        :return: dictionary of results
                 
    """
    worker = get_worker()
    start_time = time.time()
    results_dict = {}

    for _iter in range(nIterations):
        # run the alignment script
        try:
            # 1) run USalign between the two structure files, using alignment 
            #    method, outputting alignment results to stdout, and saving the
            #    translocation/rotation matrix to file
            # 2) cat the trans/rot matrix file to stdout
            # 3) delete the trans/rot matrix file so we don't get a ton of files 
            temp_string = str(uuid4())
            completed_process = subprocess.run(
                    f'~/Apps/USalign/USalign {alignments[0]} {alignments[1]} -mm {alignments[2]} -outfmt 0 -m {temp_string}.dat; cat {temp_string}.dat; rm {temp_string}.dat',
                    shell=True,
                    capture_output=True,
                    check=True)
        # if the alignment script errors out, skip the parsing steps
        except Exception as e:
            print(f'submit_pipeline, aln, {_iter}, {alignments[0]}, {alignments[1]}, {alignments[2]}', e, file=sys.stderr, flush=True)
            continue

        # parse the alignment output
        try:
            # creates a file-like-object text stream
            stdout_string = completed_process.stdout.decode() 
            stdout_file = io.StringIO(stdout_string)
            # put this text stream through the parser
            # hard coding the aln_algo as '' because we are not interested in
            # the residue mapping at the present
            parsed_results = parse_usalign_file(stdout_file)
            # save alignment results as a key-value in the results_dict
            # key: string with format of "iter_0"
            # value: list, [TMscore1, TMscore2, RMSD, SeqIDAli, Len1, Len2,
            #               LenAligned, dt, trans_vector, rot_matrix]
            key = f'iter_{_iter}'
            results_dict[key] = [float(parsed_results['TMscore1']),
                                 float(parsed_results['TMscore2']),
                                 float(parsed_results['RMSD']),
                                 float(parsed_results['SeqIDAli']),
                                 float(parsed_results['Len1']),
                                 float(parsed_results['Len2']),
                                 float(parsed_results['LenAligned']),
                                 float(parsed_results['dt']),
                                 parsed_results['trans_vector'],
                                 parsed_results['rot_matrix']]
        except Exception as e:
            print(f'submit_pipeline, parsing, {_iter}, {alignments[0]}, {alignments[1]}', e, file=sys.stderr, flush=True)
            continue
        
    stop_time = time.time()
    return platform.node(), worker.id, alignments[0], alignments[1], alignments[2], start_time, stop_time, results_dict


#######################################
### MAIN
#######################################

if __name__ == '__main__':
    # read command line arguments.
    parser = argparse.ArgumentParser(description='Structural alignment task manager')
    # essential parameters
    parser.add_argument('--timings-file', 
                        '-ts', 
                        required=True, 
                        help='CSV file for protein processing timings')
    parser.add_argument('--tskmgr-log-file', 
                        '-log', 
                        required=True, 
                        help='string that will be used to store logging info for this run')
    parser.add_argument('--alignments-pdb-list-file', 
                        '-inp', 
                        required=True, 
                        help='list file where each line consists of the paths for the two protein models that are to be aligned')
    # optional parameters
    parser.add_argument('--scheduler-file', 
                        '-s', 
                        required=False,
                        default=None, 
                        help='dask scheduler file')
    parser.add_argument('--nIterations', 
                        '-niters', 
                        required=False,
                        default=1, 
                        help='number of iterations of each alignment to perform',
                        type=int)
    args = parser.parse_args()

    # set up the main logger file and list all relevant parameters.
    main_logger = setup_logger('tskmgr_logger',args.tskmgr_log_file)
    main_logger.info(f'Starting dask pipeline and setting up logging. Time: {time.time()}')
    main_logger.info(f'Scheduler file: {args.scheduler_file}')
    main_logger.info(f'Timing file: {args.timings_file}')
    main_logger.info(f'Alignments are listed in {args.alignments_pdb_list_file}.')
    main_logger.info(f'Each alignment will be performed {args.nIterations} times.')

    dask_parameter_string = ''
    for key, value in dask.config.config.items():
        dask_parameter_string += f"'{key}': '{value}'\n"
    dask_parameter_string += '################################################################################'
    main_logger.info(f'\n################################################################################\nDask parameters:\n{dask_parameter_string}')

    # start dask client.
    # if a scheduler file is input, then the CLI is being used, where the user 
    # has defined a scheduler and number of workers. Likely has also defined 
    # what resources each worker gets.
    if args.scheduler_file:
        client = Client(scheduler_file=args.scheduler_file,
                        timeout=120,
                        name='AlignmentTaskMgr')
    # no scheduler file, so assume that the default Client behavior is expected
    # (spin up a scheduler and worker within the Client call); since tasks in 
    # this workflow are known to only need 1 thread, we hard wire the client to
    # spin up the appropriate resourced workers.
    else:
        client = Client(timeout=120,
                        name='AlignmentTaskMgr',
                        threads_per_worker=1)
   
    # set up timing log file.
    main_logger.info(f'Opening the timing file.')
    timings_file = open(args.timings_file, 'w')
    timings_csv  = csv.DictWriter(timings_file,
                                  ['hostname',
                                   'worker_id',
                                   'start_time',
                                   'stop_time', 
                                   'AlnID'])
    timings_csv.writeheader()

    # parsing the alignments_pdb_list_file
    # alignment_list represents the list of structure pairs that will be aligned 
    main_logger.info(f'Reading the alignments file.')
    with open(args.alignments_pdb_list_file,'r') as structures_file:
        alignments_list = [line.split() for line in structures_file.readlines() if line[0] != '#']
    
    main_logger.info(f'A total of {len(alignments_list)*args.nIterations:,} alignments will be performed.')
    aln_futures = client.map(submit_subprocess_pipeline, 
                             alignments_list, 
                             nIterations = args.nIterations, 
                             pure=False)

    results_dict = {}
    aln_completed = as_completed(aln_futures)
    for finished_task in aln_completed:
        hostname, worker_id, query, target, aln_method, start_time, stop_time, results = finished_task.result()
        main_logger.info(f'-mm {aln_method} alignment between {query} and {target} has finished.')
        query = Path(query).stem
        target = Path(target).stem
        identifier = f'{query}|{target}|mm{aln_method}'
        append_timings(timings_csv, 
                       timings_file, 
                       hostname, 
                       worker_id,
                       start_time, 
                       stop_time,
                       identifier)
        results_dict[identifier] = results
   
    with open(f'USalign_subprocess_results.dat','w') as out_file:
        out_file.write('Query_Target_Method,TMscore1,TMscore2,RMSD,SeqIDAli,Len1,Len2,LenAligned),Avg Time (s),Stdev Time (s)\n')
        sorted_keys = list(results_dict.keys())
        for key in sorted_keys:
            results = results_dict[key]
            dts = [results[f'iter_{i}'][7] for i in range(args.nIterations)]
            assert results['iter_0'][0] == results['iter_1'][0]
            avg_time = np.mean(dts)
            std_time = np.std(dts)
            TMscore1 = results['iter_0'][0]
            TMscore2 = results['iter_0'][1]
            RMSD     = results['iter_0'][2]
            SeqIDAli = results['iter_0'][3]
            Len1     = results['iter_0'][4]
            Len2     = results['iter_0'][5]
            LenAlign = results['iter_0'][6]
            out_file.write(f'{key},{TMscore1},{TMscore2},{RMSD},{SeqIDAli},{Len1},{Len2},{LenAlign},{avg_time:.3f},{std_time:.3f}\n')

    timings_file.close()
    main_logger.info(f'Done. Shutting down the cluster. Time: {time.time()}')
    clean_logger(main_logger)

