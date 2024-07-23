#!/usr/bin/env python3
"""
    For parsing USalign alignment results

"""
import sys
import time
from pathlib import Path
import numpy as np
import pandas as pd
import re


def parse_usalign_file(file_object,alignment_type):
    """ Parse a USalign alignment file

    #BUG: regex to get struct1 and struct2 expect local or global paths that 
          include some slashes... so this propagates to the alignment workflow 
          where the list files need to include the paths... grumble

    INPUT:
    :param file_object: file object that can be read. many methods of creating 
                        such an object; ex: io.StringIO or "with open()" 
    :param alignment_type: string parameter to set format of alignment output to 
                           parse; currently accepted: CP, SNS, FNS, or empty
                           if left empty (i.e. ''), then just the quant metrics
                           will be gathered and returned, no mapping will be 
                           gathered
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
                              alignment
                'd0_1': float, the normalized distance for query structure; used 
                        in TMscore metric calculation
                'd0_2': float, the normalized distance for target structure; used 
                        in TMscore metric calculation
                'map_1_to_2': optional, dictionary, keys are query structure's 
                              ('1') residue indices (string) that map to a 
                              tuple with query residue name, target structure's
                              ('2') residue index (string), target residue name,
                              and alignment distance (float, only available for 
                              SNS and FNS alignments); if alignment_type is not 
                              in ['CP', 'SNS', 'FNS'], this mapping will not be 
                              created
    :return start_time: float, epoch time when the task began
    :return stop_time: float, epoch time when the task stopped
    :return return_value: integer, valued 0 if function returns all quant values
                          including the alignment mapping; otherwise valued 1. 
    """
    start_time = time.time()
    results = {}
    # file_object is a TextIOBase text stream created by open() or io.StringIO 
    # or some other way; here we just collect all lines
    lines = file_object.readlines()

    ### gather relevant lines common across all USalign outputs
    # query structure lines
    struct1_lines = [line.strip() for line in lines if 'Structure_1:' in line]
    results['struct1'] = re.search(r'((?<!\w)(\.{1,2})?(?<!\/)(\/((\\\b)|[^ \b%\|:\n\"\\\/])+)+\/?)', struct1_lines[0])[0]    # finds absolute or relative paths; finds the first instance of a match. 
    results['struct1_chainID'] = re.search(r'(?<=[:])\w+',struct1_lines[0])[0]    # finds the chainID that was used in the alignment
    results['TMscore1'], results['Len1'], results['d0_1'] = [float(elem) for elem in re.findall(r'(?<=[=\s](?=\d))\d*[\.]?\d*',struct1_lines[2])]    # gathering quantitative metrics associated with TM-score, Length of query structure, and d0 value for the alignment
    
    # target structure lines
    struct2_lines = [line.strip() for line in lines if 'Structure_2:' in line]
    results['struct2'] = re.search(r'((?<!\w)(\.{1,2})?(?<!\/)(\/((\\\b)|[^ \b%\|:\n\"\\\/])+)+\/?)', struct2_lines[0])[0]    # finds absolute or relative paths; finds the first instance of a match. 
    results['struct2_chainID'] = re.search(r'(?<=[:])\w+',struct2_lines[0])[0]    # finds the chainID that was used in the alignment
    results['TMscore2'], results['Len2'], results['d0_2'] = [float(elem) for elem in re.findall(r'(?<=[=\s](?=\d))\d*[\.]?\d*',struct2_lines[2])]    # gathering quantitative metrics associated with TM-score, Length of target structure, and d0 value for the alignment
    
    # Alignment overview line
    aln_lines = [line.strip() for line in lines if 'Aligned length=' in line]
    results['LenAligned'], results['RMSD'], results['SeqIDAli'] = [float(elem) for elem in re.findall(r'(?<=[=\s](?=\d))\d*[\.]?\d*',aln_lines[0])]    # gathering quantitative metrics associated with Length of Alignment, RMSD, and sequence ID for the aligned residues

    # collect translation and rotation lines if present
    trans_rot_array = np.array([line.split()[1:] for line in lines if line[0] in ['0','1','2']],dtype=float)
    if trans_rot_array.size > 0:
        results['trans_vector'] = trans_rot_array[:,0]
        results['rot_matrix'] = trans_rot_array[:,1:]

    if alignment_type.upper() == 'CP':
        # gather the alignment mapping, strip white space on both sides
        map_lines    = [line.strip() for line in lines if ' CA ' == line[:4]]

        # dictionary of tuples; keys are mobile residue index with values being a tuple of (target residue index, mobile resname, target resname)
        results['map_1_to_2'] = {}
        for line in map_lines:
            # line format: 'CA  LEU A 111 \t CA  GLN A  97'
            # awkward white spaces:      ^^
            # wanna keep: {'111': ('97','LEU','GLN')}
            temp = [elem.strip() for elem in line.split('\t')]
            results['map_1_to_2'].update({temp[0][-4:].strip(): (temp[1][-4:].strip(),temp[0][4:7],temp[1][4:7])})
    
    elif alignment_type.upper() in ['SNS','FNS']:
        # gather the alignment mapping, strip white space on both sides
        map_lines    = [line.strip() for line in lines if ' CA ' == line[:4]]
        
        # dictionary of tuples; keys are mobile residue index with values being a tuple of (target residue index, mobile resname, target resname, distance)
        results['map_1_to_2'] = {}
        for line in map_lines:
            # line format: 'CA  LEU A 111 \t CA  GLN A  97 \t    1.568'
            # awkward white spaces:      ^^               ^^
            # wanna keep: {'111': ('97','LEU','GLN',1.568)}
            temp = [elem.strip() for elem in line.split('\t')]
            results['map_1_to_2'].update({temp[0][-4:].strip(): (temp[1][-4:].strip(),temp[0][4:7],temp[1][4:7],float(temp[2]))})

    else:
        print(f"'{alignment_type}' not expected by parser function. Only returning alignment quantitative metrics. No mapping.", file=sys.stderr, flush=True)
        stop_time = time.time()
        return results, start_time, stop_time, 1

    stop_time = time.time()
    return results, start_time, stop_time, 0


def threshold_comparison(results_dict,ranking_metric_string,metric_cutoff):
    """ check the quant metrics against a cutoff
    INPUT:
        :param results_dict: dictionary of parsed alignment file, see output 
                             from parse_usalign_file
        :param ranking_metric_string: user-set string to determine what metric
                                      to apply cutoff on; accepted values = 
                                      'TMSCORE1', 'TMSCORE2', 'MAXTMSCORE',
                                      'MINTMSCORE', 'AVGTMSCORE' (or any other 
                                      combo of capitalization)
        :param metric_cutoff: float, threshold value
    RETURNS:
        :return: return_value, will be 0 if the alignment result passes the 
                 cutoff or 1 otherwise
    """
    if ranking_metric_string.upper() not in ['TMSCORE1','TMSCORE2','MAXTMSCORE','MINTMSCORE','AVGTMSCORE']:
        print(f'Unexpected ranking metric string {ranking_metric_string}. Returning False boolean.', file=sys.stderr, flush=True)
        return_value = 0
    
    elif ranking_metric_string.upper() == 'TMSCORE1' and results_dict['TMscore1'] > metric_cutoff:
        return_value = 1
    
    elif ranking_metric_string.upper() == 'TMSCORE2' and results_dict['TMscore2'] > metric_cutoff:
        return_value = 1
    
    elif ranking_metric_string.upper() == 'MAXTMSCORE' and np.max([results_dict['TMscore1'],results_dict['TMscore2']]) > metric_cutoff:
        return_value = 1
    
    elif ranking_metric_string.upper() == 'MINTMSCORE' and np.min([results_dict['TMscore1'],results_dict['TMscore2']]) > metric_cutoff:
        return_value = 1
    
    elif ranking_metric_string.upper() == 'AVGTMSCORE' and np.mean([results_dict['TMscore1'],results_dict['TMscore2']]) > metric_cutoff:
        return_value = 1
    
    else:
        return_value = 0
    
    return return_value


if __name__ == '__main__':
    import json
    
    #####################################
    # test harness for parse_usalign_file
    #####################################
    # parsing log file for SQ method where the translational and rotational 
    # matrix is not saved
    with open('./examples/3cnaA_to_2pelA_USalign_SQ_noTRMatrix.log','r') as log_file:
        result_dict, start, stop, return_code = parse_usalign_file(log_file,'None')
    with open('./examples/parsed_3cnaA_to_2pelA_USalign_SQ_noTRMatrix.log','w') as file_out:
        json.dump(result_dict,file_out)
    assert result_dict == {'struct1': './3cna_rcsb.pdb', 'struct1_chainID': 'A', 'TMscore1': 0.47057, 'Len1': 237.0, 'd0_1': 5.71, 'struct2': './2pel_rcsb.pdb', 'struct2_chainID': 'A', 'TMscore2': 0.48028, 'Len2': 232.0, 'd0_2': 5.65, 'LenAligned': 118.0, 'RMSD': 1.63, 'SeqIDAli': 0.432}, 'Error in parsing ./examples/parsed_3cnaA_to_2pelA_USalign_SQ_noTRMatrix.log'

    # parsing log file for SQ method where the translational and rotational 
    # matrix is saved
    with open('./examples/3cnaA_to_2pelA_USalign_SQ.log','r') as log_file:
        result_dict, start, stop, return_code = parse_usalign_file(log_file,'None')
    assert np.sum(result_dict['trans_vector'] - np.array([29.61731488,  4.98254133, 36.53334062])) < 1e-5, f"Error in parsing the translational vector in ./examples/3cnaA_to_2pelA_USalign_SQ.log"
    assert np.sum(result_dict['rot_matrix'] - np.array([[ 0.08032842, -0.11836715,  0.98971539],[ 0.99614362, -0.02561858, -0.08391407],[ 0.03528777,  0.99263936,  0.11585278]])) < 1e-5, 'Error in parsing the rotational matrix in ./examples/3cnaA_to_2pelA_USalign_SQ.log'

    # parsing log file for CP method where the translational and rotational 
    # matrix as well as the mapping section are saved
    with open('./examples/3cnaA_to_2pelA_USalign_CP.log','r') as log_file:
        result_dict, start, stop, return_code = parse_usalign_file(log_file,'CP')
    assert result_dict['map_1_to_2']['123'] == ('1', 'THR', 'ALA'), 'Error in parsing alignment mapping section in ./examples/3cnaA_to_2pelA_USalign_CP.log'
    assert result_dict['map_1_to_2']['118'] == ('232', 'ASN', 'THR'), 'Error in parsing alignment mapping section in ./examples/3cnaA_to_2pelA_USalign_CP.log'

    # parsing log file for sNS method where the translational and rotational
    # matrix as well as the mapping section are saved
    with open('./examples/3cnaA_to_2pelA_USalign_sNS.log','r') as log_file:
        result_dict, start, stop, return_code = parse_usalign_file(log_file,'sNS')
    assert result_dict['map_1_to_2']['123'] == ('1', 'THR', 'ALA', 0.42), 'Error in parsing alignment mapping section in ./examples/3cnaA_to_2pelA_USalign_sNS.log'
    assert result_dict['map_1_to_2']['118'] == ('232', 'ASN', 'THR', 1.994), 'Error in parsing alignment mapping section in ./examples/3cnaA_to_2pelA_USalign_sNS.log'

    # parsing log file for fNS method where the translational and rotational
    # matrix as well as the mapping section are saved
    with open('./examples/3cnaA_to_2pelA_USalign_fNS.log','r') as log_file:
        result_dict, start, stop, return_code = parse_usalign_file(log_file,'fNS')
    assert result_dict['map_1_to_2']['123'] == ('1', 'THR', 'ALA', 0.494), 'Error in parsing alignment mapping section in ./examples/3cnaA_to_2pelA_USalign_fNS.log'
    assert result_dict['map_1_to_2']['118'] == ('232', 'ASN', 'THR', 1.99), 'Error in parsing alignment mapping section in ./examples/3cnaA_to_2pelA_USalign_fNS.log'

    #######################################
    # test harness for threshold_comparison
    #######################################
    # using result_dict from the fNS alignment
    test_boolean = threshold_comparison(result_dict,'TMscore1',0.7)
    assert test_boolean, 'error in comparing TMscore1 in ./examples/3cnaA_to_2pelA_USalign_fNS.log against threshold value of 0.7'
    
    test_boolean = threshold_comparison(result_dict,'tMsCoRe1',1.0)
    assert not test_boolean, 'error in comparing TMscore1 in ./examples/3cnaA_to_2pelA_USalign_fNS.log against threshold value of 1.0'

    test_boolean = threshold_comparison(result_dict,'tMsCoRe2',0.7)
    assert test_boolean, 'error in comparing TMscore2 in ./examples/3cnaA_to_2pelA_USalign_fNS.log against threshold value of 0.7'

    test_boolean = threshold_comparison(result_dict,'avgtmscore',0.7)
    assert test_boolean, 'error in comparing AvgTMscore in ./examples/3cnaA_to_2pelA_USalign_fNS.log against threshold value of 0.7'

    test_boolean = threshold_comparison(result_dict,'mintmscore',0.7)
    assert test_boolean, 'error in comparing MinTMscore in ./examples/3cnaA_to_2pelA_USalign_fNS.log against threshold value of 0.7'

    test_boolean = threshold_comparison(result_dict,'maxtmscore',0.7)
    assert test_boolean, 'error in comparing MaxTMscore in ./examples/3cnaA_to_2pelA_USalign_fNS.log against threshold value of 0.7'

    pass

