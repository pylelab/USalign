#!/usr/bin/env pymol
'''
PyMOL plugin for US-align

USAGE: 

    usalign mobile, fix [,args [,exe]]

INSTALLATION

    Install this script as a PyMOL plugin by 
    "Plugin" - "Plugin Manager" - "Install New Plugin"

    This plugin depends on the binary executable of US-align, which must be
    available within a directory specified by PATH. You can get the PATH
    value within PyMOL by the following command:
    
    print(os.getenv('PATH'))
'''
#This script is partly based on tmalign plugin by Thomas Holder available at
#https://github.com/Pymol-Scripts/Pymol-script-repo/blob/master/tmalign.py

from __future__ import print_function

__author__ = 'Chengxin Zhang'
__version__ = '20220924'
__license__ = 'BSD-2-Clause'

from pymol import cmd, CmdException
import subprocess
import tempfile
import os
import platform

def get_usalign_path(exe="USalign"):
    if platform.system().lower().startswith("win"):
        exe+=".exe"
    filename = os.path.join(os.path.dirname(os.path.abspath(__file__)),exe)
    if os.path.isfile(filename):
        return filename
    else:
        for p in os.getenv("PATH").split(os.pathsep):
            filename=os.path.join(p,exe)
            if os.path.isfile(filename):
                return filename
    print("ERROR! Cannot locate %s at %s or at %s"%(exe,
        os.path.dirname(os.path.abspath(__file__)),os.getenv("PATH")))
    print("Please put the USalign executable at one of the aforementioned paths")
    return exe

def usalign(mobile, target, args='', exe='', transform=1):
    '''
USAGE

    usalign mobile, target [, args [, exe ]]

ARGUMENTS

    mobile, target = string: atom selections

    args = string: Extra arguments such as -mm and -byresi

    exe = string: Path to USalign executable {default: USalign}

CITATION

    Zhang C, Shine M, Pyle AM, Zhang Y. bioRxiv 2022.04.18.488565.
    https://github.com/pylelab/USalign
    '''

    mobile_filename = tempfile.mktemp('.pdb', 'mobile')
    target_filename = tempfile.mktemp('.pdb', 'target')
    mobile_ca_sele = '(%s) and (not hetatm) and alt +A' % (mobile)
    target_ca_sele = '(%s) and (not hetatm) and alt +A' % (target)
    if not "-atom" in args:
        mobile_ca_sele+=" and (name CA or name C3')"
        target_ca_sele+=" and (name CA or name C3')"

    cmd.save(mobile_filename, mobile_ca_sele)
    cmd.save(target_filename, target_ca_sele)

    if len(exe)==0:
        exe=get_usalign_path("USalign")
    if args=='""':
        args=''
    if len(args)>2 and args[0]=='"' and args[-1]=='"':
        args=args[1:-1]
    if not "-outfmt" in args:
        args+=" -outfmt -1"
    args = ' '.join([exe, mobile_filename, target_filename, args, '-m -'])
    print(args)

    try:
        process = subprocess.Popen(args, stdout=subprocess.PIPE, shell=True,
                universal_newlines=True)
        lines = process.stdout.readlines()
    except OSError:
        print('Cannot execute "%s", please provide full path to USalign executable' % (args))
        raise CmdException
    finally:
        os.remove(mobile_filename)
        os.remove(target_filename)

    rowcount = 0
    matrix = []
    for line in iter(lines):
        print(line.rstrip())
        if line.strip().startswith('------ The rotation matrix to rotate '):
            rowcount = 1
        elif 4 >= rowcount and rowcount> 0:
            if rowcount >= 2:
                a = list(map(float, line.split()))
                matrix.extend(a[2:5])
                matrix.append(a[1])
            rowcount += 1

    assert len(matrix) == 3 * 4
    matrix.extend([0, 0, 0, 1])

    if int(transform):
        cmd.transform_selection('byobject (%s)' % (mobile), matrix, homogenous=1)
    return

# pymol commands
cmd.extend('usalign', usalign)
cmd.extend('USalign', usalign)

# autocompletion
cmd.auto_arg[0].update({ 'usalign': cmd.auto_arg[0]['align'], })
cmd.auto_arg[1].update({ 'usalign': cmd.auto_arg[1]['align'], })
cmd.auto_arg[0].update({ 'USalign': cmd.auto_arg[0]['align'], })
cmd.auto_arg[1].update({ 'USalign': cmd.auto_arg[1]['align'], })
