#!/usr/bin/env pymol
'''
PyMOL plugin for US-align

USAGE: 

    usalign mobile, fix [,args [,exe]]
    usalign_msta pdb1, pdb2, pdb3 [...] [,exe=...]

INSTALLATION

    Install this script as a PyMOL plugin by 
    "Plugin" - "Plugin Manager" - "Install New Plugin"

    This plugin depends on the binary executable of US-align, which must be
    available within a directory specified by PATH. You can get the PATH
    value within PyMOL by the following command:
    
    print(os.getenv('PATH'))

    On Mac OS, the US-align binary executable can be installed by
    conda install -c bioconda usalign
'''
#This script is partly based on tmalign plugin by Thomas Holder available at
#https://github.com/Pymol-Scripts/Pymol-script-repo/blob/master/tmalign.py

from __future__ import print_function

__author__ = 'Chengxin Zhang'
__version__ = '20241201'
__license__ = 'BSD-2-Clause'

from pymol import cmd, CmdException
import subprocess
import tempfile
import os
import platform

def check_executable(exe="USalign"):
    try:
        process = subprocess.Popen('"'+exe+'"', stdout=subprocess.PIPE, shell=True,
                universal_newlines=True)
        if "version" in process.stdout.read():
            return True
    except OSError:
        print(exe+" not executable")
        return False
    print(exe+" not executable")
    return False

def get_usalign_path(exe="USalign"):
    if platform.system().lower().startswith("win"):
        exe+=".exe"
    filename = os.path.join(os.path.dirname(os.path.abspath(__file__)),exe)
    if os.path.isfile(filename) and check_executable(filename):
        return filename
    else:
        for p in os.getenv("PATH").split(os.pathsep):
            filename=os.path.join(p,exe)
            if os.path.isfile(filename) and check_executable(filename):
                return filename
    print("ERROR! Cannot locate %s at %s or at %s"%(exe,
        os.path.dirname(os.path.abspath(__file__)),os.getenv("PATH")))
    print("Please put the USalign executable at one of the aforementioned paths")
    if platform.system().lower().startswith("darwin"):
        print("You may install US-align executable by:")
        print("conda install -c bioconda usalign")
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

    Zhang C, Shine M, Pyle AM, Zhang Y (2022) Nat Methods, 19, 1109-1115.
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
    args = ' '.join(['"'+exe+'"', mobile_filename, target_filename, args, '-m -'])
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

def usalign_msta(*args, **kwargs):
    '''
USAGE

    usalign_msta pdb1, pdb2, pdb3 [...] [,args] [,exe]

ARGUMENTS

    pdb1, pdb2, pdb3 ...  = string: atom selections

    args = string: Extra arguments such as -mm and -byresi

    exe = string: Path to USalign executable {default: USalign}

CITATION

    Zhang C, Shine M, Pyle AM, Zhang Y (2022) Nat Methods, 19, 1109-1115.
    https://github.com/pylelab/USalign
    '''

    target_list = []
    argument = ""
    for t,target in enumerate(args):
        target_list.append(target)
        print("[%d] %s"%(t,target))

    exe=''
    transform=1
    argument=''
    for key in kwargs:
        if key=="exe":
            exe=kwargs[key]
        elif key=="transform":
            transform=kwargs[key]
        elif key=="args":
            argument=kwargs[key]
        elif "self" in key:
            continue
        else:
            print("ignore keyword argument %s= %s"%(key,kwargs[key]))
    
    if len(target_list)<=1:
        print("ERROR! less than three targets two align")
        return
    elif len(target_list)==2:
        print("ERROR! only two targets; please use the 'usalign' command instead")
        return
    if "-mm " in argument and not "-mm 4" in argument:
        print("ERROR! -mm cannot be used with usalign_msta")
        return
    if "-m " in argument:
        print("ERROR! do not set -m for usalign_msta")
        return
    if "-outfmt " in argument and not "-outfmt -1" in argument and \
        not "-outfmt 0" in argument and not "-outfmt 1" in argument:
        print("ERROR! use -outfmt 0 with usalign_msta")
        return
    
    if len(exe)==0:
        exe=get_usalign_path("USalign")
        
    args = argument
    if args=='""':
        args=''
    elif len(args)>2 and args[0]=='"' and args[-1]=='"':
        args=args[1:-1]
    if not "-outfmt " in args:
        args+=" -outfmt 0"
    
    tmpdir   = tempfile.mkdtemp(prefix="tmp_msta")
    listfile = os.path.join(tmpdir,"list")
    matrixfile = os.path.join(tmpdir,"matrix.txt")
    fp = open(listfile,'w')
    filename_list = []
    for t,target in enumerate(target_list):
        target_filename = os.path.join(tmpdir, "%d.pdb"%t)
        target_ca_sele = '(%s) and (not hetatm) and alt +A' % (target)
        if not "-atom " in args:
            target_ca_sele+=" and (name CA or name C3')"
        cmd.save(target_filename, target_ca_sele)
        filename_list.append(target_filename)
        fp.write("%d.pdb\n"%t)
    fp.close()

    args = ' '.join(['"'+exe+'"', "-mm 4 -dir",tmpdir,listfile, args, "-m",matrixfile])
    print(args)
    try:
        process = subprocess.Popen(args, stdout=subprocess.PIPE, shell=True,
                universal_newlines=True)
        lines = process.stdout.readlines()
        for line in iter(lines):
            print(line.rstrip())
    except OSError:
        print('Cannot execute "%s", please provide full path to USalign executable' % (args))
        raise CmdException
    
    for target_filename in filename_list:
        os.remove(target_filename)
    os.remove(listfile)

    if not int(transform):
        return
    
    fp=open(matrixfile)
    matrixtxt=fp.read()
    fp.close()
    for block in matrixtxt.split('------ The rotation matrix to rotate '):
        matrix = []
        t = -1
        for rowcount,line in enumerate(block.splitlines()):
            if rowcount==0 and " to " in line:
                target = line.strip().split('.pdb')[0].lstrip('/').lstrip('\\')
                t = int(target)
            elif rowcount>=5:
                break
            elif rowcount>=2:
                a = list(map(float, line.split()))
                matrix.extend(a[2:5])
                matrix.append(a[1])

        if t==-1:
            continue
        if len(matrix) != 3*4:
            print("len(matrix)=%d != 3*4 for %s"%(len(matrix),target_list[t]))
            continue
        matrix.extend([0, 0, 0, 1])
        print('------ The rotation matrix to rotate '+target_list[t]+' '+block)
        cmd.transform_selection('byobject (%s)' % (target_list[t]), matrix, homogenous=1)
    return

# pymol commands
cmd.extend('usalign', usalign)
cmd.extend('USalign', usalign)
cmd.extend('usalign_msta', usalign_msta)
cmd.extend('USalign_MSTA', usalign_msta)

# autocompletion
cmd.auto_arg[0].update({ 'usalign': cmd.auto_arg[0]['align'], })
cmd.auto_arg[1].update({ 'usalign': cmd.auto_arg[1]['align'], })
cmd.auto_arg[0].update({ 'USalign': cmd.auto_arg[0]['align'], })
cmd.auto_arg[1].update({ 'USalign': cmd.auto_arg[1]['align'], })

cmd.auto_arg[0].update({ 'usalign_msta': cmd.auto_arg[0]['align'], })
cmd.auto_arg[1].update({ 'usalign_msta': cmd.auto_arg[1]['align'], })
cmd.auto_arg[0].update({ 'USalign_MSTA': cmd.auto_arg[0]['align'], })
cmd.auto_arg[1].update({ 'USalign_MSTA': cmd.auto_arg[1]['align'], })
