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
__version__ = '20220623'
__license__ = 'BSD-2-Clause'

from pymol import cmd, CmdException
import subprocess
import tempfile
import os

def usalign(mobile, target, args='', exe='USalign', transform=1):
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

    exe = cmd.exp_path(exe)
    if args=='""':
        args=''
    if len(args)>2 and args[0]=='"' and args[-1]=='"':
        args=args[1:-1]
    if not "-outfmt" in args:
        args+=" -outfmt -1"
    args = ' '.join([exe, mobile_filename, target_filename, args, '-m -'])
    print(args)

    try:
        process = subprocess.Popen(args, stdout=subprocess.PIPE,
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
        if line.strip().startswith('----'):
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

# autocompletion
cmd.auto_arg[0].update({
    'usalign': cmd.auto_arg[0]['align'],
})
cmd.auto_arg[1].update({
    'usalign': cmd.auto_arg[1]['align'],
})
