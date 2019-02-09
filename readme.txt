==============================================================================
   TM-align: protein and RNA structure alignment by TM-score superposition.

   This program was written by (in reverse chronological order)
   Chengxin Zhang, Sha Gong, Jianjie Wu, and Jianyi Yang
   at Yang Zhang lab, Department of Computational Medicine and Bioinformatics,
   University of Michigan, 100 Washtenaw Ave, Ann Arbor, MI 48109-2218.
   Please report issues to yangzhanglab@umich.edu

   References to cite:
   S Gong, C Zhang, Y Zhang. Bioinformatics (2019)
   Y Zhang, J Skolnick. Nucl Acids Res 33, 2302-9 (2005)

   DISCLAIMER:
     Permission to use, copy, modify, and distribute this program for 
     any purpose, with or without fee, is hereby granted, provided that
     the notices on the head, the reference information, and this
     copyright notice appear in all copies or substantial portions of 
     the Software. It is provided "as is" without express or implied 
     warranty.

   *************** updating history ********************************
   2012/01/24: A C/C++ code of TM-align was constructed by J Yang
   2016/05/21: Several updates of this program were made by J Wu, including
              (1) fixed several compiling bugs
              (2) made I/O of C/C++ version consistent with the Fortran version
              (3) added outputs including full-atom and ligand structures
              (4) added options of '-i', '-I' and '-m'
   2016/05/25: fixed a bug on PDB file reading
   2018/06/04: Several updates were made by C Zhang, including
              (1) Fixed bug in reading PDB files with negative residue index,
                  at the expense of the '-o' option now only being able to
                  output superposed structure instead of full rasmol script.
              (2) Implemented the fTM-align algorithm (by the '-fast' option)
                  as described in R Dong, S Pan, Z Peng, Y Zhang, J Yang
                  (2018) Nucleic acids research. gky430.
              (3) Included option to perform TM-align against a whole 
                  folder of PDB files. A full list of options not available
                  in the Fortran version can be explored by TMalign -h
   2018/07/27: Added the -byresi option for TM-score superposition without
               re-alignment as in TMscore and TMscore -c
   2018/08/07: Added the -dir option
   2018/08/14: Added the -split option
   2018/08/16: Added the -infmt1, -infmt2 options.
               TMalign can now read .gz and .bz2 compressed files.
   2018/10/20: C Zhang and S Gong updated the RNA alignment part of
               the program. Changes include:
              (1) new d0 calculation for RNA.
              (2) secondary structure assignment for RNA.
              (3) automatic detection of molecule type (protein vs RNA).
   2019/01/07: C Zhang added support for PDBx/mmCIF format.
   2019/02/09: Fixed asymmetric alignment bug.
===============================================================================

=========================
 How to install TM-align
=========================
To compile the program in your Linux computer, simply enter

 make

or

 g++ -static -O3 -ffast-math -lm -o TMalign TMalign.cpp

The '-static' flag should be removed on Mac OS, which does not support
building static executables.

=====================
 How to use TM-align
=====================
You can run the program without arguments to obtain a brief instruction

 ./TMalign structure1.pdb structure2.pdb

===================
 Fortran version
===================
You can download the fortran version of TM-align from
https://zhanglab.ccmb.med.umich.edu/TM-align/

This C++ version of TM-align implemented several features not available in the
fortran version, including RNA alignment and batch alignment of multiple 
structures. A full list of available options can be explored by:
  ./TMalign -h

02/09/2019
