===============================================================================
   This is a re-implementation of TM-align algorithm in C/C++. The code was 
   is written by Jianyi Yang and later updated by Jianjie Wu at The Yang Zhang 
   lab, Department of Computational Medicine and Bioinformatics, University of 
   Michigan, 100 Washtenaw Avenue, Ann Arbor, MI 48109-2218. Please report bugs 
   and questions to jianjiew@umich.edu or zhng@umich.edu

   DISCLAIMER:
     Permission to use, copy, modify, and distribute this program for 
     any purpose, with or without fee, is hereby granted, provided that
     the notices on the head, the reference information, and this
     copyright notice appear in all copies or substantial portions of 
     the Software. It is provided "as is" without express or implied 
     warranty.

   *************** updating history ********************************
   2012/01/24: A C/C++ code of TM-align was constructed by Jianyi Yang
   2016/05/21: Several updates of this program were made by Jianjie Wu, including
              (1) fixed several compiling bugs
	      (2) made I/O of C/C++ version consistent with the Fortran version
	      (3) added outputs including full-atom and ligand structures
	      (4) added options of '-i', '-I' and '-m'
   2016/05/25: fixed a bug on PDB file reading
===============================================================================

=========================
 How to install TM-align
=========================
to compile the program in your Linux computer, simply enter

 make

or

 g++ -static -O3 -ffast-math -lm -o TMalign TMalign.cpp

=====================
 How to use TM-align
=====================
you can run the program without arguments to obtain a brief instruction

 ./TMalign structure1.pdb structure2.pdb

====================
 About this program
====================
   This program is written by Jianyi Yang at
   Yang Zhang lab
   And it is updated by Jianjie Wu at
   Yang Zhang lab
   Department of Computational Medicine and Bioinformatics 
   University of Michigan 
   100 Washtenaw Avenue, Ann Arbor, MI 48109-2218 
           
   Please report bugs and questions to jianjiew@umich.edu or zhng@umich.edu

===================
 Fortran version
===================
You can download the fortran version of TM-align from
http://zhanglab.ccmb.med.umich.edu/TM-align/

05/21/2016
