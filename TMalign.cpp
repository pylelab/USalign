/*
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
*/
#include "basic_define.h"
#include "global_var.h"
#include "param_set.h"

using namespace std;


#include "basic_fun.h"
#include "NW.h"
#include "Kabsch.h"
#include "TMalign.h"

void print_extra_help()
{
    cout <<
"Additional options: \n"
"    -fast    Fast but slightly inaccurate alignment\n"
"\n"
"    -dir1    Use chain2 to search a list of PDB chains listed by 'chain1_list'\n"
"             under 'chain1_folder'. Note that the slash is necessary.\n"
"             $ TMalign -dir1 chain1_folder/ chain1_list chain2\n"
"\n"
"    -dir2    Use chain1 to search a list of PDB chains listed by 'chain2_list'\n"
"             under 'chain2_folder'\n"
"             $ TMalign chain1 -dir2 chain2_folder/ chain2_list\n"
"\n"
"    -suffix  (Only when -dir1 and/or -dir2 are set, default is empty)\n"
"             add file name suffix to files listed by chain1_list or chain2_list\n"
"\n"
"    -atom    4-character atom name used to represent a residue\n"
"             default is \" CA \" (note the space before and after CA)\n"
"\n"
"    -ter     Strings to mark the end of a chain\n"
"             3: (current default) 'TER', 'END', or different chain ID\n"
"             2: 'END', or different chain ID\n"
"             1: 'END'\n"
"             0: (default in the first C++ version of TMalign) end of file\n"
"\n"
"    -outfmt  Output format\n"
"             0: (default) full output\n"
"             1: fasta format compact output\n"
"             2: tabular format very compact output\n"
    <<endl;
}

void print_help(bool h_opt=false)
{
    cout <<
"\n"
"*****************************************************************************\n"
"* TM-align (Version "<< TMalign_version
                      <<    "): An algorithm for protein structure alignment *\n"
"* Based on statistics:                                                      *\n"
"*          0.0 < TM-score < 0.30, random structural similarity              *\n"
"*          0.5 < TM-score < 1.00, in about the same fold                    *\n"
"* Reference: Y Zhang and J Skolnick, Nucl Acids Res 33, 2302-9 (2005)       *\n"
"* Please email your comments and suggestions to Yang Zhang (zhng@umich.edu) *\n"
"*****************************************************************************\n"
"\n"
"Usage: TMalign PDB1.pdb PDB2.pdb [Options]\n"
"\n"
"Options:\n"
"    -u    TM-score normalized by user assigned length\n"
"          warning: it should be >= minimum length of the two structures\n"
"          otherwise, TM-score may be >1\n"
"\n"
"    -a    TM-score normalized by the average length of two structures\n"
"          T or F, (default F)\n"
"\n"
"    -i    Ask TM-align to start with an alignment, specified in fasta\n"
"          file 'align.txt'\n"
"\n"
"    -I    Ask TM-align to stick the alignment to 'align.txt'\n"
"\n"
"    -m    Output TM-align rotation matrix\n"
"\n"
"    -d    TM-score scaled by an assigned d0, e.g. 5 Angstroms\n"
"\n"
"    -o    output the superposition of chain1 to TM.sup\n"
"          $ TMalign chain1 chain2 -o TM.sup\n"
"          To view superimposed full-atom structures:\n"
"          $ pymol TM.sup chain2\n"
"\n"
"    -v    print the version of TM-align\n"
"\n"
"    -h    print help message\n"
"\n"
"    (Options -u, -a, -d -o won't change the final structure alignment)\n\n"
"Example usages:\n"
"    TMalign PDB1.pdb PDB2.pdb\n"
"    TMalign PDB1.pdb PDB2.pdb -u 100 -d 5.0\n"
"    TMalign PDB1.pdb PDB2.pdb -a T -o PDB1.sup\n"
"    TMalign PDB1.pdb PDB2.pdb -i align.txt\n"
"    TMalign PDB1.pdb PDB2.pdb -m matrix.txt\n"
    <<endl;

    if (h_opt) print_extra_help();

    exit(EXIT_SUCCESS);
}

int main(int argc, char *argv[])
{
    if (argc < 2) print_help();


    clock_t t1, t2;
    t1 = clock();

    /**********************/
    /*    get argument    */
    /**********************/
    char xname[MAXLEN], yname[MAXLEN],  Lnorm_ave[MAXLEN];
    bool A_opt, B_opt, h_opt=false;
    bool v_opt = false;
    int ter_opt = 3; // TER, END, or different chainID
    int outfmt_opt=0;  // set -outfmt to full output
    A_opt = B_opt = o_opt = a_opt = u_opt = d_opt = false;
    i_opt = false;// set -i flag to be false
    m_opt = false;// set -m flag to be false
    char fname_lign[MAXLEN] = "";
    char fname_matrix[MAXLEN] = "";// set names to ""
    I_opt = false;// set -I flag to be false
    fast_opt = false;// set -fast flag to be false
    string atom_opt=" CA "; // use C alpha atom to represent a residue
    string suffix_opt=""; // set -suffix to empty
    string dir1_opt="";   // set -dir1 to empty
    string dir2_opt="";   // set -dir2 to empty
    vector<string> chain1_list; // only when -dir1 is set
    vector<string> chain2_list; // only when -dir2 is set

    int nameIdx = 0;
    for(int i = 1; i < argc; i++)
    {
        if ( !strcmp(argv[i],"-o") && i < (argc-1) )
        {
            strcpy(out_reg, argv[i + 1]);      o_opt = true; i++;
        }
        else if ( !strcmp(argv[i],"-u") && i < (argc-1) )
        {
            Lnorm_ass = atof(argv[i + 1]); u_opt = true; i++;
        }
        else if ( !strcmp(argv[i],"-a") && i < (argc-1) )
        {
            strcpy(Lnorm_ave, argv[i + 1]);     a_opt = true; i++;
        }
        else if ( !strcmp(argv[i],"-d") && i < (argc-1) )
        {
            d0_scale = atof(argv[i + 1]); d_opt = true; i++;
        }
        else if ( !strcmp(argv[i],"-v") )
        {
            v_opt = true;
        }
        else if ( !strcmp(argv[i],"-h") )
        {
            h_opt = true;
        }
        else if ( !strcmp(argv[i],"-i") && i < (argc-1) )
        {
            strcpy(fname_lign, argv[i + 1]);      i_opt = true; i++;
        }
        else if (!strcmp(argv[i], "-m") && i < (argc-1) )
        {
            strcpy(fname_matrix, argv[i + 1]);      m_opt = true; i++;
        }// get filename for rotation matrix
        else if (!strcmp(argv[i], "-I") && i < (argc-1) )
        {
            strcpy(fname_lign, argv[i + 1]);      I_opt = true; i++;
        }
        else if (!strcmp(argv[i], "-fast"))
        {
            fast_opt = true;
        }
        else if ( !strcmp(argv[i],"-ter") && i < (argc-1) )
        {
            ter_opt=atoi(argv[i + 1]); i++;
        }
        else if ( !strcmp(argv[i],"-atom") && i < (argc-1) )
        {
            atom_opt=argv[i + 1]; i++;
            if (atom_opt.size()!=4)
                PrintErrorAndQuit("ERROR! atom name must be 4 characters, including space.");
        }
        else if ( !strcmp(argv[i],"-dir1") && i < (argc-1) )
        {
            dir1_opt=argv[i + 1]; i++;
        }
        else if ( !strcmp(argv[i],"-dir2") && i < (argc-1) )
        {
            dir2_opt=argv[i + 1]; i++;
        }
        else if ( !strcmp(argv[i],"-suffix") && i < (argc-1) )
        {
            suffix_opt=argv[i + 1]; i++;
        }
        else if ( !strcmp(argv[i],"-outfmt") && i < (argc-1) )
        {
            outfmt_opt=atoi(argv[i + 1]); i++;
        }
        else
        {
            if (nameIdx == 0)
            {
                strcpy(xname, argv[i]);      B_opt = true;
            }
            else if (nameIdx == 1)
            {
                strcpy(yname, argv[i]);      A_opt = true;
            }
            nameIdx++;
        }
    }

    if(!B_opt || !A_opt)
    {

        if( h_opt ) print_help(h_opt);

        if(v_opt)
        {
            print_version();
            exit(EXIT_FAILURE);
        }
    }

    if( !A_opt )
        PrintErrorAndQuit("Please provide structure A");
    if( !B_opt )
        PrintErrorAndQuit("Please provide structure B");

    if (suffix_opt.size() && dir1_opt.size()==0 && dir2_opt.size()==0)
        PrintErrorAndQuit("-suffix is only valid if -dir1 or -dir2 is set");
    if ((dir1_opt.size() || dir2_opt.size()) && (m_opt || o_opt))
        PrintErrorAndQuit("-m or -o cannot be set with -dir1 or -dir2");

    if( a_opt )
    {
        if(!strcmp(Lnorm_ave, "T"))
        {
        }
        else if(!strcmp(Lnorm_ave, "F"))
        {
            a_opt=false;
        }
        else
        {
            cout << "Wrong value for option -a!  It should be T or F" << endl;
            exit(EXIT_FAILURE);
        }
    }
    if( u_opt )
    {
        if(Lnorm_ass<=0)
        {
            cout << "Wrong value for option -u!  It should be >0" << endl;
            exit(EXIT_FAILURE);
        }
    }
    if( d_opt )
    {
        if(d0_scale<=0)
        {
            cout << "Wrong value for option -d!  It should be >0" << endl;
            exit(EXIT_FAILURE);
        }
    }
    ////// read initial alignment file from 'alignment.txt' //////
    string basename = string(argv[0]);
    int idx = basename.find_last_of("\\");
    basename = basename.substr(0, idx + 1);
    if (i_opt || I_opt)// Ask TM-align to start with an alignment,specified in fasta file 'align.txt'
    {
        if (fname_lign == "")
        {
            cout << "Please provide a file name for option -i!" << endl;
            exit(EXIT_FAILURE);
        }
        // open alignment file
        int n_p = 0;// number of structures in alignment file
        string line;

        string fullpath = basename + fname_lign;
        char path1[1000];
        strcpy(path1, fullpath.c_str());
        ifstream fileIn(path1);
        if (fileIn.is_open())
        {
            bool bContinue = true;
            while (fileIn.good() && bContinue)
            {
                getline(fileIn, line);
                if (line.compare(0, 1, ">") == 0)// Flag for a new structure
                {
                    strcpy(sequence[n_p], "");
                    n_p++;
                    if (n_p > 2)
                        bContinue = false;
                }
                else// Read data
                {
                    if (n_p > 0 && line!="")
                    {
                        strcat(sequence[n_p-1], line.c_str());
                    }
                }
            }
            fileIn.close();
        }
        else
        {
            cout << "\nAlignment file does not exist.\n";
            exit(EXIT_FAILURE);
        }

        if (n_p < 2)
        {
            cout << "\nERROR: Fasta format is wrong, two proteins should be included.\n";
            exit(EXIT_FAILURE);
        }
        else
        {
            if (strlen(sequence[0]) != strlen(sequence[1]))
            {
                cout << "\nWarning: FASTA format may be wrong, the length in alignment should be equal respectively to the aligned proteins.\n";
                exit(EXIT_FAILURE);
            }
        }
    }

    if (m_opt)// Output TM - align rotation matrix: matrix.txt
    {
        if (fname_matrix == "")
        {
            cout << "Please provide a file name for option -m!" << endl;
            exit(EXIT_FAILURE);
        }
    }

    /* parse file list */
    if (dir1_opt.size()==0)
        chain1_list.push_back(xname);
    else
    {
        ifstream fp(xname);
        if (! fp.is_open())
        {
		    char message[5000];
		    sprintf(message, "Can not open file: %s\n", xname);
            PrintErrorAndQuit(message);
        }
        string line;
        string filename;
        while (fp.good())
        {
            getline(fp, line);
            if (! line.size()) continue;
            filename=dir1_opt+Trim(line)+suffix_opt;
            if (isfile_openable(filename))
                chain1_list.push_back(filename);
            else
                cerr<<"Warning! Skipped inaccesible file"<<filename<<endl;
        }
        fp.close();
        line.clear();
        filename.clear();
    }

    if (dir2_opt.size()==0)
        chain2_list.push_back(yname);
    else
    {
        ifstream fp(yname);
        if (! fp.is_open())
        {
		    char message[5000];
		    sprintf(message, "Can not open file: %s\n", yname);
            PrintErrorAndQuit(message);
        }
        string line;
        string filename;
        while (fp.good())
        {
            getline(fp, line);
            if (! line.size()) continue;
            filename=dir2_opt+Trim(line)+suffix_opt;
            if (isfile_openable(filename))
                chain2_list.push_back(filename);
            else
                cerr<<"Warning! Skipped inaccesible file"<<filename<<endl;
        }
        fp.close();
        line.clear();
        filename.clear();
    }

    /* loop over file names */
    if (outfmt_opt==2)
        cout<<"#PDBchain1\tPDBchain2\tTM1\tTM2\tRMSD\tID1\tID2\tIDali\tL1\tL2\tLali"<<endl;

    for (int i=0;i<chain1_list.size();i++)
    {
        strcpy(xname,chain1_list[i].c_str());
        for (int j=0;j<chain2_list.size();j++)
        {
            strcpy(yname,chain2_list[j].c_str());

            /* load data */
            load_PDB_allocate_memory(xname, yname, ter_opt, atom_opt);

            /* entry function for structure alignment */
            TMalign_main(xname, yname, fname_matrix, ter_opt, 
                dir1_opt, dir2_opt, outfmt_opt);

            /* Done! Free memory */
            free_memory();
        }
    }

    chain1_list.clear();
    chain2_list.clear();

    t2 = clock();
    float diff = ((float)t2 - (float)t1)/CLOCKS_PER_SEC;
    printf("Total running time is %5.2f seconds\n", diff);
    return 0;
}
