/*
===============================================================================
   This is a re-implementation of TM-align algorithm in C/C++. The code was
   is written by Jianyi Yang and later updated by Jianjie Wu at The Yang Zhang
   lab, Department of Computational Medicine and Bioinformatics, University of
   Michigan, 100 Washtenaw Avenue, Ann Arbor, MI 48109-2218. Please report bugs
   and questions to zhng@umich.edu

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
   2018/06/04: Fixed a bug in PDB file with negative residue number. Added
               options -fast, -dir1, -dir2, -suffix, -atom, -ter, -outfmt.
               Re-write the file reading function to reduce the number of
               times a PDB file need to be read.
===============================================================================
*/
#include "TMalign.h"

using namespace std;

void print_version()
{
    cout << 
"\n"
" *****************************************************************************\n"
" * TM-align (Version 20180604): A protein structural alignment algorithm     *\n"
" * Reference: Y Zhang and J Skolnick, Nucl Acids Res 33, 2302-9 (2005)       *\n"
" * Please email your comments and suggestions to Yang Zhang (zhng@umich.edu) *\n"
" *****************************************************************************"
    << endl;
}

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
"            -1: full output, but without version and citation information\n"
    <<endl;
}

void print_help(bool h_opt=false)
{
    print_version();
    cout <<
"\n"
"Usage: TMalign PDB1.pdb PDB2.pdb [Options]\n"
"\n"
"Options:\n"
"    -u    TM-score normalized by user assigned length (the same as -L)\n"
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
"    (Options -u, -a, -d, -o won't change the final structure alignment)\n\n"
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
    string xname       = "";
    string yname       = "";
    string fname_super = ""; // file name for superposed structure
    string fname_lign  = ""; // file name for user alignment
    string fname_matrix= ""; // file name for output matrix
    vector<string> sequence; // get value from alignment file
    double Lnorm_ass, d0_scale;

    bool A_opt = false; // marker for whether structure A is specified
    bool B_opt = false; // marker for whether structure B is specified
    bool h_opt = false; // print full help message
    bool v_opt = false; // print version
    bool m_opt = false; // flag for -m, output rotation matrix
    bool i_opt = false; // flag for -i, with user given initial alignment
    bool I_opt = false; // flag for -I, stick to user given alignment
    bool o_opt = false; // flag for -o, output superposed structure
    bool a_opt = false; // flag for -a, normalized by average length
    bool u_opt = false; // flag for -u, normalized by user specified length
    bool d_opt = false; // flag for -d, user specified d0

    int ter_opt = 3;    // TER, END, or different chainID
    int outfmt_opt=0;   // set -outfmt to full output
    bool fast_opt = false;  // flags for -fast, fTM-align algorithm
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
            fname_super = argv[i + 1];     o_opt = true; i++;
        }
        else if ( (!strcmp(argv[i],"-u") || 
                   !strcmp(argv[i],"-L")) && i < (argc-1) )
        {
            Lnorm_ass = atof(argv[i + 1]); u_opt = true; i++;
        }
        else if ( !strcmp(argv[i],"-a") && i < (argc-1) )
        {
            if (!strcmp(argv[i + 1], "T"))      a_opt=true;
            else if (!strcmp(argv[i + 1], "F")) a_opt=false;
            else PrintErrorAndQuit("Wrong value for option -a! It should be T or F");
            i++;
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
            fname_lign = argv[i + 1];      i_opt = true; i++;
        }
        else if (!strcmp(argv[i], "-m") && i < (argc-1) )
        {
            fname_matrix = argv[i + 1];    m_opt = true; i++;
        }// get filename for rotation matrix
        else if (!strcmp(argv[i], "-I") && i < (argc-1) )
        {
            fname_lign = argv[i + 1];      I_opt = true; i++;
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
                xname=argv[i]; B_opt = true;
            }
            else if (nameIdx == 1)
            {
                yname=argv[i]; A_opt = true;
            }
            nameIdx++;
        }
    }

    if(!B_opt || !A_opt)
    {

        if (h_opt) print_help(h_opt);
        if (v_opt)
        {
            print_version();
            exit(EXIT_FAILURE);
        }
    }

    if (!A_opt) PrintErrorAndQuit("Please provide structure A");
    if (!B_opt) PrintErrorAndQuit("Please provide structure B");

    if (suffix_opt.size() && dir1_opt.size()==0 && dir2_opt.size()==0)
        PrintErrorAndQuit("-suffix is only valid if -dir1 or -dir2 is set");
    if ((dir1_opt.size() || dir2_opt.size()) && (m_opt || o_opt))
        PrintErrorAndQuit("-m or -o cannot be set with -dir1 or -dir2");
    if (atom_opt.size()!=4)
        PrintErrorAndQuit("ERROR! atom name must have 4 characters, including space.");

    if (i_opt && I_opt)
        PrintErrorAndQuit("ERROR! -I and -i cannot be used together");
    if (u_opt && Lnorm_ass<=0)
        PrintErrorAndQuit("Wrong value for option -u!  It should be >0");
    if (d_opt && d0_scale<=0)
        PrintErrorAndQuit("Wrong value for option -d!  It should be >0");
    if (outfmt_opt>=2 && (a_opt || u_opt || d_opt))
        PrintErrorAndQuit("-outfmt 2 cannot be used with -a, -u, -L, -d");

    /* read initial alignment file from 'align.txt' */
    string basename = string(argv[0]);
    int idx = basename.find_last_of("\\");
    basename = basename.substr(0, idx + 1);
    if (i_opt || I_opt)// Ask TM-align to start with an alignment,specified in fasta file 'align.txt'
    {
        if (fname_lign == "")
            PrintErrorAndQuit("Please provide a file name for option -i!");
        // open alignment file
        int n_p = 0;// number of structures in alignment file
        string line;

        string fullpath = basename + fname_lign;
        ifstream fileIn(fullpath.c_str());
        if (fileIn.is_open())
        {
            while (fileIn.good())
            {
                getline(fileIn, line);
                if (line.compare(0, 1, ">") == 0)// Flag for a new structure
                {
                    if (n_p >= 2) break;
                    sequence.push_back("");
                    n_p++;
                }
                else if (n_p > 0 && line!="") sequence.back()+=line;
            }
            fileIn.close();
        }
        else
            PrintErrorAndQuit("ERROR! Alignment file does not exist.");

        if (n_p < 2)
            PrintErrorAndQuit("ERROR: Fasta format is wrong, two proteins should be included.");
        if (sequence[0].size() != sequence[1].size())
            PrintErrorAndQuit("ERROR! FASTA file is wrong. The length in alignment should be equal respectively to the two aligned proteins.");
        if (I_opt)
        {
            int aligned_resNum=0;
            for (int i=0;i<sequence[0].size();i++) 
                aligned_resNum+=(sequence[0][i]!='-' && sequence[1][i]!='-');
            if (aligned_resNum<3)
                PrintErrorAndQuit("ERROR! Superposition is undefined for <3 aligned residues.");
        }
    }

    if (m_opt && fname_matrix == "") // Output rotation matrix: matrix.txt
        PrintErrorAndQuit("ERROR! Please provide a file name for option -m!");

    /* parse file list */
    if (dir1_opt.size()==0)
        chain1_list.push_back(xname);
    else
    {
        ifstream fp(xname.c_str());
        if (! fp.is_open())
        {
            char message[5000];
            sprintf(message, "Can not open file: %s\n", xname.c_str());
            PrintErrorAndQuit(message);
        }
        string line;
        while (fp.good())
        {
            getline(fp, line);
            if (! line.size()) continue;
            chain1_list.push_back(dir1_opt+Trim(line)+suffix_opt);
        }
        fp.close();
        line.clear();
    }

    if (dir2_opt.size()==0)
        chain2_list.push_back(yname);
    else
    {
        ifstream fp(yname.c_str());
        if (! fp.is_open())
        {
            char message[5000];
            sprintf(message, "Can not open file: %s\n", yname.c_str());
            PrintErrorAndQuit(message);
        }
        string line;
        while (fp.good())
        {
            getline(fp, line);
            if (! line.size()) continue;
            chain2_list.push_back(dir2_opt+Trim(line)+suffix_opt);
        }
        fp.close();
        line.clear();
    }

    if (outfmt_opt==2)
        cout<<"#PDBchain1\tPDBchain2\tTM1\tTM2\t"
            <<"RMSD\tID1\tID2\tIDali\tL1\tL2\tLali"<<endl;

    /* declare previously global variables */
    vector<string> PDB_lines1; // text of chain1
    vector<string> PDB_lines2; // text of chain2
    int    xlen, ylen;         // chain length
    char   *seqx, *seqy;       // for the protein sequence 
    int    *secx, *secy;       // for the secondary structure 
    int    *xresno, *yresno;   // residue number for fragment gapless threading
    double **xa, **ya;         // for input vectors xa[0...xlen-1][0..2] and
                               // ya[0...ylen-1][0..2], in general,
                               // ya is regarded as native structure 
                               // --> superpose xa onto ya

    /* loop over file names */
    for (int i=0;i<chain1_list.size();i++)
    {
        /* parse chain 1 */
        xname=chain1_list[i];
        xlen=get_PDB_lines(xname.c_str(),PDB_lines1,ter_opt,atom_opt);
        if (!xlen)
        {
            cerr<<"Warning! Cannot parse file: "<<xname
                <<". Chain length 0."<<endl;
            PDB_lines1.clear();
            continue;
        }
        NewArray(&xa, xlen, 3);
        seqx = new char[xlen + 1];
        secx = new int[xlen];
        xresno = new int[xlen];
        xlen = read_PDB(PDB_lines1, xa, seqx, xresno);
        make_sec(xa, xlen, secx); // secondary structure assignment

        for (int j=0;j<chain2_list.size();j++)
        {
            /* parse chain 2 */
            if (PDB_lines2.size()==0)
            {
                yname=chain2_list[j];
                ylen=get_PDB_lines(yname.c_str(),PDB_lines2,ter_opt,atom_opt);
                if (!ylen)
                {
                    cerr<<"Warning! Cannot parse file: "<<yname
                        <<". Chain length 0."<<endl;
                    PDB_lines2.clear();
                    continue;
                }
                NewArray(&ya, ylen, 3);
                seqy = new char[ylen + 1];
                yresno = new int[ylen];
                secy = new int[ylen];
                ylen = read_PDB(PDB_lines2, ya, seqy, yresno);
                make_sec(ya, ylen, secy);
            }

            /* declare variable specific to this pair of TMalign */
            double t0[3], u0[3][3];
            double TM1, TM2;
            double TM3, TM4, TM5;     // for a_opt, u_opt, d_opt
            double d0_0, TM_0;
            double d0A, d0B, d0u, d0a;
            double d0_out=5.0;
            string seqM, seqxA, seqyA;// for output alignment
            double rmsd0 = 0.0;
            int L_ali;                // Aligned length in standard_TMscore
            double Liden=0;
            double TM_ali, rmsd_ali;  // TMscore and rmsd in standard_TMscore
            int n_ali=0;
            int n_ali8=0;

            /* entry function for structure alignment */
            TMalign_main(
                xa, ya, xresno, yresno, seqx, seqy, secx, secy,
                t0, u0, TM1, TM2, TM3, TM4, TM5,
                d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                seqM, seqxA, seqyA,
                rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                xlen, ylen, sequence, Lnorm_ass, d0_scale,
                i_opt, I_opt, a_opt, u_opt, d_opt, fast_opt);

            /* print result */
            if (outfmt_opt==0) print_version();
            output_results(
                xname.substr(dir1_opt.size()).c_str(),
                yname.substr(dir2_opt.size()).c_str(),
                xlen, ylen, t0, u0, TM1, TM2, 
                TM3, TM4, TM5, rmsd0, d0_out,
                seqM.c_str(), seqxA.c_str(), seqyA.c_str(), Liden,
                n_ali8, n_ali, L_ali, TM_ali, rmsd_ali,
                TM_0, d0_0, d0A, d0B,
                Lnorm_ass, d0_scale, d0a, d0u, fname_matrix.c_str(),
                outfmt_opt, ter_opt, fname_super.c_str(),
                i_opt, I_opt, o_opt, a_opt, u_opt, d_opt);

            /* Done! Free memory */
            seqM.clear();
            seqxA.clear();
            seqyA.clear();
            if (chain2_list.size()>1)
            {
                yname.clear();
                PDB_lines2.clear();
                DeleteArray(&ya, ylen);
                delete [] seqy;
                delete [] secy;
                delete [] yresno;
            }
        }
        xname.clear();
        PDB_lines1.clear();
        DeleteArray(&xa, xlen);
        delete [] seqx;
        delete [] secx;
        delete [] xresno;
    }
    if (chain2_list.size()==1)
    {
        yname.clear();
        PDB_lines2.clear();
        DeleteArray(&ya, ylen);
        delete [] seqy;
        delete [] secy;
        delete [] yresno;
    }
    chain1_list.clear();
    chain2_list.clear();
    sequence.clear();

    t2 = clock();
    float diff = ((float)t2 - (float)t1)/CLOCKS_PER_SEC;
    printf("Total CPU time is %5.2f seconds\n", diff);
    return 0;
}
