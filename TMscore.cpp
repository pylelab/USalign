#include "TMscore.h"

using namespace std;

void print_version()
{
    cout << 
"\n"
" *************************************************************************\n"
" *                                 TM-SCORE                              *\n"
" * A scoring function to assess the similarity of protein structures     *\n"
" * Based on statistics:                                                  *\n"
" *       0.0 < TM-score < 0.17, random structural similarity             *\n"
" *       0.5 < TM-score < 1.00, in about the same fold                   *\n"
" * Reference: Yang Zhang and Jeffrey Skolnick, Proteins 2004 57: 702-710 *\n"
" * For comments, please email to: yangzhanglab@umich.edu                 *\n"
" *************************************************************************"
    << endl;
}

void print_extra_help()
{
    cout <<
"Additional options:\n"
"    -a       TM-score normalized by the average length of two structures\n"
"             T or F, (default F)\n"
"\n"
"    -m       Output TM-score rotation matrix\n"
"\n"
"    -d       TM-score scaled by an assigned d0, e.g. 5 Angstroms\n"
"\n"
"    -fast    Fast but slightly inaccurate alignment\n"
"\n"
"    -dir     Perform all-against-all alignment among the list of PDB\n"
"             chains listed by 'chain_list' under 'chain_folder'. Note\n"
"             that the slash is necessary.\n"
"             $ TMscore -dir chain_folder/ chain_list\n"
"\n"
"    -dir1    Use chain2 to search a list of PDB chains listed by 'chain1_list'\n"
"             under 'chain1_folder'. Note that the slash is necessary.\n"
"             $ TMscore -dir1 chain1_folder/ chain1_list chain2\n"
"\n"
"    -dir2    Use chain1 to search a list of PDB chains listed by 'chain2_list'\n"
"             under 'chain2_folder'\n"
"             $ TMscore chain1 -dir2 chain2_folder/ chain2_list\n"
"\n"
"    -suffix  (Only when -dir1 and/or -dir2 are set, default is empty)\n"
"             add file name suffix to files listed by chain1_list or chain2_list\n"
"\n"
"    -atom    4-character atom name used to represent a residue.\n"
"             Default is \" C3'\" for RNA/DNA and \" CA \" for proteins\n"
"             (note the spaces before and after CA).\n"
"\n"
"    -mol     Molecule type: RNA or protein\n"
"             Default is detect molecule type automatically\n"
"\n"
"    -ter     Strings to mark the end of a chain\n"
"             3: (default) TER, ENDMDL, END or different chain ID\n"
"             2: ENDMDL, END, or different chain ID\n"
"             1: ENDMDL or END\n"
"             0: (default in the first C++ TMalign) end of file\n"
"\n"
"    -split   Whether to split PDB file into multiple chains\n"
"             0: (default) treat the whole structure as one single chain\n"
"             1: treat each MODEL as a separate chain (-ter should be 0)\n"
"             2: treat each chain as a seperate chain (-ter should be <=1)\n"
"\n"
"    -outfmt  Output format\n"
"             0: (default) full output\n"
"             1: fasta format compact output\n"
"             2: tabular format very compact output\n"
"            -1: full output, but without version or citation information\n"
"\n"
"    -mirror  Whether to align the mirror image of input structure\n"
"             0: (default) do not align mirrored structure\n"
"             1: align mirror of chain1 to origin chain2\n"
"\n"
"    -het     Whether to align residues marked as 'HETATM' in addition to 'ATOM  '\n"
"             0: (default) only align 'ATOM  ' residues\n"
"             1: align both 'ATOM  ' and 'HETATM' residues\n"
"\n"
"    -infmt1  Input format for chain1\n"
"    -infmt2  Input format for chain2\n"
"            -1: (default) automatically detect PDB or PDBx/mmCIF format\n"
"             0: PDB format\n"
"             1: SPICKER format\n"
"             2: xyz format\n"
"             3: PDBx/mmCIF format\n"
    <<endl;
}

void print_help(bool h_opt=false)
{
    //print_version();
    cout <<
"\n"
" Brief instruction for running TM-score program:\n"
" (For detail: Zhang & Skolnick, Proteins, 2004 57:702-10)\n"
"\n"
" 1. Run TM-score to compare 'model' and 'native':\n"
"     $ TMscore model.pdb native.pdb\n"
"\n"
" 2. Run TM-score to compare two complex structures with multiple chains\n"
"     $ TMscore -c model.pdb native.pdb\n"
"\n"
" 2. TM-score normalized with an assigned scale d0 e.g. 5 A:\n"
"     $ TMscore model.pdb native.pdb -d 5\n"
"\n"
" 3. TM-score normalized by a specific length, e.g. 120 residues:\n"
"     $ TMscore model.pdb native.pdb -l 120\n"
"\n"
" 4. TM-score with superposition output, e.g. 'TM_sup.pdb':\n"
"     $ TMscore model.pdb native.pdb -o TM_sup.pdb\n"
"    To view superimposed atomic model by PyMOL:\n"
"     $ pymol TM_sup.pdb native.pdb\n"
"\n"
" 5. By default, this program assumes that residue pair with the same\n"
"    residue index accross the two structure files are equivalent. This\n"
"    often requires that the residue index in the input structures are\n"
"    renumbered beforehand. Alternatively, residue equivalence can be\n"
"    established by sequence alignment:\n"
"     $ TMscore model.pdb native.pdb -seq\n"
"\n"
    <<endl;

    if (h_opt) print_extra_help();

    exit(EXIT_SUCCESS);
}

int main(int argc, char *argv[])
{
    if (argc < 2) print_help();

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

    bool h_opt = false; // print full help message
    bool v_opt = false; // print version
    bool m_opt = false; // flag for -m, output rotation matrix
    bool o_opt = false; // flag for -o, output superposed structure
    int  a_opt = 0;     // flag for -a, do not normalized by average length
    bool u_opt = false; // flag for -u, normalized by user specified length
    bool d_opt = false; // flag for -d, user specified d0

    double TMcut     =-1;
    int    infmt1_opt=-1;    // PDB or PDBx/mmCIF format for chain_1
    int    infmt2_opt=-1;    // PDB or PDBx/mmCIF format for chain_2
    int    ter_opt   =3;     // TER, END, or different chainID
    int    split_opt =0;     // do not split chain
    int    outfmt_opt=0;     // set -outfmt to full output
    bool   fast_opt  =false; // flags for -fast, fTM-align algorithm
    int    mirror_opt=0;     // do not align mirror
    int    het_opt=0;        // do not read HETATM residues
    string atom_opt  ="auto";// use C alpha atom for protein and C3' for RNA
    string mol_opt   ="auto";// auto-detect the molecule type as protein/RNA
    string suffix_opt="";    // set -suffix to empty
    string dir_opt   ="";    // set -dir to empty
    string dir1_opt  ="";    // set -dir1 to empty
    string dir2_opt  ="";    // set -dir2 to empty
    int    byresi_opt=1;     // TM-score without -c
    vector<string> chain1_list; // only when -dir1 is set
    vector<string> chain2_list; // only when -dir2 is set
    vector<string> chain2parse1;
    vector<string> chain2parse2;
    vector<string> model2parse1;
    vector<string> model2parse2;

    for(int i = 1; i < argc; i++)
    {
        if ( !strcmp(argv[i],"-o") && i < (argc-1) )
        {
            fname_super = argv[i + 1];     o_opt = true; i++;
        }
        else if ( (!strcmp(argv[i],"-u") || !strcmp(argv[i],"-l") ||
                   !strcmp(argv[i],"-L")) && i < (argc-1) )
        {
            Lnorm_ass = atof(argv[i + 1]); u_opt = true; i++;
        }
        else if ( !strcmp(argv[i],"-a") && i < (argc-1) )
        {
            if (!strcmp(argv[i + 1], "T"))      a_opt=true;
            else if (!strcmp(argv[i + 1], "F")) a_opt=false;
            else 
            {
                a_opt=atoi(argv[i + 1]);
                if (a_opt!=-2 && a_opt!=-1 && a_opt!=1)
                    PrintErrorAndQuit("-a must be -2, -1, 1, T or F");
            }
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
        else if (!strcmp(argv[i], "-chain1") )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -chain1");
            split(argv[i+1],chain2parse1,',');
            i++;
        }
        else if (!strcmp(argv[i], "-chain2") )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -chain2");
            split(argv[i+1],chain2parse2,',');
            i++;
        }
        else if (!strcmp(argv[i], "-model1") )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -model1");
            split(argv[i+1],model2parse1,',');
            i++;
        }
        else if (!strcmp(argv[i], "-model2") )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -model2");
            split(argv[i+1],model2parse2,',');
            i++;
        }
        else if (!strcmp(argv[i], "-m") && i < (argc-1) )
        {
            fname_matrix = argv[i + 1];    m_opt = true; i++;
        }// get filename for rotation matrix
        else if (!strcmp(argv[i], "-fast"))
        {
            fast_opt = true;
        }
        else if ( !strcmp(argv[i],"-infmt1") && i < (argc-1) )
        {
            infmt1_opt=atoi(argv[i + 1]); i++;
        }
        else if ( !strcmp(argv[i],"-infmt2") && i < (argc-1) )
        {
            infmt2_opt=atoi(argv[i + 1]); i++;
        }
        else if ( !strcmp(argv[i],"-ter") && i < (argc-1) )
        {
            ter_opt=atoi(argv[i + 1]); i++;
        }
        else if ( !strcmp(argv[i],"-split") && i < (argc-1) )
        {
            split_opt=atoi(argv[i + 1]); i++;
        }
        else if ( !strcmp(argv[i],"-atom") && i < (argc-1) )
        {
            atom_opt=argv[i + 1]; i++;
        }
        else if ( !strcmp(argv[i],"-mol") && i < (argc-1) )
        {
            mol_opt=argv[i + 1]; i++;
        }
        else if ( !strcmp(argv[i],"-dir") && i < (argc-1) )
        {
            dir_opt=argv[i + 1]; i++;
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
        else if ( !strcmp(argv[i],"-c") )
        {
            byresi_opt=2;
        }
        else if ( !strcmp(argv[i],"-seq") )
        {
            byresi_opt=5;
        }
        else if ( !strcmp(argv[i],"-mirror") && i < (argc-1) )
        {
            mirror_opt=atoi(argv[i + 1]); i++;
        }
        else if ( !strcmp(argv[i],"-het") && i < (argc-1) )
        {
            het_opt=atoi(argv[i + 1]); i++;
        }
        else if (xname.size() == 0) xname=argv[i];
        else if (yname.size() == 0) yname=argv[i];
        else PrintErrorAndQuit(string("ERROR! Undefined option ")+argv[i]);
    }

    if(xname.size()==0 || (yname.size()==0 && dir_opt.size()==0) || 
                          (yname.size()    && dir_opt.size()))
    {
        if (h_opt) print_help(h_opt);
        if (v_opt)
        {
            print_version();
            exit(EXIT_FAILURE);
        }
        if (xname.size()==0)
            PrintErrorAndQuit("Please provide input structures");
        else if (yname.size()==0 && dir_opt.size()==0)
            PrintErrorAndQuit("Please provide structure B");
        else if (yname.size() && dir_opt.size())
            PrintErrorAndQuit("Please provide only one file name if -dir is set");
    }

    if (suffix_opt.size() && dir_opt.size()+dir1_opt.size()+dir2_opt.size()==0)
        PrintErrorAndQuit("-suffix is only valid if -dir, -dir1 or -dir2 is set");
    if ((dir_opt.size() || dir1_opt.size() || dir2_opt.size()))
    {
        if (m_opt || o_opt)
            PrintErrorAndQuit("-m or -o cannot be set with -dir, -dir1 or -dir2");
        else if (dir_opt.size() && (dir1_opt.size() || dir2_opt.size()))
            PrintErrorAndQuit("-dir cannot be set with -dir1 or -dir2");
    }

    bool autojustify=(atom_opt=="auto" || atom_opt=="PC4'"); // auto re-pad atom name
    if (atom_opt.size()!=4)
        PrintErrorAndQuit("ERROR! Atom name must have 4 characters, including space.");
    if (mol_opt!="auto" && mol_opt!="protein" && mol_opt!="RNA")
        PrintErrorAndQuit("ERROR! Molecule type must be either RNA or protein.");
    else if (mol_opt=="protein" && atom_opt=="auto")
        atom_opt=" CA ";
    else if (mol_opt=="RNA" && atom_opt=="auto")
        atom_opt=" C3'";

    if (u_opt && Lnorm_ass<=0)
        PrintErrorAndQuit("Wrong value for option -u!  It should be >0");
    if (d_opt && d0_scale<=0)
        PrintErrorAndQuit("Wrong value for option -d!  It should be >0");
    if (outfmt_opt>=2 && (a_opt || u_opt || d_opt))
        PrintErrorAndQuit("-outfmt 2 cannot be used with -a, -u, -L, -d");
    if (byresi_opt>=2 && byresi_opt<=3 && ter_opt>=2)
        PrintErrorAndQuit("-c should be used with -ter <=1");
    if (split_opt==1 && ter_opt!=0)
        PrintErrorAndQuit("-split 1 should be used with -ter 0");
    else if (split_opt==2 && ter_opt!=0 && ter_opt!=1)
        PrintErrorAndQuit("-split 2 should be used with -ter 0 or 1");
    if (split_opt<0 || split_opt>2)
        PrintErrorAndQuit("-split can only be 0, 1 or 2");

    if (m_opt && fname_matrix == "") // Output rotation matrix: matrix.txt
        PrintErrorAndQuit("ERROR! Please provide a file name for option -m!");

    /* parse file list */
    if (dir1_opt.size()+dir_opt.size()==0) chain1_list.push_back(xname);
    else file2chainlist(chain1_list, xname, dir_opt+dir1_opt, suffix_opt);

    if (dir_opt.size())
        for (int i=0;i<chain1_list.size();i++)
            chain2_list.push_back(chain1_list[i]);
    else if (dir2_opt.size()==0) chain2_list.push_back(yname);
    else file2chainlist(chain2_list, yname, dir2_opt, suffix_opt);

    if (byresi_opt>=4)
        cerr<<"WARNING! The residue correspondence between the two structures"
            <<" are automatically established by sequence alignment. Results"
            <<" may be unreliable."<<endl;

    if (outfmt_opt==2)
        cout<<"#PDBchain1\tPDBchain2\tTM1\tTM2\t"
            <<"RMSD\tID1\tID2\tIDali\tL1\tL2\tLali"<<endl;

    /* declare previously global variables */
    vector<vector<string> >PDB_lines1; // text of chain1
    vector<vector<string> >PDB_lines2; // text of chain2
    vector<int> mol_vec1;              // molecule type of chain1, RNA if >0
    vector<int> mol_vec2;              // molecule type of chain2, RNA if >0
    vector<string> chainID_list1;      // list of chainID1
    vector<string> chainID_list2;      // list of chainID2
    int    i,j;                // file index
    int    chain_i,chain_j;    // chain index
    int    r;                  // residue index
    int    xlen, ylen;         // chain length
    int    xchainnum,ychainnum;// number of chains in a PDB file
    char   *seqx, *seqy;       // for the protein sequence 
    double **xa, **ya;         // for input vectors xa[0...xlen-1][0..2] and
                               // ya[0...ylen-1][0..2], in general,
                               // ya is regarded as native structure 
                               // --> superpose xa onto ya
    vector<string> resi_vec1;  // residue index for chain1
    vector<string> resi_vec2;  // residue index for chain2

    /* loop over file names */
    for (i=0;i<chain1_list.size();i++)
    {
        /* parse chain 1 */
        xname=chain1_list[i];
        xchainnum=get_PDB_lines(xname, PDB_lines1, chainID_list1, mol_vec1,
            ter_opt, infmt1_opt, atom_opt, autojustify, split_opt, het_opt,
            chain2parse1, model2parse1);
        if (!xchainnum)
        {
            cerr<<"Warning! Cannot parse file: "<<xname
                <<". Chain number 0."<<endl;
            continue;
        }
        for (chain_i=0;chain_i<xchainnum;chain_i++)
        {
            xlen=PDB_lines1[chain_i].size();
            if (mol_opt=="RNA") mol_vec1[chain_i]=1;
            else if (mol_opt=="protein") mol_vec1[chain_i]=-1;
            if (!xlen)
            {
                cerr<<"Warning! Cannot parse file: "<<xname
                    <<". Chain length 0."<<endl;
                continue;
            }
            else if (xlen<3)
            {
                cerr<<"Sequence is too short <3!: "<<xname<<endl;
                continue;
            }
            NewArray(&xa, xlen, 3);
            seqx = new char[xlen + 1];
            xlen = read_PDB(PDB_lines1[chain_i], xa, seqx, 
                resi_vec1, byresi_opt);
            if (mirror_opt) for (r=0;r<xlen;r++) xa[r][2]=-xa[r][2];

            for (j=(dir_opt.size()>0)*(i+1);j<chain2_list.size();j++)
            {
                /* parse chain 2 */
                if (PDB_lines2.size()==0)
                {
                    yname=chain2_list[j];
                    ychainnum=get_PDB_lines(yname, PDB_lines2, chainID_list2,
                        mol_vec2, ter_opt, infmt2_opt, atom_opt, autojustify,
                        split_opt, het_opt, chain2parse2, model2parse2);
                    if (!ychainnum)
                    {
                        cerr<<"Warning! Cannot parse file: "<<yname
                            <<". Chain number 0."<<endl;
                        continue;
                    }
                }
                for (chain_j=0;chain_j<ychainnum;chain_j++)
                {
                    ylen=PDB_lines2[chain_j].size();
                    if (mol_opt=="RNA") mol_vec2[chain_j]=1;
                    else if (mol_opt=="protein") mol_vec2[chain_j]=-1;
                    if (!ylen)
                    {
                        cerr<<"Warning! Cannot parse file: "<<yname
                            <<". Chain length 0."<<endl;
                        continue;
                    }
                    else if (ylen<3)
                    {
                        cerr<<"Sequence is too short <3!: "<<yname<<endl;
                        continue;
                    }
                    NewArray(&ya, ylen, 3);
                    seqy = new char[ylen + 1];
                    ylen = read_PDB(PDB_lines2[chain_j], ya, seqy,
                        resi_vec2, byresi_opt);

                    if (byresi_opt) extract_aln_from_resi(sequence,
                        seqx,seqy,resi_vec1,resi_vec2,byresi_opt);

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

                    double rmsd_d0_out=0;
                    int L_lt_d=0;
                    double GDT_list[5]={0,0,0,0,0}; // 0.5, 1, 2, 4, 8
                    double maxsub=0;
                    TM1=TM2=TM3=TM4=TM5=0;

                    /* entry function for structure alignment */
                    TMscore_main(
                        xa, ya, seqx, seqy,
                        t0, u0, TM1, TM2, TM3, TM4, TM5,
                        d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                        seqM, seqxA, seqyA,
                        rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                        xlen, ylen, sequence, Lnorm_ass, d0_scale,
                        a_opt, u_opt, d_opt, fast_opt,
                        mol_vec1[chain_i]+mol_vec2[chain_j],
                        GDT_list,maxsub,TMcut);

                    /* print result */
                    if (outfmt_opt==0) print_version();
                    output_TMscore_results(
                        xname.substr(dir1_opt.size()+dir_opt.size()),
                        yname.substr(dir2_opt.size()+dir_opt.size()),
                        chainID_list1[chain_i],
                        chainID_list2[chain_j],
                        xlen, ylen, t0, u0, TM1, TM2, 
                        TM3, TM4, TM5, rmsd0, d0_out,
                        seqM.c_str(), seqxA.c_str(), seqyA.c_str(), Liden,
                        n_ali8, L_ali, TM_ali, rmsd_ali,
                        TM_0, d0_0, d0A, d0B,
                        Lnorm_ass, d0_scale, d0a, d0u, 
                        (m_opt?fname_matrix:"").c_str(),
                        outfmt_opt, ter_opt, 
                        (o_opt?fname_super:"").c_str(),
                        a_opt, u_opt, d_opt, mirror_opt,
                        L_lt_d, rmsd_d0_out, GDT_list, maxsub,
                        split_opt, resi_vec1, resi_vec2);

                    /* Done! Free memory */
                    seqM.clear();
                    seqxA.clear();
                    seqyA.clear();
                    DeleteArray(&ya, ylen);
                    delete [] seqy;
                    resi_vec2.clear();
                } // chain_j
                if (chain2_list.size()>1)
                {
                    yname.clear();
                    for (chain_j=0;chain_j<ychainnum;chain_j++)
                        PDB_lines2[chain_j].clear();
                    PDB_lines2.clear();
                    chainID_list2.clear();
                    mol_vec2.clear();
                }
            } // j
            PDB_lines1[chain_i].clear();
            DeleteArray(&xa, xlen);
            delete [] seqx;
            resi_vec1.clear();
        } // chain_i
        xname.clear();
        PDB_lines1.clear();
        chainID_list1.clear();
        mol_vec1.clear();
    } // i
    if (chain2_list.size()==1)
    {
        yname.clear();
        for (chain_j=0;chain_j<ychainnum;chain_j++)
            PDB_lines2[chain_j].clear();
        PDB_lines2.clear();
        resi_vec2.clear();
        chainID_list2.clear();
        mol_vec2.clear();
    }
    chain1_list.clear();
    chain2_list.clear();
    sequence.clear();
    vector<string>().swap(chain2parse1);
    vector<string>().swap(chain2parse2);
    vector<string>().swap(model2parse1);
    vector<string>().swap(model2parse2);
    return 0;
}
