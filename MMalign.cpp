/* command line argument parsing and document of MMalign main program */

#include "MMalign.h"

using namespace std;

void print_version()
{
    cout << 
"\n"
" **********************************************************************\n"
" * MM-align (Version 20231222): complex structure alignment           *\n"
" * References: S Mukherjee, Y Zhang. Nucl Acids Res 37(11):e83 (2009) *\n"
" * Please email comments and suggestions to yangzhanglab@umich.edu    *\n"
" **********************************************************************"
    << endl;
}

void print_extra_help()
{
    cout <<
"Additional options:\n"
"    -fast    Fast but slightly inaccurate alignment\n"
"\n"
"    -dir1    Use a list of PDB chains listed by 'chain1_list' under\n"
"             'chain1_folder' as all chains for the first complex.\n"
"             Note that the slash is necessary.\n"
"             $ MMalign -dir1 chain1_folder/ chain1_list complex2\n"
"\n"
"    -dir2    Use a list of PDB chains listed by'chain2_list'\n"
"             under 'chain2_folder' as all chains for the second complex.\n"
"             $ MMalign complex1 -dir2 chain2_folder/ chain2_list\n"
"\n"
"    -suffix  (Only when -dir1 and/or -dir2 are set, default is empty)\n"
"             add file name suffix to files listed by chain1_list or chain2_list\n"
"\n"
"    -atom    4-character atom name used to represent a residue.\n"
"             Default is \" C3'\" for RNA/DNA and \" CA \" for proteins\n"
"             (note the spaces before and after CA).\n"
"\n"
"    -mol     Types of molecules to align\n""Molecule type: RNA or protein\n"
"             auto   : (default) align both proteins and nucleic acids\n"
"             protein: only align proteins\n"
"             RNA    : only align nucleic acids (RNA and DNA)\n"
"\n"
"    -split   Whether to split PDB file into multiple chains\n"
"             2: (default) treat each chain as a seperate chain (-ter should be <=1)\n"
"             1: treat each MODEL as a separate chain (-ter should be 0)\n"
"                and joins all chains in a MODEL into a single chain.\n"
"\n"
"    -outfmt  Output format\n"
"             0: (default) full output\n"
"             1: fasta format compact output\n"
"             2: tabular format very compact output\n"
"            -1: full output, but without version or citation information\n"
"\n"
"    -TMcut   -1: (default) do not consider TMcut\n"
"             Values in [0.5,1): Do not proceed with TM-align for this\n"
"                 structure pair if TM-score is unlikely to reach TMcut.\n"
"                 TMcut is normalized is set by -a option:\n"
"                 -2: normalized by longer structure length\n"
"                 -1: normalized by shorter structure length\n"
"                  0: (default, same as F) normalized by second structure\n"
"                  1: same as T, normalized by average structure length\n"
"\n"
"    -mirror  Whether to align the mirror image of input structure\n"
"             0: (default) do not align mirrored structure\n"
"             1: align mirror of chain1 to origin chain2\n"
"\n"
"    -het     Whether to align residues marked as 'HETATM' in addition to 'ATOM  '\n"
"             0: (default) only align 'ATOM  ' residues\n"
"             1: align both 'ATOM  ' and 'HETATM' residues\n"
"\n"
"    -infmt1  Input format for complex1\n"
"    -infmt2  Input format for complex2\n"
"            -1: (default) automatically detect PDB or PDBx/mmCIF format\n"
"             0: PDB format\n"
"             1: SPICKER format\n"
"             2: xyz format\n"
"             3: PDBx/mmCIF format\n"
    <<endl;
}

void print_help(bool h_opt=false)
{
    print_version();
    cout <<
"\n"
"Usage: MMalign complex1.pdb complex2.pdb [Options]\n"
"\n"
"Options:\n"
"    -a    TM-score normalized by the average length of two structures\n"
"          T or F, (default F)\n"
"\n"
"    -m    Output MM-align rotation matrix\n"
"\n"
"    -d    TM-score scaled by an assigned d0, e.g. 5 Angstroms\n"
"\n"
"    -o    Output the superposition of complex1.pdb to MM_sup.pdb\n"
"          $ MMalign complex1.pdb complex2.pdb -o MM_sup.pdb\n"
"          To view superposed full-atom structures:\n"
"          $ pymol MM_sup.pdb complex2.pdb\n"
"\n"
"    -full Whether to show full alignment result, including alignment of\n"
"          individual chains. T or F, (default F)\n"
"\n"
"    -ter  Whether to read all MODELs in a multi-model structure file\n"
"          1: (default) only read the first model, recommended for alignment\n"
"             of asymetric units.\n"
"          0: read all MODEL, recomended for alignment of biological\n"
"             assemblies, i.e., biological units (biounits).\n"
"\n"
"    -v    Print the version of MM-align\n"
"\n"
"    -h    Print the full help message\n"
"\n"
"    (Options -a, -d, -m, -o won't change the final structure alignment)\n\n"
"Example usages:\n"
"    MMalign complex1.pdb complex2.pdb\n"
"    MMalign complex1.pdb complex2.pdb -d 5.0\n"
"    MMalign complex1.pdb complex2.pdb -a T -o complex1.sup\n"
"    MMalign complex1.pdb complex2.pdb -m matrix.txt\n"
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
    double d0_scale    =0;

    bool h_opt = false; // print full help message
    bool v_opt = false; // print version
    bool m_opt = false; // flag for -m, output rotation matrix
    bool o_opt = false; // flag for -o, output superposed structure
    int  a_opt = 0;     // flag for -a, do not normalized by average length
    bool d_opt = false; // flag for -d, user specified d0

    bool   full_opt  = false;// do not show chain level alignment
    double TMcut     =-1;
    int    infmt1_opt=-1;    // PDB or PDBx/mmCIF format for chain_1
    int    infmt2_opt=-1;    // PDB or PDBx/mmCIF format for chain_2
    int    ter_opt   =1;     // ENDMDL or END
    int    split_opt =2;     // split by chain
    int    outfmt_opt=0;     // set -outfmt to full output
    bool   fast_opt  =false; // flags for -fast, fTM-align algorithm
    int    mirror_opt=0;     // do not align mirror
    int    het_opt   =0;     // do not read HETATM residues
    string atom_opt  ="auto";// use C alpha atom for protein and C3' for RNA
    string mol_opt   ="auto";// auto-detect the molecule type as protein/RNA
    string suffix_opt="";    // set -suffix to empty
    string dir1_opt  ="";    // set -dir1 to empty
    string dir2_opt  ="";    // set -dir2 to empty
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
        else if ( !strcmp(argv[i],"-full") && i < (argc-1) )
        {
            if (!strcmp(argv[i + 1], "T"))      full_opt=true;
            else if (!strcmp(argv[i + 1], "F")) full_opt=false;
            else PrintErrorAndQuit("-full must be T or F");
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
        else if ( !strcmp(argv[i],"-TMcut") && i < (argc-1) )
        {
            TMcut=atof(argv[i + 1]); i++;
        }
        else if ( !strcmp(argv[i],"-het") && i < (argc-1) )
        {
            het_opt=atoi(argv[i + 1]); i++;
        }
        else if (xname.size() == 0) xname=argv[i];
        else if (yname.size() == 0) yname=argv[i];
        else PrintErrorAndQuit(string("ERROR! Undefined option ")+argv[i]);
    }

    if(yname.size()==0)
    {
        if (h_opt) print_help(h_opt);
        if (v_opt)
        {
            print_version();
            exit(EXIT_FAILURE);
        }
        if (xname.size()==0)
            PrintErrorAndQuit("Please provide input structures");
        PrintErrorAndQuit("Please provide the second input structure");
    }

    if (suffix_opt.size() && dir1_opt.size()+dir2_opt.size()==0)
        PrintErrorAndQuit("-suffix is only valid if -dir1 or -dir2 is set");
    if ((dir1_opt.size() || dir2_opt.size()) && (m_opt || o_opt))
        PrintErrorAndQuit("-m or -o cannot be set with -dir1 or -dir2");
    if (atom_opt.size()!=4)
        PrintErrorAndQuit("ERROR! Atom name must have 4 characters, including space.");
    if (mol_opt!="auto" && mol_opt!="protein" && mol_opt!="RNA")
        PrintErrorAndQuit("ERROR! Molecule type must be either RNA or protein.");
    else if (mol_opt=="protein" && atom_opt=="auto")
        atom_opt=" CA ";
    else if (mol_opt=="RNA" && atom_opt=="auto")
        atom_opt=" C3'";

    if (d_opt && d0_scale<=0)
        PrintErrorAndQuit("Wrong value for option -d!  It should be >0");
    if (outfmt_opt>=2 && (a_opt || d_opt))
        PrintErrorAndQuit("-outfmt 2 cannot be used with -a, -d");
    if (ter_opt!=0 && ter_opt!=1)
        PrintErrorAndQuit("-ter should be 1 or 0");
    if (split_opt!=1 && split_opt!=2)
        PrintErrorAndQuit("-split should be 1 or 2");
    else if (split_opt==1 && ter_opt!=0)
        PrintErrorAndQuit("-split 1 should be used with -ter 0");

    if (m_opt && fname_matrix == "") // Output rotation matrix: matrix.txt
        PrintErrorAndQuit("ERROR! Please provide a file name for option -m!");

    /* parse file list */
    if (dir1_opt.size()==0) chain1_list.push_back(xname);
    else file2chainlist(chain1_list, xname, dir1_opt, suffix_opt);

    if (dir2_opt.size()==0) chain2_list.push_back(yname);
    else file2chainlist(chain2_list, yname, dir2_opt, suffix_opt);

    if (outfmt_opt==2)
        cout<<"#PDBchain1\tPDBchain2\tTM1\tTM2\t"
            <<"RMSD\tID1\tID2\tIDali\tL1\tL2\tLali"<<endl;

    /* declare previously global variables */
    vector<vector<vector<double> > > xa_vec; // structure of complex1
    vector<vector<vector<double> > > ya_vec; // structure of complex2
    vector<vector<char> >seqx_vec; // sequence of complex1
    vector<vector<char> >seqy_vec; // sequence of complex2
    vector<vector<char> >secx_vec; // secondary structure of complex1
    vector<vector<char> >secy_vec; // secondary structure of complex2
    vector<int> mol_vec1;          // molecule type of complex1, RNA if >0
    vector<int> mol_vec2;          // molecule type of complex2, RNA if >0
    vector<string> chainID_list1;  // list of chainID1
    vector<string> chainID_list2;  // list of chainID2
    vector<int> xlen_vec;          // length of complex1
    vector<int> ylen_vec;          // length of complex2
    int    i,j;                    // chain index
    int    xlen, ylen;             // chain length
    double **xa, **ya;             // structure of single chain
    char   *seqx, *seqy;           // for the protein sequence 
    char   *secx, *secy;           // for the secondary structure 
    int    xlen_aa,ylen_aa;        // total length of protein
    int    xlen_na,ylen_na;        // total length of RNA/DNA
    vector<string> resi_vec1;  // residue index for chain1
    vector<string> resi_vec2;  // residue index for chain2

    /* parse complex */
    parse_chain_list(chain1_list, xa_vec, seqx_vec, secx_vec, mol_vec1,
        xlen_vec, chainID_list1, ter_opt, split_opt, mol_opt, infmt1_opt,
        atom_opt, false, mirror_opt, het_opt, xlen_aa, xlen_na, o_opt,
        resi_vec1, chain2parse1, model2parse1);
    if (xa_vec.size()==0) PrintErrorAndQuit("ERROR! 0 chain in complex 1");
    parse_chain_list(chain2_list, ya_vec, seqy_vec, secy_vec, mol_vec2,
        ylen_vec, chainID_list2, ter_opt, split_opt, mol_opt, infmt2_opt,
        atom_opt, false, 0, het_opt, ylen_aa, ylen_na, o_opt,
        resi_vec2, chain2parse2, model2parse2);
    if (ya_vec.size()==0) PrintErrorAndQuit("ERROR! 0 chain in complex 2");
    int len_aa=getmin(xlen_aa,ylen_aa);
    int len_na=getmin(xlen_na,ylen_na);
    if (a_opt)
    {
        len_aa=(xlen_aa+ylen_aa)/2;
        len_na=(xlen_na+ylen_na)/2;
    }

    map<int,int> chainmap;

    /* perform monomer alignment if there is only one chain */
    if (xa_vec.size()==1 && ya_vec.size()==1)
    {
        xlen = xlen_vec[0];
        ylen = ylen_vec[0];
        seqx = new char[xlen+1];
        seqy = new char[ylen+1];
        secx = new char[xlen+1];
        secy = new char[ylen+1];
        NewArray(&xa, xlen, 3);
        NewArray(&ya, ylen, 3);
        copy_chain_data(xa_vec[0],seqx_vec[0],secx_vec[0], xlen,xa,seqx,secx);
        copy_chain_data(ya_vec[0],seqy_vec[0],secy_vec[0], ylen,ya,seqy,secy);
        
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
        vector<double> do_vec;

        /* entry function for structure alignment */
        TMalign_main(xa, ya, seqx, seqy, secx, secy,
            t0, u0, TM1, TM2, TM3, TM4, TM5,
            d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
            seqM, seqxA, seqyA, do_vec,
            rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
            xlen, ylen, sequence, 0, d0_scale,
            0, a_opt, false, d_opt, fast_opt,
            mol_vec1[0]+mol_vec2[0],TMcut);

        /* print result */
        output_results(
            xname.substr(dir1_opt.size()),
            yname.substr(dir2_opt.size()),
            chainID_list1[0], chainID_list2[0],
            xlen, ylen, t0, u0, TM1, TM2, TM3, TM4, TM5, rmsd0, d0_out,
            seqM.c_str(), seqxA.c_str(), seqyA.c_str(), Liden,
            n_ali8, L_ali, TM_ali, rmsd_ali, TM_0, d0_0, d0A, d0B,
            0, d0_scale, d0a, d0u, (m_opt?fname_matrix:"").c_str(),
            outfmt_opt, ter_opt, true, split_opt, o_opt, fname_super,
            0, a_opt, false, d_opt, mirror_opt, resi_vec1, resi_vec2);

        /* clean up */
        seqM.clear();
        seqxA.clear();
        seqyA.clear();
        delete[]seqx;
        delete[]seqy;
        delete[]secx;
        delete[]secy;
        DeleteArray(&xa,xlen);
        DeleteArray(&ya,ylen);
        chain1_list.clear();
        chain2_list.clear();
        sequence.clear();
        do_vec.clear();

        vector<vector<vector<double> > >().swap(xa_vec); // structure of complex1
        vector<vector<vector<double> > >().swap(ya_vec); // structure of complex2
        vector<vector<char> >().swap(seqx_vec); // sequence of complex1
        vector<vector<char> >().swap(seqy_vec); // sequence of complex2
        vector<vector<char> >().swap(secx_vec); // secondary structure of complex1
        vector<vector<char> >().swap(secy_vec); // secondary structure of complex2
        mol_vec1.clear();       // molecule type of complex1, RNA if >0
        mol_vec2.clear();       // molecule type of complex2, RNA if >0
        chainID_list1.clear();  // list of chainID1
        chainID_list2.clear();  // list of chainID2
        xlen_vec.clear();       // length of complex1
        ylen_vec.clear();       // length of complex2

        t2 = clock();
        float diff = ((float)t2 - (float)t1)/CLOCKS_PER_SEC;
        printf("#Total CPU time is %5.2f seconds\n", diff);
        return 0;
    }

    /* declare TM-score tables */
    int chain1_num=xa_vec.size();
    int chain2_num=ya_vec.size();
    vector<string> tmp_str_vec(chain2_num,"");
    double **TMave_mat;
    double **ut_mat; // rotation matrices for all-against-all alignment
    int ui,uj,ut_idx;
    NewArray(&TMave_mat,chain1_num,chain2_num);
    NewArray(&ut_mat,chain1_num*chain2_num,4*3);
    vector<vector<string> >seqxA_mat(chain1_num,tmp_str_vec);
    vector<vector<string> > seqM_mat(chain1_num,tmp_str_vec);
    vector<vector<string> >seqyA_mat(chain1_num,tmp_str_vec);

    double maxTMmono=-1;
    int maxTMmono_i,maxTMmono_j;

    /* get all-against-all alignment */
    if (len_aa+len_na>500) fast_opt=true;
    for (i=0;i<chain1_num;i++)
    {
        xlen=xlen_vec[i];
        if (xlen<3)
        {
            for (j=0;j<chain2_num;j++) TMave_mat[i][j]=-1;
            continue;
        }
        seqx = new char[xlen+1];
        secx = new char[xlen+1];
        NewArray(&xa, xlen, 3);
        copy_chain_data(xa_vec[i],seqx_vec[i],secx_vec[i],
            xlen,xa,seqx,secx);

        for (j=0;j<chain2_num;j++)
        {
            ut_idx=i*chain2_num+j;
            for (ui=0;ui<4;ui++)
                for (uj=0;uj<3;uj++) ut_mat[ut_idx][ui*3+uj]=0;
            ut_mat[ut_idx][0]=1;
            ut_mat[ut_idx][4]=1;
            ut_mat[ut_idx][8]=1;

            if (mol_vec1[i]*mol_vec2[j]<0) //no protein-RNA alignment
            {
                TMave_mat[i][j]=-1;
                continue;
            }

            ylen=ylen_vec[j];
            if (ylen<3)
            {
                TMave_mat[i][j]=-1;
                continue;
            }
            seqy = new char[ylen+1];
            secy = new char[ylen+1];
            NewArray(&ya, ylen, 3);
            copy_chain_data(ya_vec[j],seqy_vec[j],secy_vec[j],
                ylen,ya,seqy,secy);

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
            vector<double> do_vec;
            int Lnorm_tmp=len_aa;
            if (mol_vec1[i]+mol_vec2[j]>0) Lnorm_tmp=len_na;

            /* entry function for structure alignment */
            TMalign_main(xa, ya, seqx, seqy, secx, secy,
                t0, u0, TM1, TM2, TM3, TM4, TM5,
                d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                seqM, seqxA, seqyA, do_vec,
                rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                xlen, ylen, sequence, Lnorm_tmp, d0_scale,
                0, false, true, false, fast_opt,
                mol_vec1[i]+mol_vec2[j],TMcut);

            /* store result */
            for (ui=0;ui<3;ui++)
                for (uj=0;uj<3;uj++) ut_mat[ut_idx][ui*3+uj]=u0[ui][uj];
            for (uj=0;uj<3;uj++) ut_mat[ut_idx][9+uj]=t0[uj];
            seqxA_mat[i][j]=seqxA;
            seqyA_mat[i][j]=seqyA;
            TMave_mat[i][j]=TM4*Lnorm_tmp;
            if (TMave_mat[i][j]>maxTMmono)
            {
                maxTMmono=TMave_mat[i][j];
                maxTMmono_i=i;
                maxTMmono_j=j;
            }

            /* clean up */
            seqM.clear();
            seqxA.clear();
            seqyA.clear();

            delete[]seqy;
            delete[]secy;
            DeleteArray(&ya,ylen);
            do_vec.clear();
        }

        delete[]seqx;
        delete[]secx;
        DeleteArray(&xa,xlen);
    }

    /* calculate initial chain-chain assignment */
    int *assign1_list; // value is index of assigned chain2
    int *assign2_list; // value is index of assigned chain1
    assign1_list=new int[chain1_num];
    assign2_list=new int[chain2_num];
    double total_score=enhanced_greedy_search(TMave_mat, assign1_list,
        assign2_list, chain1_num, chain2_num);
    if (total_score<=0) PrintErrorAndQuit("ERROR! No assignable chain");

    /* refine alignment for large oligomers */
    int aln_chain_num=count_assign_pair(assign1_list,chain1_num);
    bool is_oligomer=(aln_chain_num>=3);
    if (aln_chain_num==2) // dimer alignment
    {
        int na_chain_num1,na_chain_num2,aa_chain_num1,aa_chain_num2;
        count_na_aa_chain_num(na_chain_num1,aa_chain_num1,mol_vec1);
        count_na_aa_chain_num(na_chain_num2,aa_chain_num2,mol_vec2);

        /* align protein-RNA hybrid dimer to another hybrid dimer */
        if (na_chain_num1==1 && na_chain_num2==1 && 
            aa_chain_num1==1 && aa_chain_num2==1) is_oligomer=false;
        /* align pure protein dimer or pure RNA dimer */
        else if ((getmin(na_chain_num1,na_chain_num2)==0 && 
                    aa_chain_num1==2 && aa_chain_num2==2) ||
                 (getmin(aa_chain_num1,aa_chain_num2)==0 && 
                    na_chain_num1==2 && na_chain_num2==2))
        {
            adjust_dimer_assignment(xa_vec,ya_vec,xlen_vec,ylen_vec,mol_vec1,
                mol_vec2,assign1_list,assign2_list,seqxA_mat,seqyA_mat);
            is_oligomer=false; // cannot refiner further
        }
        else is_oligomer=true; /* align oligomers to dimer */
    }

    if (aln_chain_num>=3 || is_oligomer) // oligomer alignment
    {
        /* extract centroid coordinates */
        double **xcentroids;
        double **ycentroids;
        NewArray(&xcentroids, chain1_num, 3);
        NewArray(&ycentroids, chain2_num, 3);
        double d0MM=getmin(
            calculate_centroids(xa_vec, chain1_num, xcentroids),
            calculate_centroids(ya_vec, chain2_num, ycentroids));

        /* refine enhanced greedy search with centroid superposition */
        //double het_deg=check_heterooligomer(TMave_mat, chain1_num, chain2_num);
        homo_refined_greedy_search(TMave_mat, assign1_list,
            assign2_list, chain1_num, chain2_num, xcentroids,
            ycentroids, d0MM, len_aa+len_na, ut_mat);
        hetero_refined_greedy_search(TMave_mat, assign1_list,
            assign2_list, chain1_num, chain2_num, xcentroids,
            ycentroids, d0MM, len_aa+len_na);
        
        /* clean up */
        DeleteArray(&xcentroids, chain1_num);
        DeleteArray(&ycentroids, chain2_num);
    }

    /* store initial assignment */
    int init_pair_num=count_assign_pair(assign1_list,chain1_num);
    int *assign1_init, *assign2_init;
    assign1_init=new int[chain1_num];
    assign2_init=new int[chain2_num];
    double **TMave_init;
    NewArray(&TMave_init,chain1_num,chain2_num);
    vector<vector<string> >seqxA_init(chain1_num,tmp_str_vec);
    vector<vector<string> >seqyA_init(chain1_num,tmp_str_vec);
    vector<string> sequence_init;
    copy_chain_assign_data(chain1_num, chain2_num, sequence_init,
        seqxA_mat,  seqyA_mat,  assign1_list, assign2_list, TMave_mat,
        seqxA_init, seqyA_init, assign1_init, assign2_init, TMave_init);

    /* perform iterative alignment */
    double max_total_score=0; // ignore old total_score because previous
                              // score was from monomeric chain superpositions
    int max_iter=5-(int)((len_aa+len_na)/200);
    if (max_iter<2) max_iter=2;
    MMalign_iter(max_total_score, max_iter, xa_vec, ya_vec,
        seqx_vec, seqy_vec, secx_vec, secy_vec, mol_vec1, mol_vec2, xlen_vec,
        ylen_vec, xa, ya, seqx, seqy, secx, secy, len_aa, len_na, chain1_num,
        chain2_num, TMave_mat, seqxA_mat, seqyA_mat, assign1_list, assign2_list,
        sequence, d0_scale, fast_opt, chainmap);

    if (aln_chain_num>=4 && is_oligomer && chainmap.size()==0) // oligomer alignment
    {
        MMalign_final(xname.substr(dir1_opt.size()), yname.substr(dir2_opt.size()),
            chainID_list1, chainID_list2,
            fname_super, fname_lign, fname_matrix,
            xa_vec, ya_vec, seqx_vec, seqy_vec,
            secx_vec, secy_vec, mol_vec1, mol_vec2, xlen_vec, ylen_vec,
            xa, ya, seqx, seqy, secx, secy, len_aa, len_na,
            chain1_num, chain2_num, TMave_mat,
            seqxA_mat, seqM_mat, seqyA_mat, assign1_list, assign2_list, sequence,
            d0_scale, 1, 0, 5, ter_opt, split_opt,
            0, 0, true, true, mirror_opt, resi_vec1, resi_vec2);


        /* extract centroid coordinates */
        double **xcentroids;
        double **ycentroids;
        NewArray(&xcentroids, chain1_num, 3);
        NewArray(&ycentroids, chain2_num, 3);
        double d0MM=getmin(
            calculate_centroids(xa_vec, chain1_num, xcentroids),
            calculate_centroids(ya_vec, chain2_num, ycentroids));

        /* refine enhanced greedy search with centroid superposition */
        //double het_deg=check_heterooligomer(TMave_mat, chain1_num, chain2_num);
        homo_refined_greedy_search(TMave_mat, assign1_list,
            assign2_list, chain1_num, chain2_num, xcentroids,
            ycentroids, d0MM, len_aa+len_na, ut_mat);

        hetero_refined_greedy_search(TMave_mat, assign1_list,
            assign2_list, chain1_num, chain2_num, xcentroids,
            ycentroids, d0MM, len_aa+len_na);

        /* clean up */
        DeleteArray(&xcentroids, chain1_num);
        DeleteArray(&ycentroids, chain2_num);
    }

    /* sometime MMalign_iter is even worse than monomer alignment */
    if (max_total_score<maxTMmono)
    {
        copy_chain_assign_data(chain1_num, chain2_num, sequence,
            seqxA_init, seqyA_init, assign1_init, assign2_init, TMave_init,
            seqxA_mat, seqyA_mat, assign1_list, assign2_list, TMave_mat);
        for (i=0;i<chain1_num;i++)
        {
            if (i!=maxTMmono_i) assign1_list[i]=-1;
            else assign1_list[i]=maxTMmono_j;
        }
        for (j=0;j<chain2_num;j++)
        {
            if (j!=maxTMmono_j) assign2_list[j]=-1;
            else assign2_list[j]=maxTMmono_i;
        }
        sequence[0]=seqxA_mat[maxTMmono_i][maxTMmono_j];
        sequence[1]=seqyA_mat[maxTMmono_i][maxTMmono_j];
        max_total_score=maxTMmono;
        MMalign_iter(max_total_score, max_iter, xa_vec, ya_vec, seqx_vec, seqy_vec,
            secx_vec, secy_vec, mol_vec1, mol_vec2, xlen_vec, ylen_vec,
            xa, ya, seqx, seqy, secx, secy, len_aa, len_na, chain1_num, chain2_num,
            TMave_mat, seqxA_mat, seqyA_mat, assign1_list, assign2_list, sequence,
            d0_scale, fast_opt, chainmap);
    }

    /* perform cross chain alignment
     * in some cases, this leads to dramatic improvement, esp for homodimer */
    int iter_pair_num=count_assign_pair(assign1_list,chain1_num);
    if (iter_pair_num>=init_pair_num) copy_chain_assign_data(
        chain1_num, chain2_num, sequence_init,
        seqxA_mat, seqyA_mat, assign1_list, assign2_list, TMave_mat,
        seqxA_init, seqyA_init, assign1_init,  assign2_init,  TMave_init);
    double max_total_score_cross=max_total_score;

    //if (init_pair_num!=2 && is_oligomer==false) MMalign_cross(
        //max_total_score_cross, max_iter, xa_vec, ya_vec, seqx_vec, seqy_vec,
        //secx_vec, secy_vec, mol_vec1, mol_vec2, xlen_vec, ylen_vec,
        //xa, ya, seqx, seqy, secx, secy, len_aa, len_na, chain1_num, chain2_num,
        //TMave_init, seqxA_init, seqyA_init, assign1_init, assign2_init, sequence_init,
        //d0_scale, true);
    //else 
    if (len_aa+len_na<10000)
    {
        MMalign_dimer(max_total_score_cross, xa_vec, ya_vec, seqx_vec, seqy_vec,
            secx_vec, secy_vec, mol_vec1, mol_vec2, xlen_vec, ylen_vec,
            xa, ya, seqx, seqy, secx, secy, len_aa, len_na, chain1_num, chain2_num,
            TMave_init, seqxA_init, seqyA_init, assign1_init, assign2_init,
            sequence_init, d0_scale, fast_opt);
        if (max_total_score_cross>max_total_score) 
        {
            max_total_score=max_total_score_cross;
            copy_chain_assign_data(chain1_num, chain2_num, sequence,
                seqxA_init, seqyA_init, assign1_init, assign2_init, TMave_init,
                seqxA_mat,  seqyA_mat,  assign1_list, assign2_list, TMave_mat);
        }
    } 

    /* final alignment */
    if (outfmt_opt==0) print_version();
    MMalign_final(xname.substr(dir1_opt.size()), yname.substr(dir2_opt.size()),
        chainID_list1, chainID_list2,
        fname_super, fname_lign, fname_matrix,
        xa_vec, ya_vec, seqx_vec, seqy_vec,
        secx_vec, secy_vec, mol_vec1, mol_vec2, xlen_vec, ylen_vec,
        xa, ya, seqx, seqy, secx, secy, len_aa, len_na,
        chain1_num, chain2_num, TMave_mat,
        seqxA_mat, seqM_mat, seqyA_mat, assign1_list, assign2_list, sequence,
        d0_scale, m_opt, o_opt, outfmt_opt, ter_opt, split_opt,
        a_opt, d_opt, fast_opt, full_opt, mirror_opt, resi_vec1, resi_vec2);

    /* clean up everything */
    delete [] assign1_list;
    delete [] assign2_list;
    DeleteArray(&TMave_mat,chain1_num);
    DeleteArray(&ut_mat,   chain1_num*chain2_num);
    vector<vector<string> >().swap(seqxA_mat);
    vector<vector<string> >().swap(seqM_mat);
    vector<vector<string> >().swap(seqyA_mat);
    vector<string>().swap(tmp_str_vec);

    delete [] assign1_init;
    delete [] assign2_init;
    DeleteArray(&TMave_init,chain1_num);
    vector<vector<string> >().swap(seqxA_init);
    vector<vector<string> >().swap(seqyA_init);

    vector<vector<vector<double> > >().swap(xa_vec); // structure of complex1
    vector<vector<vector<double> > >().swap(ya_vec); // structure of complex2
    vector<vector<char> >().swap(seqx_vec); // sequence of complex1
    vector<vector<char> >().swap(seqy_vec); // sequence of complex2
    vector<vector<char> >().swap(secx_vec); // secondary structure of complex1
    vector<vector<char> >().swap(secy_vec); // secondary structure of complex2
    mol_vec1.clear();       // molecule type of complex1, RNA if >0
    mol_vec2.clear();       // molecule type of complex2, RNA if >0
    vector<string>().swap(chainID_list1);  // list of chainID1
    vector<string>().swap(chainID_list2);  // list of chainID2
    xlen_vec.clear();       // length of complex1
    ylen_vec.clear();       // length of complex2
    vector<string>().swap(chain1_list);
    vector<string>().swap(chain2_list);
    vector<string>().swap(sequence);
    vector<string>().swap(resi_vec1);  // residue index for chain1
    vector<string>().swap(resi_vec2);  // residue index for chain2
    vector<string>().swap(chain2parse1);
    vector<string>().swap(chain2parse2);
    vector<string>().swap(model2parse1);
    vector<string>().swap(model2parse2);

    t2 = clock();
    float diff = ((float)t2 - (float)t1)/CLOCKS_PER_SEC;
    printf("#Total CPU time is %5.2f seconds\n", diff);
    return 0;
}
