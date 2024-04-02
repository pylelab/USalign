#include "NWalign.h"

using namespace std;

void print_extra_help()
{
    cout <<
"Additional options:\n"
"    -dir     Perform all-against-all alignment among the list of PDB\n"
"             chains listed by 'chain_list' under 'chain_folder'. Note\n"
"             that the slash is necessary.\n"
"             $ NWalign -dir chain_folder/ chain_list\n"
"\n"
"    -dir1    Use chain2 to search a list of PDB chains listed by 'chain1_list'\n"
"             under 'chain1_folder'. Note that the slash is necessary.\n"
"             $ NWalign -dir1 chain1_folder/ chain1_list chain2\n"
"\n"
"    -dir2    Use chain1 to search a list of PDB chains listed by 'chain2_list'\n"
"             under 'chain2_folder'\n"
"             $ NWalign chain1 -dir2 chain2_folder/ chain2_list\n"
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
"             For FASTA intput (-infmt1/-infmt2 4), -ter 0 means read all\n"
"             sequences; -ter >=1 means read the first sequence only."
"\n"
"    -split   Whether to split PDB file into multiple chains\n"
"             0: (default) treat the whole structure as one single chain\n"
"             1: treat each MODEL as a separate chain (-ter should be 0)\n"
"             2: treat each chain as a seperate chain (-ter should be <=1)\n"
"\n"
"             For FASTA intput, -split 0 means concatenate all sequences into\n"
"             one read all sequence; -split >=1 means each sequence is an\n"
"             individual entry."
"\n"
"    -het     Whether to align residues marked as 'HETATM' in addition to 'ATOM  '\n"
"             0: (default) only align 'ATOM  ' residues\n"
"             1: align both 'ATOM  ' and 'HETATM' residues\n"
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
"Pairwise sequence alignment between two sequences.\n"
"\n"
"Usage: NWalign PDB1.pdb PDB2.pdb [Options]\n"
"\n"
"Options:\n"
"    -h       Print the full help message\n"
"\n"
"    -glocal  Global or local alignment\n"
"             0: (default) Needleman-Wunsch algorithm for global alignment\n"
"             1: glocal-query alignment\n"
"             2: glocal-both alignment\n"
"             3: Smith-Waterman algorithm for local alignment\n"
"\n"
"    -infmt1  Input format for chain1\n"
"    -infmt2  Input format for chain2\n"
"            -1: (default) automatically detect PDB or PDBx/mmCIF format\n"
"             0: PDB format\n"
"             2: xyz format\n"
"             3: PDBx/mmCIF format\n"
"             4: FASTA format sequence\n"
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
    string xname     ="";
    string yname     ="";
    bool   h_opt     =false; // print full help message
    int    infmt1_opt=-1;    // FASTA sequence
    int    infmt2_opt=-1;    // FASTA sequence
    int    ter_opt   =3;     // TER, END, or different chainID
    int    split_opt =0;     // do not split chain
    int    outfmt_opt=0;     // set -outfmt to full output
    int    het_opt=0;        // do not read HETATM residues
    string atom_opt  ="auto";// use C alpha atom for protein and C3' for RNA
    string mol_opt   ="auto";// auto-detect the molecule type as protein/RNA
    string suffix_opt="";    // set -suffix to empty
    string dir_opt   ="";    // set -dir to empty
    string dir1_opt  ="";    // set -dir1 to empty
    string dir2_opt  ="";    // set -dir2 to empty
    vector<string> chain1_list; // only when -dir1 is set
    vector<string> chain2_list; // only when -dir2 is set
    vector<string> chain2parse1;
    vector<string> chain2parse2;
    vector<string> model2parse1;
    vector<string> model2parse2;
    int    glocal    =0;

    for(int i = 1; i < argc; i++)
    {
        if ( !strcmp(argv[i],"-h") )
        {
            h_opt = true;
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
        else if ( !strcmp(argv[i],"-glocal") && i < (argc-1) )
        {
            glocal=atoi(argv[i + 1]); i++;
        }
        else if ( !strcmp(argv[i],"-het") && i < (argc-1) )
        {
            het_opt=atoi(argv[i + 1]); i++;
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
        else if (xname.size() == 0) xname=argv[i];
        else if (yname.size() == 0) yname=argv[i];
        else PrintErrorAndQuit(string("ERROR! Undefined option ")+argv[i]);
    }

    if(xname.size()==0 || (yname.size()==0 && dir_opt.size()==0) || 
                          (yname.size()    && dir_opt.size()))
    {
        if (h_opt) print_help(h_opt);
        if (xname.size()==0)
            PrintErrorAndQuit("Please provide input sequences");
        else if (yname.size()==0 && dir_opt.size()==0)
            PrintErrorAndQuit("Please provide sequence 2");
        else if (yname.size() && dir_opt.size())
            PrintErrorAndQuit("Please provide only one file name if -dir is set");
    }

    if (suffix_opt.size() && dir_opt.size()+dir1_opt.size()+dir2_opt.size()==0)
        PrintErrorAndQuit("-suffix is only valid if -dir, -dir1 or -dir2 is set");
    if (dir_opt.size() && (dir1_opt.size() || dir2_opt.size()))
        PrintErrorAndQuit("-dir cannot be set with -dir1 or -dir2");

    bool autojustify=(atom_opt=="auto" || atom_opt=="PC4'"); // auto re-pad atom name
    if (atom_opt.size()!=4)
        PrintErrorAndQuit("ERROR! Atom name must have 4 characters, including space.");
    if (mol_opt!="auto" && mol_opt!="protein" && mol_opt!="RNA")
        PrintErrorAndQuit("ERROR! Molecule type must be either RNA or protein.");
    else if (mol_opt=="protein" && atom_opt=="auto")
        atom_opt=" CA ";
    else if (mol_opt=="RNA" && atom_opt=="auto")
        atom_opt=" C3'";

    if (split_opt==1 && ter_opt!=0)
        PrintErrorAndQuit("-split 1 should be used with -ter 0");
    else if (split_opt==2 && ter_opt!=0 && ter_opt!=1)
        PrintErrorAndQuit("-split 2 should be used with -ter 0 or 1");
    if (split_opt<0 || split_opt>2)
        PrintErrorAndQuit("-split can only be 0, 1 or 2");

    /* parse file list */
    if (dir1_opt.size()+dir_opt.size()==0) chain1_list.push_back(xname);
    else file2chainlist(chain1_list, xname, dir_opt+dir1_opt, suffix_opt);

    if (dir_opt.size())
        for (int i=0;i<chain1_list.size();i++)
            chain2_list.push_back(chain1_list[i]);
    else if (dir2_opt.size()==0) chain2_list.push_back(yname);
    else file2chainlist(chain2_list, yname, dir2_opt, suffix_opt);

    if (outfmt_opt==2)
        cout<<"#sequence1\tsequence2\tID1\tID2\tIDali\tL1\tL2\tLali"<<endl;

    /* declare previously global variables */
    vector<vector<string> >PDB_lines1; // text of chain1
    vector<vector<string> >PDB_lines2; // text of chain2
    vector<int> mol_vec1;              // molecule type of chain1, RNA if >0
    vector<int> mol_vec2;              // molecule type of chain2, RNA if >0
    vector<string> chainID_list1;      // list of chainID1
    vector<string> chainID_list2;      // list of chainID2
    int  i,j;                // file index
    int  chain_i,chain_j;    // chain index
    int  xlen, ylen;         // chain length
    int  xchainnum,ychainnum;// number of chains in a PDB file
    char *seqx, *seqy;       // for the protein sequence 
    int  l;                  // residue index

    /* loop over file names */
    for (i=0;i<chain1_list.size();i++)
    {
        /* parse chain 1 */
        xname=chain1_list[i];
        if (infmt1_opt>=4) xchainnum=get_FASTA_lines(xname, PDB_lines1, 
                chainID_list1, mol_vec1, ter_opt, split_opt);
        else xchainnum=get_PDB_lines(xname, PDB_lines1, chainID_list1, mol_vec1,
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
            if (infmt1_opt>=4) xlen=PDB_lines1[chain_i][0].size();
            else xlen=PDB_lines1[chain_i].size();
            if (mol_opt=="RNA") mol_vec1[chain_i]=1;
            else if (mol_opt=="protein") mol_vec1[chain_i]=-1;
            if (!xlen)
            {
                cerr<<"Warning! Cannot parse file: "<<xname
                    <<". Chain length 0."<<endl;
                continue;
            }
            seqx = new char[xlen + 1];
            if (infmt1_opt>=4) strcpy(seqx,PDB_lines1[chain_i][0].c_str());
            else for (l=0;l<xlen;l++)
                seqx[l]=AAmap(PDB_lines1[chain_i][l].substr(17,3));
            seqx[xlen]=0;
            
            for (j=(dir_opt.size()>0)*(i+1);j<chain2_list.size();j++)
            {
                /* parse chain 2 */
                if (PDB_lines2.size()==0)
                {
                    yname=chain2_list[j];
                    if (infmt2_opt>=4)
                         ychainnum=get_FASTA_lines(yname, PDB_lines2,
                            chainID_list2, mol_vec2, ter_opt, split_opt);
                    else ychainnum=get_PDB_lines(yname, PDB_lines2,
                            chainID_list2, mol_vec2, ter_opt, infmt2_opt,
                            atom_opt, autojustify, split_opt, het_opt,
                            chain2parse2, model2parse2);
                    if (!ychainnum)
                    {
                        cerr<<"Warning! Cannot parse file: "<<yname
                            <<". Chain number 0."<<endl;
                        continue;
                    }
                }
                for (chain_j=0;chain_j<ychainnum;chain_j++)
                {
                    if (infmt2_opt>=4) ylen=PDB_lines2[chain_j][0].size();
                    else ylen=PDB_lines2[chain_j].size();
                    if (mol_opt=="RNA") mol_vec2[chain_j]=1;
                    else if (mol_opt=="protein") mol_vec2[chain_j]=-1;
                    if (!ylen)
                    {
                        cerr<<"Warning! Cannot parse file: "<<yname
                            <<". Chain length 0."<<endl;
                        continue;
                    }
                    seqy = new char[ylen + 1];
                    if (infmt2_opt>=4) 
                        strcpy(seqy,PDB_lines2[chain_j][0].c_str());
                    else for (l=0;l<ylen;l++)
                        seqy[l]=AAmap(PDB_lines2[chain_j][l].substr(17,3));
                    seqy[ylen]=0;

                    int L_ali;                // Aligned length
                    double Liden=0;
                    string seqM, seqxA, seqyA;// for output alignment
                    int *invmap = new int[ylen+1];
                    
                    int aln_score=NWalign_main(seqx, seqy, xlen, ylen,
                        seqxA, seqyA, mol_vec1[chain_i]+mol_vec2[chain_j],
                        invmap, (outfmt_opt>=2)?1:0, glocal);
                    
                    if (outfmt_opt>=2) get_seqID(invmap, seqx, seqy, 
                        ylen, Liden, L_ali);
                    else get_seqID(seqxA, seqyA, seqM, Liden, L_ali);

                    output_NWalign_results(
                        xname.substr(dir1_opt.size()+dir_opt.size()),
                        yname.substr(dir2_opt.size()+dir_opt.size()),
                        chainID_list1[chain_i].c_str(),
                        chainID_list2[chain_j].c_str(),
                        xlen, ylen, seqM.c_str(), seqxA.c_str(),
                        seqyA.c_str(), Liden, L_ali, aln_score, outfmt_opt);

                    /* Done! Free memory */
                    seqM.clear();
                    seqxA.clear();
                    seqyA.clear();
                    delete [] seqy;
                    delete [] invmap;
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
            delete [] seqx;
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
        chainID_list2.clear();
        mol_vec2.clear();
    }
    chain1_list.clear();
    chain2_list.clear();
    vector<string>().swap(chain2parse1);
    vector<string>().swap(chain2parse2);
    vector<string>().swap(model2parse1);
    vector<string>().swap(model2parse2);
    return 0;
}
