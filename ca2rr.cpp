#include "ca2rr.h"

using namespace std;

void print_help()
{
    cout <<
"Calculate C beta atom contact map from C alpha atoms in PDB file(s).\n"
"\n"
"Usage: ca2rr pdb.pdb > pdb.rr\n"
"\n"
"    -atom    4-character atom name used to represent a residue.\n"
"             \" CA \": (default) infer C beta contact from C alpha atoms\n"
"             \" CB \": calculate C beta contact\n"
"                       (note the spaces before and after CA and CB).\n"
"\n"
"    -ter     Strings to mark the end of a chain\n"
"             3: (default) TER, ENDMDL, END or different chain ID\n"
"             2: ENDMDL, END, or different chain ID\n"
"             1: ENDMDL or END\n"
"             0: end of file\n"
"\n"
"    -split   Whether to split PDB file into multiple chains\n"
"             0: (default) treat the whole structure as one single chain\n"
"             1: treat each MODEL as a separate chain (-ter should be 0)\n"
"             2: treat each chain as a seperate chain (-ter should be <=1)\n"
"\n"
"    -infmt   Input format for chain\n"
"            -1: (default) automatically detect PDB or PDBx/mmCIF format\n"
"             0: PDB format\n"
"             2: xyz format\n"
"             3: PDBx/mmCIF format\n"
"\n"
"    -outfmt  Output format\n"
"             0: PDB format CA and predicted CB atoms of full length structure\n"
"             1: (default) CASP RR format contact map\n"
"             2: FU-score for contact-based domain partition\n"
"             3: Domain structure\n"
"\n"
"    -sep     Minimal sequence separation. default: 6. must be non-negative\n"
"\n"
"    -mdl     Minimal domain length. default: 30\n"
"\n"
"    -hinge   Maximum number of domains allowed. default: 9\n"
"\n"
"    -het     Whether to read residues marked as 'HETATM' in addition to 'ATOM  '\n"
"             0: (default) only align 'ATOM  ' residues\n"
"             1: align both 'ATOM  ' and 'HETATM' residues\n"
"\n"
    <<endl;
    exit(EXIT_SUCCESS);
}

int main(int argc, char *argv[])
{
    if (argc < 2) print_help();


    /**********************/
    /*    get argument    */
    /**********************/
    string xname     = "";
    int    ter_opt   =3;     // TER, END, or different chainID
    int    infmt_opt =-1;    // PDB format
    int    outfmt_opt=1;     // set -outfmt to CASP RR output
    int    sep_opt   =6;     // minimal sequence separation
    int    mdl_opt   =30;    // minimal domain length
    int    split_opt =0;     // do not split chain
    int    het_opt   =0;     // do not read HETATM residues
    int    hinge_opt =9;     // maximum number of hinge allowed for flexible
    string atom_opt  =" CA ";// use C alpha atom for protein and C3' for RNA
    string mol_opt   ="protein";// auto-detect the molecule type as protein/RNA
    vector<string> chain2parse;
    vector<string> model2parse;

    int nameIdx = 0;
    for(int i = 1; i < argc; i++)
    {
        if ( !strcmp(argv[i],"-ter") && i < (argc-1) )
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
        else if ( !strcmp(argv[i],"-infmt") && i < (argc-1) )
        {
            infmt_opt=atoi(argv[i + 1]); i++;
        }
        else if ( !strcmp(argv[i],"-outfmt") )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -outfmt");
            outfmt_opt=atoi(argv[i + 1]); i++;
        }
        else if ( !strcmp(argv[i],"-sep") )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -sep");
            sep_opt=atoi(argv[i + 1]); i++;
        }
        else if ( !strcmp(argv[i],"-mdl") )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -mdl");
            mdl_opt=atoi(argv[i + 1]); i++;
        }
        else if ( !strcmp(argv[i],"-het") && i < (argc-1) )
        {
            het_opt=atoi(argv[i + 1]); i++;
        }
        else if ( !strcmp(argv[i],"-hinge") )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -hinge");
            hinge_opt = atoi(argv[i + 1]); i++;
        }
        else if (!strcmp(argv[i], "-chain") )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -chain");
            split(argv[i+1],chain2parse,',');
            i++;
        }
        else if (!strcmp(argv[i], "-model") )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -model");
            split(argv[i+1],model2parse,',');
            i++;
        }
        else xname=argv[i];
    }

    if(xname.size()==0||xname=="-h") print_help();

    if (!(atom_opt==" CA " || atom_opt==" CB "))
        PrintErrorAndQuit("ERROR! Atom name must be \" CA \" or \" CB \"");
    if (split_opt==1 && ter_opt!=0)
        PrintErrorAndQuit("-split 1 should be used with -ter 0");
    else if (split_opt==2 && ter_opt!=0 && ter_opt!=1)
        PrintErrorAndQuit("-split 2 should be used with -ter 0 or 1");
    if (split_opt<0 || split_opt>2)
        PrintErrorAndQuit("-split can only be 0, 1 or 2");
    if (sep_opt<0)
        PrintErrorAndQuit("-sep must be non-negative");
    if (mdl_opt<3)
        PrintErrorAndQuit("-mdl must be >=3");
    if (hinge_opt<=1)// || hinge_opt>10)
        PrintErrorAndQuit("ERROR! -hinge must be >=2");

    /* declare previously global variables */
    vector<vector<string> >PDB_lines; // text of chain
    vector<int> mol_vec;              // molecule type of chain
    vector<string> chainID_list;      // list of chainID1
    int    i;                         // file index
    int    r,r1,r2;                   // residue index
    int    chain_i;                   // chain index
    int    xlen;                      // chain length
    int    xchainnum;                 // number of chains in a PDB file
    char   *seqx;                     // for the protein sequence 
    char   *secx;                     // for the secondary structure 
    double **xa;                      // CA atom
    double **ya;                      // CB atom
    bool   **ct;                      // contact map
    vector<string> resi_vec;          // residue index for chain
    int    l;
    long   N1,N2,N12;

    /* loop over file names */
    xchainnum=get_PDB_lines(xname, PDB_lines, chainID_list, mol_vec,
        ter_opt, infmt_opt, atom_opt, false, split_opt, het_opt,
        chain2parse, model2parse);
    if (!xchainnum)
    {
        cerr<<"Warning! Cannot parse file: "<<xname
            <<". Chain number 0."<<endl;
    }
    for (chain_i=0;chain_i<xchainnum;chain_i++)
    {
        xlen=PDB_lines[chain_i].size();
        if (!xlen)
        {
            cerr<<"Warning! Cannot parse file: "<<xname
                <<". Chain length 0."<<endl;
            continue;
        }
        if (xlen<=1 && outfmt_opt>=1) continue;
        NewArray(&xa, xlen, 3);
        NewArray(&ct, xlen, xlen);
        seqx = new char[xlen + 1];
        secx = new char[xlen + 1];
        xlen = read_PDB(PDB_lines[chain_i], xa, seqx, resi_vec, 0);
        if (outfmt_opt>=2)
        {
            for (r1=0;r1<xlen;r1++)
                for (r2=0;r2<xlen;r2++)
                    ct[r1][r2]=0;
        }

        if (atom_opt==" CB ")
        {
            if (outfmt_opt==0)
            {
                for (r=0;r<xlen;r++)
                    cout<<PDB_lines[chain_i][r]<<endl;
            }
            else if (outfmt_opt==1)
            {
                for (r1=0;r1<xlen;r1++)
                    for (r2=r1+sep_opt;r2<xlen;r2++)
                        if (dist(xa[r1],xa[r2])<=64) // dij<=8
                            cout<<Trim(PDB_lines[chain_i][r1].substr(22,5))
                                <<(split_opt?chainID_list[chain_i]:string())
                                <<' '<<Trim(PDB_lines[chain_i][r2].substr(22,5))
                                <<(split_opt?chainID_list[chain_i]:string())
                                <<" 0 8 1"<<endl;
            }
            else if (outfmt_opt==2)
            {
                for (r1=0;r1<xlen;r1++)
                    for (r2=r1+sep_opt;r2<xlen;r2++)
                        ct[r2][r1]=ct[r1][r2]=(dist(xa[r1],xa[r2])<=64);
            }
        }
        else
        {
            NewArray(&ya, xlen, 3); // CB atom
            make_sec(xa, xlen, secx); // protein
            ca2cb(xa, seqx, secx, xlen, ya);
            
            if (outfmt_opt==0)
            {
                for (r=0;r<xlen;r++)
                {
                    cout<<PDB_lines[chain_i][r]<<endl;
                    if (seqx[r]=='G') continue;
                    cout<<PDB_lines[chain_i][r].substr(0,12)<<" CB "
                        <<PDB_lines[chain_i][r].substr(16,14)
                        <<setiosflags(ios::fixed)<<setprecision(3)
                        <<setw(8)<<ya[r][0]
                        <<setw(8)<<ya[r][1]
                        <<setw(8)<<ya[r][2]
                        <<PDB_lines[chain_i][r].substr(54)<<endl;
                }
            }
            else if (outfmt_opt==1)
            {
                for (r1=0;r1<xlen;r1++)
                    for (r2=r1+sep_opt;r2<xlen;r2++)
                        if (dist(ya[r1],ya[r2])<=64) // dij<=8
                            cout<<Trim(PDB_lines[chain_i][r1].substr(22,5))
                                <<(split_opt?chainID_list[chain_i]:string())
                                <<' '<<Trim(PDB_lines[chain_i][r2].substr(22,5))
                                <<(split_opt?chainID_list[chain_i]:string())
                                <<" 0 8 1"<<endl;
            }
            else if (outfmt_opt>=2)
            {
                for (r1=0;r1<xlen;r1++)
                    for (r2=r1+sep_opt;r2<xlen;r2++)
                        ct[r2][r1]=ct[r1][r2]=(dist(ya[r1],ya[r2])<=64);
            }
        }

        if (outfmt_opt>=2)
        {
            vector<int> l_vec;
            iterative_calFUscore(ct, l_vec, xlen, hinge_opt, mdl_opt, true);
            if (outfmt_opt==3)
            {
                if (l_vec.size()==0)
                    cout<<"#cannot partition "<<xname<<endl;
                else
                {
                    int d_start=0;
                    int d_end  =xlen;
                    for (r1=0;r1<=l_vec.size();r1++)
                    {
                        d_start=(r1==0)?0:l_vec[r1-1];
                        d_end  =(r1==l_vec.size())?xlen:l_vec[r1];
                        string suffix=to_string(r1);
                        string filename=xname+'_'+suffix+".pdb";
                        cout<<'('<<1+d_start<<','<<d_end<<") "<<filename<<endl;
                        ofstream fout(filename);
                        for (r=d_start;r<d_end;r++)
                        {
                            fout<<PDB_lines[chain_i][r]<<endl;
                            if (atom_opt==" CA ") fout
                                <<PDB_lines[chain_i][r].substr(0,12)<<" CB "
                                <<PDB_lines[chain_i][r].substr(16,14)
                                <<setiosflags(ios::fixed)<<setprecision(3)
                                <<setw(8)<<ya[r][0]
                                <<setw(8)<<ya[r][1]
                                <<setw(8)<<ya[r][2]
                                <<PDB_lines[chain_i][r].substr(54)<<endl;

                        }
                        fout.close();
                    }
                }
            }
        }
        if (atom_opt==" CA ") DeleteArray(&ya, xlen);
        DeleteArray(&ct, xlen);
        PDB_lines[chain_i].clear();
        DeleteArray(&xa, xlen);
        delete [] seqx;
        delete [] secx;
    } // chain_i
    xname.clear();
    PDB_lines.clear();
    resi_vec.clear();
    mol_vec.clear();
    vector<string>().swap(chain2parse);
    vector<string>().swap(model2parse);
    return 0;
}
