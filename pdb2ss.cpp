#include "TMalign.h"

using namespace std;

// secondary structure    01234
//const char* SSmapProtein=" CHTE";
//const char* SSmapRNA    =" .<>";

void print_help()
{
    cout <<
"Converting PDB file(s) into FASTA format secondary structure sequence.\n"
"Proteins have four states: H E C T (helix, strand, coil, turn)\n"
"RNA have three states: < > . (paired with 3', paired with 5', unpaired)\n"
"\n"
"Usage: pdb2ss pdb.pdb > seq.ss\n"
"\n"
"    -dir     Convert all chains listed by 'chain_list' under 'chain_folder'.\n"
"             Note that the slash is necessary.\n"
"             $ pdb2xyz -dir chain_folder/ chain_list\n"
"\n"
"    -suffix  (Only when -dir is set, default is empty)\n"
"             add file name suffix to files listed by chain_list\n"
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
    int    split_opt =0;     // do not split chain
    int    het_opt=0;        // do not read HETATM residues
    string atom_opt  ="auto";// use C alpha atom for protein and C3' for RNA
    string mol_opt   ="auto";// auto-detect the molecule type as protein/RNA
    string suffix_opt="";    // set -suffix to empty
    string dir_opt   ="";    // set -dir to empty
    vector<string> chain_list; // only when -dir1 is set
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
        else if ( !strcmp(argv[i],"-mol") && i < (argc-1) )
        {
            mol_opt=argv[i + 1]; i++;
        }
        else if ( !strcmp(argv[i],"-dir") && i < (argc-1) )
        {
            dir_opt=argv[i + 1]; i++;
        }
        else if ( !strcmp(argv[i],"-suffix") && i < (argc-1) )
        {
            suffix_opt=argv[i + 1]; i++;
        }
        else if ( !strcmp(argv[i],"-infmt") && i < (argc-1) )
        {
            infmt_opt=atoi(argv[i + 1]); i++;
        }
        else if ( !strcmp(argv[i],"-het") && i < (argc-1) )
        {
            het_opt=atoi(argv[i + 1]); i++;
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

    if (suffix_opt.size() && dir_opt.size()==0)
        PrintErrorAndQuit("-suffix is only valid if -dir is set");

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
    if (dir_opt.size()==0)
        chain_list.push_back(xname);
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
            chain_list.push_back(dir_opt+Trim(line)+suffix_opt);
        }
        fp.close();
        line.clear();
    }

    /* declare previously global variables */
    vector<vector<string> >PDB_lines; // text of chain
    vector<int> mol_vec;              // molecule type of chain
    vector<string> chainID_list;      // list of chainID1
    int    i;                         // file index
    int    l;                         // residue index
    int    chain_i;                   // chain index
    int    xlen;                      // chain length
    int    xchainnum;                 // number of chains in a PDB file
    char   *seqx;                     // for the protein sequence 
    char   *secx;                     // for the secondary structure 
    double **xa;                      // for input vectors xa[0...xlen-1][0..2] and
    vector<string> resi_vec;          // residue index for chain

    /* loop over file names */
    for (i=0;i<chain_list.size();i++)
    {
        xname=chain_list[i];
        xchainnum=get_PDB_lines(xname, PDB_lines, chainID_list, mol_vec,
            ter_opt, infmt_opt, atom_opt, autojustify, split_opt, het_opt,
            chain2parse, model2parse);
        if (!xchainnum)
        {
            cerr<<"Warning! Cannot parse file: "<<xname
                <<". Chain number 0."<<endl;
            continue;
        }
        for (chain_i=0;chain_i<xchainnum;chain_i++)
        {
            xlen=PDB_lines[chain_i].size();
            if (mol_opt=="RNA") mol_vec[chain_i]=1;
            else if (mol_opt=="protein") mol_vec[chain_i]=-1;
            if (!xlen)
            {
                cerr<<"Warning! Cannot parse file: "<<xname
                    <<". Chain length 0."<<endl;
                continue;
            }
            NewArray(&xa, xlen, 3);
            seqx = new char[xlen + 1];
            secx = new char[xlen + 1];
            xlen = read_PDB(PDB_lines[chain_i], xa, seqx, resi_vec, 0);
            if (mol_vec[chain_i]>0) make_sec(seqx,xa, xlen, secx,atom_opt);
            else make_sec(xa, xlen, secx); // protein
            
            cout<<'>'<<xname.substr(dir_opt.size(),
                xname.size()-dir_opt.size()-suffix_opt.size())
                <<chainID_list[chain_i]<<'\t'<<xlen<<'\n'<<secx<<endl;

            PDB_lines[chain_i].clear();
            DeleteArray(&xa, xlen);
            delete [] seqx;
            delete [] secx;
        } // chain_i
        xname.clear();
        PDB_lines.clear();
        resi_vec.clear();
        mol_vec.clear();
    } // i
    chain_list.clear();
    vector<string>().swap(chain2parse);
    vector<string>().swap(model2parse);
    return 0;
}
