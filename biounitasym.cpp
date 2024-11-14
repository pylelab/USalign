#include "basic_fun.h"

using namespace std;

void print_help()
{
    cout <<
"check if the input PDB format multmodel structure is a biounit\n"
"(i.e., biological assembly) or asymmetric unit\n"
"\n"
"Usage: biounitasym pdb.pdb\n"
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
"    -mol     Type of molecule(s) to align.\n"
"             auto: (default) align both protein and nucleic acids.\n"
"             prot: only align proteins in a structure.\n"
"             RNA : only align RNA and DNA in a structure.\n"
"\n"
"    -het     Whether to read residues marked as 'HETATM' in addition to 'ATOM  '\n"
"             0: (default) only align 'ATOM  ' residues\n"
"             1: align both 'ATOM  ' and 'HETATM' residues\n"
"\n"
"    -infmt   Input format for chain\n"
"            -1: (default) automatically detect PDB or PDBx/mmCIF format\n"
"             0: PDB format\n"
"             3: PDBx/mmCIF format\n"
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
    int    ter_opt   =0;     // all models
    int    infmt_opt =-1;    // PDB or PDBx/mmCIF format
    int    split_opt =1;     // do not split chain
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
        if ( !strcmp(argv[i],"-atom") && i < (argc-1) )
        {
            atom_opt=argv[i + 1]; i++;
        }
        else if ( !strcmp(argv[i],"-mol") )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -mol");
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
        else xname=argv[i];
    }

    if(xname.size()==0||xname=="-h") print_help();

    if (suffix_opt.size() && dir_opt.size()==0)
        PrintErrorAndQuit("-suffix is only valid if -dir is set");
    if (atom_opt.size()!=4)
        PrintErrorAndQuit("ERROR! Atom name must have 4 characters, including space.");
    if (mol_opt=="prot") mol_opt="protein";
    else if (mol_opt=="DNA") mol_opt="RNA";
    if (mol_opt!="auto" && mol_opt!="protein" && mol_opt!="RNA")
        PrintErrorAndQuit("ERROR! Molecule type must be one of the"
            "following:\nauto, prot (the same as 'protein'), and "
            "RNA (the same as 'DNA').");
    
    bool autojustify=(atom_opt=="auto" || atom_opt=="PC4'"); // auto re-pad atom name
    if (mol_opt=="protein" && atom_opt=="auto")
        atom_opt=" CA ";
    else if (mol_opt=="RNA" && atom_opt=="auto")
        atom_opt=" C3'";

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
    int    chain_i,chain_j;           // chain index
    int    xlen,ylen;                 // chain length
    int    chainnum;       // number of chains in a PDB file
    char   *seqx, *seqy;       // for the protein sequence 
    double **xa, **ya;         // for input vectors xa[0...xlen-1][0..2] and
                               // ya[0...ylen-1][0..2], in general,
                               // ya is regarded as native structure 
                               // --> superpose xa onto ya
    vector<string> resi_vec1;  // residue index for chain1
    vector<string> resi_vec2;  // residue index for chain2
    vector<double> clashratio_vec;
    double clashcount=0;
    int r1,r2;
    double d2;

    /* loop over file names */
    for (i=0;i<chain_list.size();i++)
    {
        xname=chain_list[i];
        chainnum=get_PDB_lines(xname, PDB_lines, chainID_list, mol_vec,
            ter_opt, infmt_opt, atom_opt, autojustify, split_opt, het_opt,
            chain2parse, model2parse);
        if (!chainnum)
        {
            cerr<<"Warning! Cannot parse file: "<<xname
                <<". Chain number 0."<<endl;
            continue;
        }
        if (chainnum<=1)
        {
            cout<<xname<<"\tsingle model file, use -ter 1 for oligomer alignment"<<endl;
            xname.clear();
            PDB_lines.clear();
            mol_vec.clear();
            clashratio_vec.clear();
            continue;
        }
        for (chain_i=0;chain_i<chainnum;chain_i++)
        {
            xlen=PDB_lines[chain_i].size();
            if (!xlen)
            {
                cerr<<"Warning! Cannot parse file: "<<xname
                    <<". Chain length 0."<<endl;
                continue;
            }
            NewArray(&xa, xlen, 3);
            seqx = new char[xlen + 1];
            xlen = read_PDB(PDB_lines[chain_i], xa, seqx, resi_vec1, 0);
            for (chain_j=chain_i+1;chain_j<chainnum;chain_j++)
            {
                ylen=PDB_lines[chain_j].size();
                if (!ylen)
                {
                    cerr<<"Warning! Cannot parse file: "<<xname
                        <<". Chain length 0."<<endl;
                    continue;
                }
                NewArray(&ya, ylen, 3);
                seqy = new char[ylen + 1];
                ylen = read_PDB(PDB_lines[chain_j], ya, seqy, resi_vec2, 0);
                clashcount=0;
                for (r1=0;r1<xlen;r1++)
                {
                    for (r2=0;r2<ylen;r2++)
                    {
                        d2=dist(xa[r1],ya[r2]);
                        clashcount+=d2<16;
                    }
                }
                if (xlen<=ylen) clashratio_vec.push_back(clashcount/xlen);
                else            clashratio_vec.push_back(clashcount/ylen);
            
                DeleteArray(&ya, ylen);
                delete [] seqy;
                vector<string>().swap(resi_vec2);
            }
            DeleteArray(&xa, xlen);
            delete [] seqx;
            vector<string>().swap(resi_vec1);
        } // chain_i
        clashcount=0;
        for (chain_i=0;chain_i<clashratio_vec.size();chain_i++)
            clashcount+=clashratio_vec[chain_i];
        if (clashratio_vec.size()) clashcount/=clashratio_vec.size();
        if (clashcount>=0.1)
            cout<<xname<<"\tmultimodel asymmetric unit, use -ter 1 for oligomer alignment\t"
                <<"portion of residues with inter-model clash="
                <<setiosflags(ios::fixed)<<setprecision(3)
                <<clashcount<<endl;
        else
            cout<<xname<<"\tbiounit split over multiple model, use -ter 0 for oligomer alignment\t"
                <<"portion of residues with inter-model clash="
                <<setiosflags(ios::fixed)<<setprecision(3)
                <<clashcount<<endl;
        
        xname.clear();
        PDB_lines.clear();
        mol_vec.clear();
        clashratio_vec.clear();
    } // i
    chain_list.clear();
    return 0;
}
