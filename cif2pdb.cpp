#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <iterator>
#include <algorithm>
#include <string>
#include <iomanip>
#include <map>

#include "pstream.h" // For reading gzip and bz2 compressed files

using namespace std;

void print_help()
{
    cout <<
"Converting mmCIF file to PDB file(s)\n"
"\n"
"Usage: cif2pdb input.cif output.pdb\n"
"\n"
"    -chain   Specify auth chain ID to convert:\n"
"             $ cif2pdb input.cif output.pdb -chain A\n"
"\n"
"    -mol     macromolecule type. default is all macromolecules.\n"
"             1: protein only\n"
"             2: RNA only\n"
"             4: DNA only\n"
"\n"
"    -split   Whether to split PDB file into multiple chains\n"
"             0: (default) do not split; output a single PDB\n"
"             1: output one PDB file per chain\n"
"\n"
"    -het     Whether to read residues marked as 'HETATM' in addition to 'ATOM  '\n"
"             0: only 'ATOM  ' residues\n"
"             1: (default) 'ATOM  ' and 'HETATM' for MSE\n"
"             2: 'ATOM  ' and all 'HETATM', excluding HOH\n"
"             3: 'ATOM  ' and all 'HETATM', including HOH\n"
"             If -het >=1, MSE will be converted to MET\n"
    <<endl;
    exit(EXIT_SUCCESS);
}

void PrintErrorAndQuit(const string sErrorString)
{
    cout << sErrorString << endl;
    exit(1);
}

/* strip white space at the begining or end of string */
string Trim(const string &inputString)
{
    string result = inputString;
    int idxBegin = inputString.find_first_not_of(" \n\r\t");
    int idxEnd = inputString.find_last_not_of(" \n\r\t");
    if (idxBegin >= 0 && idxEnd >= 0)
        result = inputString.substr(idxBegin, idxEnd + 1 - idxBegin);
    return result;
}

/* split a long string into vectors by whitespace 
 * line          - input string
 * line_vec      - output vector 
 * delimiter     - delimiter */
void split(const string &line, vector<string> &line_vec,
    const char delimiter=' ')
{
    bool within_word = false;
    for (size_t pos=0;pos<line.size();pos++)
    {
        if (line[pos]==delimiter)
        {
            within_word = false;
            continue;
        }
        if (!within_word)
        {
            within_word = true;
            line_vec.push_back("");
        }
        line_vec.back()+=line[pos];
    }
}

void write_mmcif_to_pdb(const string filename,
    const vector<vector<string> >&PDB_lines,
    const vector<string> &chainID_list, const int split_opt)
{
    size_t c,r;
    
    ofstream fout;
    if (split_opt)
    {
        for (c=0;c<PDB_lines.size();c++)
        {
            if (PDB_lines[c].size()==0) continue;
            if (filename=="-")
            {
                cout<<"REMARK cif2pdb "<<PDB_lines[c][0][21]<<" "<<chainID_list[c]<<endl;
                for (r=0;r<PDB_lines[c].size();r++) cout<<PDB_lines[c][r];
                cout<<"TER"<<endl;
                continue;    
            }
            cout<<     filename+Trim(chainID_list[c])+".pdb"<<endl;
            fout.open((filename+Trim(chainID_list[c])+".pdb").c_str());
            fout<<"REMARK cif2pdb "<<PDB_lines[c][0][21]<<" "<<chainID_list[c]<<endl;
            for (r=0;r<PDB_lines[c].size();r++) fout<<PDB_lines[c][r];
            fout<<"TER"<<endl;
            fout.close();
        }
    }
    else if (filename=="-")
    {
        for (c=0;c<PDB_lines.size();c++)
        {
            if (PDB_lines[c].size()==0) continue;
            cout<<"REMARK cif2pdb "<<PDB_lines[c][0][21]<<" "<<chainID_list[c]<<endl;
        }
        for (c=0;c<PDB_lines.size();c++)
        {
            if (PDB_lines[c].size()==0) continue;
            for (r=0;r<PDB_lines[c].size();r++) cout<<PDB_lines[c][r];
            cout<<"TER"<<endl;
        }
        cout<<"END"<<endl;
    }
    else
    {
        fout.open(filename.c_str());
        for (c=0;c<PDB_lines.size();c++)
        {
            if (PDB_lines[c].size()==0) continue;
            fout<<"REMARK cif2pdb "<<PDB_lines[c][0][21]<<" "<<chainID_list[c]<<endl;
        }
        for (c=0;c<PDB_lines.size();c++)
        {
            if (PDB_lines[c].size()==0) continue;
            for (r=0;r<PDB_lines[c].size();r++) fout<<PDB_lines[c][r];
            fout<<"TER"<<endl;
        }
        fout<<"END"<<endl;
        fout.close();
    }
    return;
}

size_t resolve_chainID_for_mmcif(vector<vector<string> >&PDB_lines,
    const vector<string> &chainID_list)
{
    size_t changed_chains=0;
    size_t c,r,i;
    string chainID;
    
    string chainID_string="ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789";
    vector<bool> chainID_taken(chainID_string.size(),false);
    vector<bool> chainID_accept(chainID_list.size(),false);

    /* accept all single character chain ID */
    for (c=0;c<PDB_lines.size();c++)
    {
        if (PDB_lines[c].size()==0) continue;
        chainID=PDB_lines[c][0][21];
        if (chainID!=chainID_list[c]) continue;
        chainID_accept[c]=true;
        for (i=0;i<chainID_string.size();i++)
        {
            if (chainID_string[i]!=chainID[0]) continue;
            if (chainID_taken[i]) chainID_accept[c]=false;
            else chainID_taken[i]=true;
            break;
        }
    }

    /* accept all remaining non-conflicting chain ID */
    for (c=0;c<PDB_lines.size();c++)
    {
        if (PDB_lines[c].size()==0 || chainID_accept[c]) continue;
        chainID=PDB_lines[c][0][21];
        chainID_accept[c]=true;
        for (i=0;i<chainID_string.size();i++)
        {
            if (chainID_string[i]!=chainID[0]) continue;
            if (chainID_taken[i]) chainID_accept[c]=false;
            else chainID_taken[i]=true;
            break;
        }
    }

    /* resolve remaining chain ID */
    for (c=0;c<PDB_lines.size();c++)
    {
        if (PDB_lines[c].size()==0 || chainID_accept[c]) continue;
        chainID="";
        for (i=0;i<chainID_taken.size();i++)
        {
            if (chainID_taken[i]) continue;
            chainID=chainID_string[i];
            chainID_taken[i]=true;
            break;
        }
        if (chainID=="")
        {
            cerr<<"WARNING! Cannot parse "<<chainID_list[c]<<" with "
                <<PDB_lines[c].size()<<" atoms due to chain ID conflict. "
                <<"Please consider -split 1"<<endl;
            vector<string>().swap(PDB_lines[c]);
        }
        else
        {
            for (r=0;r<PDB_lines[c].size();r++) PDB_lines[c][r]=
                PDB_lines[c][r].substr(0,21)+chainID+PDB_lines[c][r].substr(22);
            changed_chains++;
        }
    }
    if (changed_chains)
        cerr<<"WARNING! Changed "<<changed_chains<<" chain ID(s)"<<endl;
    
    /* clean up*/
    chainID.clear();
    string().swap(chainID_string);
    vector<bool>().swap(chainID_taken);
    vector<bool>().swap(chainID_accept);
    return changed_chains;
}

size_t get_all_mmcif_lines(const string filename, const string chain_opt,
    vector<vector<string> >&PDB_lines, vector<string> &chainID_list,
    const bool dna_opt, const bool rna_opt, const bool protein_opt,
    const bool hoh_opt,  const bool lig_opt, const bool mse_opt)
{
    size_t a=0; // atom index
    string line;
    bool select_atom=false;
    size_t model_idx=0;
    vector<string> tmp_str_vec;
    
    int compress_type=0; // uncompressed file
    ifstream fin;
    redi::ipstream fin_gz; // if file is compressed
    if (filename.size()>=3 && 
        filename.substr(filename.size()-3,3)==".gz")
    {
        fin_gz.open("gunzip -c '"+filename+"'");
        compress_type=1;
    }
    else if (filename.size()>=4 && 
        filename.substr(filename.size()-4,4)==".bz2")
    {
        fin_gz.open("bzcat '"+filename+"'");
        compress_type=2;
    }
    else
    {
        if (filename=="-") compress_type=-1;
        else fin.open(filename.c_str());
    }

    bool loop_ = false; // not reading following content
    map<string,int> _atom_site;
    int atom_site_pos;
    vector<string> line_vec;
    string group_PDB="ATOM  ";
    string alt_id=" ";  // alternative location indicator
    string asym_id="."; // this is similar to chainID, except that
                        // chainID is char while asym_id is a string
                       // with possibly multiple char
    string prev_asym_id="";
    string resn="";       // residue name
    string resi="";
    string atom="";
    string model_index=""; // the same as model_idx but type is string
    stringstream i8_stream;
    map<string, string> alt_id_dict; // resi -> alt_id
    string resi_chain;
    while ((compress_type==-1)?cin.good():(compress_type?fin_gz.good():fin.good()))
    {
        if  (compress_type==-1) getline(cin, line);
        else if (compress_type) getline(fin_gz, line);
        else                    getline(fin, line);
        if (line.size()==0) continue;
        if (loop_) loop_ = (line.size()>=2)?(line.compare(0,2,"# ")):(line.compare(0,1,"#"));
        if (!loop_)
        {
            if (line.compare(0,5,"loop_")) continue;
            while(1)
            {
                if (compress_type==-1)
                {
                    if (cin.good()) getline(cin, line);
                    else PrintErrorAndQuit("ERROR! Unexpected end of "+filename);
                }
                else if (compress_type)
                {
                    if (fin_gz.good()) getline(fin_gz, line);
                    else PrintErrorAndQuit("ERROR! Unexpected end of "+filename);
                }
                else
                {
                    if (fin.good()) getline(fin, line);
                    else PrintErrorAndQuit("ERROR! Unexpected end of "+filename);
                }
                if (line.size()) break;
            }
            if (line.compare(0,11,"_atom_site.")) continue; 
           
            loop_=true;
            _atom_site.clear();
            atom_site_pos=0;
            _atom_site[Trim(line.substr(11))]=atom_site_pos;

            while(1)
            {
                if  (compress_type==-1) getline(cin, line);
                else if (compress_type) getline(fin_gz, line);
                else                    getline(fin, line);
                if (line.size()==0) continue;
                if (line.compare(0,11,"_atom_site.")) break;
                _atom_site[Trim(line.substr(11))]=++atom_site_pos;
            }

            if (_atom_site.count("group_PDB")*
                _atom_site.count("label_atom_id")*
                _atom_site.count("label_comp_id")*
               (_atom_site.count("auth_asym_id")+
                _atom_site.count("label_asym_id"))*
               (_atom_site.count("auth_seq_id")+
                _atom_site.count("label_seq_id"))*
                _atom_site.count("Cartn_x")*
                _atom_site.count("Cartn_y")*
                _atom_site.count("Cartn_z")==0)
            {
                loop_ = false;
                cerr<<"Warning! Missing one of the following _atom_site data items: group_PDB, label_atom_id, label_comp_id, auth_asym_id/label_asym_id, auth_seq_id/label_seq_id, Cartn_x, Cartn_y, Cartn_z"<<endl;
                continue;
            }
        }

        line_vec.clear();
        split(line,line_vec);
        atom     =line_vec[_atom_site["label_atom_id"]];
        resn     =line_vec[_atom_site["label_comp_id"]];
        group_PDB=line_vec[_atom_site["group_PDB"]];
        if (group_PDB=="ATOM") group_PDB="ATOM  ";
        if (mse_opt && resn=="MSE")
        {
            group_PDB="ATOM  ";
            if (atom=="SE") atom="SD";
        }
        if (group_PDB=="HETATM")
        {
            if (!lig_opt) continue; 
            if (!hoh_opt && resn=="HOH") continue;
            if (asym_id!=prev_asym_id && prev_asym_id.size()) 
                asym_id=prev_asym_id; // no separate chain for ligand
        }
        else if (group_PDB!="ATOM  ") continue;
            
        //alt_id=".";
        //if (_atom_site.count("label_alt_id")) // in 39.4 % of entries
            //alt_id=line_vec[_atom_site["label_alt_id"]];
        //if (alt_id!="." && alt_id!="A") continue;

        if (resn.size()==1)
        {
            if (!rna_opt && group_PDB=="ATOM  ") continue;
            resn="  "+resn;
        }
        else if (resn.size()==2)
        {
            if (!dna_opt && resn[0]=='D' && group_PDB=="ATOM  ") continue;
            resn=" " +resn;
        }
        else if (resn.size()==3 && !protein_opt && group_PDB=="ATOM  ") continue;
        else if (resn.size()>=4) resn=resn.substr(0,3);

        if (atom[0]=='"') atom=atom.substr(1);
        if (atom.size() && atom[atom.size()-1]=='"')
            atom=atom.substr(0,atom.size()-1);
        if      (atom.size()==0) continue;
        else if (atom.size()==1) atom=" "+atom+"  ";
        else if (atom.size()==2) atom=" "+atom+" "; // wrong for sidechain H
        else if (atom.size()==3) atom=" "+atom;
        else if (atom.size()>=5) atom=atom.substr(0,4);
        
        if (_atom_site.count("auth_seq_id"))
             resi=line_vec[_atom_site["auth_seq_id"]];
        else resi=line_vec[_atom_site["label_seq_id"]];
        if (_atom_site.count("pdbx_PDB_ins_code") && 
            line_vec[_atom_site["pdbx_PDB_ins_code"]]!="?")
            resi+=line_vec[_atom_site["pdbx_PDB_ins_code"]][0];
        else resi+=" ";
        if (resi.size()>5)
        {
            cerr<<"WARNING! Cannot parse line due to long residue index\n"<<line<<endl;
            continue;
        }

        if (_atom_site.count("auth_asym_id"))
             asym_id=line_vec[_atom_site["auth_asym_id"]];
        else asym_id=line_vec[_atom_site["label_asym_id"]];
        if (asym_id==".") asym_id=" ";
        if (chain_opt.size() && asym_id!=chain_opt &&
            !(asym_id==" " && (chain_opt=="_" || chain_opt=="."))) continue;
            
        if (_atom_site.count("pdbx_PDB_model_num") && 
            model_index!=line_vec[_atom_site["pdbx_PDB_model_num"]])
        {
            if (PDB_lines.size()) break;
            model_index=line_vec[_atom_site["pdbx_PDB_model_num"]];
            map<string, string>().swap(alt_id_dict);
        }

        if (_atom_site.count("label_alt_id")) // in 39.4 % of entries
        {
            alt_id=line_vec[_atom_site["label_alt_id"]];
            if (alt_id==".") alt_id=" ";
            else
            {
                resi_chain=asym_id+resi;
                if (alt_id_dict.count(resi_chain)==0)
                    alt_id_dict[resi_chain]=alt_id;
                else if (alt_id_dict.count(resi_chain) &&
                       alt_id!=alt_id_dict[resi_chain]) continue;
                //if (alt_id!="." && alt_id!="A") continue;
            }
        }

        if (prev_asym_id!=asym_id)
        {
            PDB_lines.push_back(tmp_str_vec);
            chainID_list.push_back(asym_id);
            prev_asym_id=asym_id;
        }

        a++;
        a%=100000;
        i8_stream<<group_PDB
            <<setw(5)<<a<<" "<<atom<<alt_id<<resn<<" "<<asym_id[asym_id.size()-1]
            <<setw(5)<<resi<<"   "
            <<setw(8)<<line_vec[_atom_site["Cartn_x"]].substr(0,8)
            <<setw(8)<<line_vec[_atom_site["Cartn_y"]].substr(0,8)
            <<setw(8)<<line_vec[_atom_site["Cartn_z"]].substr(0,8);
        if (_atom_site.count("B_iso_or_equiv"))
        {
            i8_stream<<"  1.00"<<setw(6)<<line_vec[_atom_site["B_iso_or_equiv"]].substr(0,6);
            if (_atom_site.count("type_symbol"))
                i8_stream<<setw(12)<<line_vec[_atom_site["type_symbol"]].substr(0,12);
        }
        i8_stream<<endl;
        PDB_lines.back().push_back(i8_stream.str());
        i8_stream.str(string());
    }
    group_PDB.clear();
    _atom_site.clear();
    line_vec.clear();
    alt_id.clear();
    asym_id.clear();
    resn.clear();
    map<string, string>().swap(alt_id_dict);
    resi_chain.clear();

    if (compress_type>=0)
    {
        if (compress_type) fin_gz.close();
        else               fin.close();
    }
    line.clear();
    chainID_list.push_back("");
    return PDB_lines.size();
}


int main(int argc, char *argv[])
{
    if (argc < 2) print_help();

    /**********************/
    /*    get argument    */
    /**********************/
    string xname       = "";
    string yname       = "";

    int    split_opt =0;     // do not split chain
    int    het_opt   =0;     // do not read HETATM residues
    int    mol_opt   =7;     // auto-detect the molecule type as protein/RNA
    string chain_opt ="";    // read all chains

    for(int i = 1; i < argc; i++)
    {
        if ( !strcmp(argv[i],"-split") && i < (argc-1) )
        {
            split_opt=atoi(argv[i + 1]); i++;
        }
        else if ( !strcmp(argv[i],"-mol") && i < (argc-1) )
        {
            mol_opt=atoi(argv[i + 1]); i++;
        }
        else if ( !strcmp(argv[i],"-chain") && i < (argc-1) )
        {
            chain_opt=argv[i + 1]; i++;
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
        if (xname.size()==0)
            PrintErrorAndQuit("Please provide input structures");
        else if (yname.size()==0) yname="-";
    }

    bool dna_opt=(mol_opt>=4);
    mol_opt %= 4;
    bool rna_opt=(mol_opt>=2);
    mol_opt %= 2;
    bool protein_opt=(mol_opt>=1);

    if (split_opt<0 || split_opt>1)
        PrintErrorAndQuit("-split can only be 0 or 1");
    if (het_opt<0 || het_opt>3)
        PrintErrorAndQuit("-het can only be 0, 1, 2, or 3");

    bool hoh_opt=(het_opt==3);
    bool lig_opt=(het_opt>=2);
    bool mse_opt=(het_opt>=1);

    /* parse structure */
    vector<vector<string> >PDB_lines;
    vector<string> chainID_list;
    get_all_mmcif_lines(xname, chain_opt, PDB_lines, chainID_list,
        dna_opt, rna_opt, protein_opt, hoh_opt, lig_opt, mse_opt);
    if (!split_opt) resolve_chainID_for_mmcif(PDB_lines,chainID_list);
    write_mmcif_to_pdb(yname, PDB_lines, chainID_list, split_opt);
    
    /* clean up */
    vector<vector<string> >().swap(PDB_lines);
    vector<string>().swap(chainID_list);
    chain_opt.clear();
    return 0;
}
