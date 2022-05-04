/* File parsing and basic geometry operations */
#ifndef TMalign_basic_fun_h
#define TMalign_basic_fun_h 1

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
//#include <malloc.h>

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


void PrintErrorAndQuit(const string sErrorString)
{
    cout << sErrorString << endl;
    exit(1);
}

template <typename T> inline T getmin(const T &a, const T &b)
{
    return b<a?b:a;
}

template <class A> void NewArray(A *** array, int Narray1, int Narray2)
{
    *array=new A* [Narray1];
    for(int i=0; i<Narray1; i++) *(*array+i)=new A [Narray2];
}

template <class A> void DeleteArray(A *** array, int Narray)
{
    for(int i=0; i<Narray; i++)
        if(*(*array+i)) delete [] *(*array+i);
    if(Narray) delete [] (*array);
    (*array)=NULL;
}

string AAmap(char A)
{
    if (A=='A') return "ALA";
    if (A=='B') return "ASX";
    if (A=='C') return "CYS";
    if (A=='D') return "ASP";
    if (A=='E') return "GLU";
    if (A=='F') return "PHE";
    if (A=='G') return "GLY";
    if (A=='H') return "HIS";
    if (A=='I') return "ILE";
    if (A=='K') return "LYS";
    if (A=='L') return "LEU";
    if (A=='M') return "MET";
    if (A=='N') return "ASN";
    if (A=='O') return "PYL";
    if (A=='P') return "PRO";
    if (A=='Q') return "GLN";
    if (A=='R') return "ARG";
    if (A=='S') return "SER";
    if (A=='T') return "THR";
    if (A=='U') return "SEC";
    if (A=='V') return "VAL";
    if (A=='W') return "TRP";    
    if (A=='Y') return "TYR";
    if (A=='Z') return "GLX";
    if ('a'<=A && A<='z') return "  "+string(1,char(toupper(A)));
    return "UNK";
}

char AAmap(const string &AA)
{
    if (AA.compare("ALA")==0 || AA.compare("DAL")==0) return 'A';
    if (AA.compare("ASX")==0) return 'B';
    if (AA.compare("CYS")==0 || AA.compare("DCY")==0) return 'C';
    if (AA.compare("ASP")==0 || AA.compare("DAS")==0) return 'D';
    if (AA.compare("GLU")==0 || AA.compare("DGL")==0) return 'E';
    if (AA.compare("PHE")==0 || AA.compare("DPN")==0) return 'F';
    if (AA.compare("GLY")==0) return 'G';
    if (AA.compare("HIS")==0 || AA.compare("DHI")==0) return 'H';
    if (AA.compare("ILE")==0 || AA.compare("DIL")==0) return 'I';
    if (AA.compare("LYS")==0 || AA.compare("DLY")==0) return 'K';
    if (AA.compare("LEU")==0 || AA.compare("DLE")==0) return 'L';
    if (AA.compare("MET")==0 || AA.compare("MED")==0 ||
        AA.compare("MSE")==0) return 'M';
    if (AA.compare("ASN")==0 || AA.compare("DSG")==0) return 'N';
    if (AA.compare("PYL")==0) return 'O';
    if (AA.compare("PRO")==0 || AA.compare("DPR")==0) return 'P';
    if (AA.compare("GLN")==0 || AA.compare("DGN")==0) return 'Q';
    if (AA.compare("ARG")==0 || AA.compare("DAR")==0) return 'R';
    if (AA.compare("SER")==0 || AA.compare("DSN")==0) return 'S';
    if (AA.compare("THR")==0 || AA.compare("DTH")==0) return 'T';
    if (AA.compare("SEC")==0) return 'U';
    if (AA.compare("VAL")==0 || AA.compare("DVA")==0) return 'V';
    if (AA.compare("TRP")==0 || AA.compare("DTR")==0) return 'W';    
    if (AA.compare("TYR")==0 || AA.compare("DTY")==0) return 'Y';
    if (AA.compare("GLX")==0) return 'Z';

    if (AA.compare(0,2," D")==0) return tolower(AA[2]);
    if (AA.compare(0,2,"  ")==0) return tolower(AA[2]);
    return 'X';
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

size_t get_PDB_lines(const string filename,
    vector<vector<string> >&PDB_lines, vector<string> &chainID_list,
    vector<int> &mol_vec, const int ter_opt, const int infmt_opt,
    const string atom_opt, const int split_opt, const int het_opt)
{
    size_t i=0; // resi i.e. atom index
    string line;
    char chainID=0;
    string resi="";
    bool select_atom=false;
    size_t model_idx=0;
    vector<string> tmp_str_vec;
    
    int compress_type=0; // uncompressed file
    ifstream fin;
#ifndef REDI_PSTREAM_H_SEEN
    ifstream fin_gz;
#else
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
#endif
    {
        if (filename=="-") compress_type=-1;
        else fin.open(filename.c_str());
    }

    if (infmt_opt==0||infmt_opt==-1) // PDB format
    {
        while ((compress_type==-1)?cin.good():(compress_type?fin_gz.good():fin.good()))
        {
            if  (compress_type==-1) getline(cin, line);
            else if (compress_type) getline(fin_gz, line);
            else                    getline(fin, line);
            if (infmt_opt==-1 && line.compare(0,5,"loop_")==0) // PDBx/mmCIF
                return get_PDB_lines(filename,PDB_lines,chainID_list,
                    mol_vec, ter_opt, 3, atom_opt, split_opt,het_opt);
            if (i > 0)
            {
                if      (ter_opt>=1 && line.compare(0,3,"END")==0) break;
                else if (ter_opt>=3 && line.compare(0,3,"TER")==0) break;
            }
            if (split_opt && line.compare(0,3,"END")==0) chainID=0;
            if (line.size()>=54 && (line[16]==' ' || line[16]=='A') && (
                (line.compare(0, 6, "ATOM  ")==0) || 
                (line.compare(0, 6, "HETATM")==0 && het_opt==1) ||
                (line.compare(0, 6, "HETATM")==0 && het_opt==2 && 
                 line.compare(17,3, "MSE")==0)))
            {
                if (atom_opt=="auto")
                {
                    if (line[17]==' ' && (line[18]=='D'||line[18]==' '))
                         select_atom=(line.compare(12,4," C3'")==0);
                    else select_atom=(line.compare(12,4," CA ")==0);
                }
                else if (atom_opt=="PC4'")
                {
                    if (line[17]==' ' && (line[18]=='D'||line[18]==' '))
                         select_atom=(line.compare(12,4," P  ")==0
                                  )||(line.compare(12,4," C4'")==0);
                    else select_atom=(line.compare(12,4," CA ")==0);
                }
                else     select_atom=(line.compare(12,4,atom_opt)==0);
                if (select_atom)
                {
                    if (!chainID)
                    {
                        chainID=line[21];
                        model_idx++;
                        stringstream i8_stream;
                        i=0;
                        if (split_opt==2) // split by chain
                        {
                            if (chainID==' ')
                            {
                                if (ter_opt>=1) i8_stream << ":_";
                                else i8_stream<<':'<<model_idx<<",_";
                            }
                            else
                            {
                                if (ter_opt>=1) i8_stream << ':' << chainID;
                                else i8_stream<<':'<<model_idx<<','<<chainID;
                            }
                            chainID_list.push_back(i8_stream.str());
                        }
                        else if (split_opt==1) // split by model
                        {
                            i8_stream << ':' << model_idx;
                            chainID_list.push_back(i8_stream.str());
                        }
                        PDB_lines.push_back(tmp_str_vec);
                        mol_vec.push_back(0);
                    }
                    else if (ter_opt>=2 && chainID!=line[21]) break;
                    if (split_opt==2 && chainID!=line[21])
                    {
                        chainID=line[21];
                        i=0;
                        stringstream i8_stream;
                        if (chainID==' ')
                        {
                            if (ter_opt>=1) i8_stream << ":_";
                            else i8_stream<<':'<<model_idx<<",_";
                        }
                        else
                        {
                            if (ter_opt>=1) i8_stream << ':' << chainID;
                            else i8_stream<<':'<<model_idx<<','<<chainID;
                        }
                        chainID_list.push_back(i8_stream.str());
                        PDB_lines.push_back(tmp_str_vec);
                        mol_vec.push_back(0);
                    }

                    if (resi==line.substr(22,5) && atom_opt!="PC4'")
                        cerr<<"Warning! Duplicated residue "<<resi<<endl;
                    resi=line.substr(22,5); // including insertion code

                    PDB_lines.back().push_back(line);
                    if (line[17]==' ' && (line[18]=='D'||line[18]==' ')) mol_vec.back()++;
                    else mol_vec.back()--;
                    i++;
                }
            }
        }
    }
    else if (infmt_opt==1) // SPICKER format
    {
        size_t L=0;
        float x,y,z;
        stringstream i8_stream;
        while ((compress_type==-1)?cin.good():(compress_type?fin_gz.good():fin.good()))
        {
            if  (compress_type==-1)
            {
                cin>>L>>x>>y>>z;
                getline(cin, line);
                if (!cin.good()) break;
            }
            else if (compress_type)
            {
                fin_gz>>L>>x>>y>>z;
                getline(fin_gz, line);
                if (!fin_gz.good()) break;
            }
            else
            {
                fin   >>L>>x>>y>>z;
                getline(fin, line);
                if (!fin.good()) break;
            }
            model_idx++;
            stringstream i8_stream;
            i8_stream << ':' << model_idx;
            chainID_list.push_back(i8_stream.str());
            PDB_lines.push_back(tmp_str_vec);
            mol_vec.push_back(0);
            for (i=0;i<L;i++)
            {
                if  (compress_type==-1) cin>>x>>y>>z;
                else if (compress_type) fin_gz>>x>>y>>z;
                else                    fin   >>x>>y>>z;
                i8_stream<<"ATOM   "<<setw(4)<<i+1<<"  CA  UNK  "<<setw(4)
                    <<i+1<<"    "<<setiosflags(ios::fixed)<<setprecision(3)
                    <<setw(8)<<x<<setw(8)<<y<<setw(8)<<z;
                line=i8_stream.str();
                i8_stream.str(string());
                PDB_lines.back().push_back(line);
            }
            if  (compress_type==-1) getline(cin, line);
            else if (compress_type) getline(fin_gz, line);
            else                    getline(fin, line);
        }
    }
    else if (infmt_opt==2) // xyz format
    {
        size_t L=0;
        stringstream i8_stream;
        while ((compress_type==-1)?cin.good():(compress_type?fin_gz.good():fin.good()))
        {
            if (compress_type==-1)  getline(cin, line);
            else if (compress_type) getline(fin_gz, line);
            else                    getline(fin, line);
            L=atoi(line.c_str());
            if (compress_type==-1)  getline(cin, line);
            else if (compress_type) getline(fin_gz, line);
            else                    getline(fin, line);
            for (i=0;i<line.size();i++)
                if (line[i]==' '||line[i]=='\t') break;
            if (!((compress_type==-1)?cin.good():(compress_type?fin_gz.good():fin.good()))) break;
            chainID_list.push_back(':'+line.substr(0,i));
            PDB_lines.push_back(tmp_str_vec);
            mol_vec.push_back(0);
            for (i=0;i<L;i++)
            {
                if (compress_type==-1)  getline(cin, line);
                else if (compress_type) getline(fin_gz, line);
                else                    getline(fin, line);
                i8_stream<<"ATOM   "<<setw(4)<<i+1<<"  CA  "
                    <<AAmap(line[0])<<"  "<<setw(4)<<i+1<<"    "
                    <<line.substr(2,8)<<line.substr(11,8)<<line.substr(20,8);
                line=i8_stream.str();
                i8_stream.str(string());
                PDB_lines.back().push_back(line);
                if (line[0]>='a' && line[0]<='z') mol_vec.back()++; // RNA
                else mol_vec.back()--;
            }
        }
    }
    else if (infmt_opt==3) // PDBx/mmCIF format
    {
        bool loop_ = false; // not reading following content
        map<string,int> _atom_site;
        int atom_site_pos;
        vector<string> line_vec;
        string alt_id=".";  // alternative location indicator
        string asym_id="."; // this is similar to chainID, except that
                            // chainID is char while asym_id is a string
                            // with possibly multiple char
        string prev_asym_id="";
        string AA="";       // residue name
        string atom="";
        string prev_resi="";
        string model_index=""; // the same as model_idx but type is string
        stringstream i8_stream;
        while ((compress_type==-1)?cin.good():(compress_type?fin_gz.good():fin.good()))
        {
            if (compress_type==-1)  getline(cin, line);
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
                        else PrintErrorAndQuit("ERROR! Unexpected end of -");
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
            if ((line_vec[_atom_site["group_PDB"]]!="ATOM" &&
                 line_vec[_atom_site["group_PDB"]]!="HETATM") ||
                (line_vec[_atom_site["group_PDB"]]=="HETATM" &&
                 (het_opt==0 || 
                 (het_opt==2 && line_vec[_atom_site["label_comp_id"]]!="MSE")))
                ) continue;
            
            alt_id=".";
            if (_atom_site.count("label_alt_id")) // in 39.4 % of entries
                alt_id=line_vec[_atom_site["label_alt_id"]];
            if (alt_id!="." && alt_id!="A") continue;

            atom=line_vec[_atom_site["label_atom_id"]];
            if (atom[0]=='"') atom=atom.substr(1);
            if (atom.size() && atom[atom.size()-1]=='"')
                atom=atom.substr(0,atom.size()-1);
            if (atom.size()==0) continue;
            if      (atom.size()==1) atom=" "+atom+"  ";
            else if (atom.size()==2) atom=" "+atom+" "; // wrong for sidechain H
            else if (atom.size()==3) atom=" "+atom;
            else if (atom.size()>=5) continue;

            AA=line_vec[_atom_site["label_comp_id"]]; // residue name
            if      (AA.size()==1) AA="  "+AA;
            else if (AA.size()==2) AA=" " +AA;
            else if (AA.size()>=4) continue;

            if (atom_opt=="auto")
            {
                if (AA[0]==' ' && (AA[1]=='D'||AA[1]==' ')) // DNA || RNA
                     select_atom=(atom==" C3'");
                else select_atom=(atom==" CA ");
            }
            else if (atom_opt=="PC4'")
            {
                if (line[17]==' ' && (line[18]=='D'||line[18]==' '))
                     select_atom=(line.compare(12,4," P  ")==0
                              )||(line.compare(12,4," C4'")==0);
                else select_atom=(line.compare(12,4," CA ")==0);
            }
            else     select_atom=(atom==atom_opt);

            if (!select_atom) continue;

            if (_atom_site.count("auth_asym_id"))
                 asym_id=line_vec[_atom_site["auth_asym_id"]];
            else asym_id=line_vec[_atom_site["label_asym_id"]];
            if (asym_id==".") asym_id=" ";
            
            if (_atom_site.count("pdbx_PDB_model_num") && 
                model_index!=line_vec[_atom_site["pdbx_PDB_model_num"]])
            {
                model_index=line_vec[_atom_site["pdbx_PDB_model_num"]];
                if (PDB_lines.size() && ter_opt>=1) break;
                if (PDB_lines.size()==0 || split_opt>=1)
                {
                    PDB_lines.push_back(tmp_str_vec);
                    mol_vec.push_back(0);
                    prev_asym_id=asym_id;

                    if (split_opt==1 && ter_opt==0) chainID_list.push_back(
                        ':'+model_index);
                    else if (split_opt==2 && ter_opt==0)
                        chainID_list.push_back(':'+model_index+','+asym_id);
                    else //if (split_opt==2 && ter_opt==1)
                        chainID_list.push_back(':'+asym_id);
                    //else
                        //chainID_list.push_back("");
                }
            }

            if (prev_asym_id!=asym_id)
            {
                if (prev_asym_id!="" && ter_opt>=2) break;
                if (split_opt>=2)
                {
                    PDB_lines.push_back(tmp_str_vec);
                    mol_vec.push_back(0);

                    if (split_opt==1 && ter_opt==0) chainID_list.push_back(
                        ':'+model_index);
                    else if (split_opt==2 && ter_opt==0)
                        chainID_list.push_back(':'+model_index+','+asym_id);
                    else //if (split_opt==2 && ter_opt==1)
                        chainID_list.push_back(':'+asym_id);
                    //else
                        //chainID_list.push_back("");
                }
            }
            if (prev_asym_id!=asym_id) prev_asym_id=asym_id;

            if (AA[0]==' ' && (AA[1]=='D'||AA[1]==' ')) mol_vec.back()++;
            else mol_vec.back()--;

            if (_atom_site.count("auth_seq_id"))
                 resi=line_vec[_atom_site["auth_seq_id"]];
            else resi=line_vec[_atom_site["label_seq_id"]];
            if (_atom_site.count("pdbx_PDB_ins_code") && 
                line_vec[_atom_site["pdbx_PDB_ins_code"]]!="?")
                resi+=line_vec[_atom_site["pdbx_PDB_ins_code"]][0];
            else resi+=" ";

            if (prev_resi==resi && atom_opt!="PC4'")
                cerr<<"Warning! Duplicated residue "<<resi<<endl;
            prev_resi=resi;

            i++;
            i8_stream<<"ATOM  "
                <<setw(5)<<i<<" "<<atom<<" "<<AA<<" "<<asym_id[0]
                <<setw(5)<<resi.substr(0,5)<<"   "
                <<setw(8)<<line_vec[_atom_site["Cartn_x"]].substr(0,8)
                <<setw(8)<<line_vec[_atom_site["Cartn_y"]].substr(0,8)
                <<setw(8)<<line_vec[_atom_site["Cartn_z"]].substr(0,8);
            PDB_lines.back().push_back(i8_stream.str());
            i8_stream.str(string());
        }
        _atom_site.clear();
        line_vec.clear();
        alt_id.clear();
        asym_id.clear();
        AA.clear();
    }

    if      (compress_type>=1) fin_gz.close();
    else if (compress_type==0) fin.close();
    line.clear();
    if (!split_opt) chainID_list.push_back("");
    return PDB_lines.size();
}

/* read fasta file from filename. sequence is stored into FASTA_lines
 * while sequence name is stored into chainID_list.
 * if ter_opt >=1, only read the first sequence.
 * if ter_opt ==0, read all sequences.
 * if split_opt >=1 and ter_opt ==0, each sequence is a separate entry.
 * if split_opt ==0 and ter_opt ==0, all sequences are combined into one */
size_t get_FASTA_lines(const string filename,
    vector<vector<string> >&FASTA_lines, vector<string> &chainID_list,
    vector<int> &mol_vec, const int ter_opt=3, const int split_opt=0)
{
    string line;
    vector<string> tmp_str_vec;
    size_t l;
    
    int compress_type=0; // uncompressed file
    ifstream fin;
#ifndef REDI_PSTREAM_H_SEEN
    ifstream fin_gz;
#else
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
#endif
    {
        if (filename=="-") compress_type=-1;
        else fin.open(filename.c_str());
    }

    while ((compress_type==-1)?cin.good():(compress_type?fin_gz.good():fin.good()))
    {
        if  (compress_type==-1) getline(cin, line);
        else if (compress_type) getline(fin_gz, line);
        else                    getline(fin, line);

        if (line.size()==0 || line[0]=='#') continue;

        if (line[0]=='>')
        {
            if (FASTA_lines.size())
            {
                if (ter_opt) break;
                if (split_opt==0) continue;
            }
            FASTA_lines.push_back(tmp_str_vec);
            FASTA_lines.back().push_back("");
            mol_vec.push_back(0);
            if (ter_opt==0 && split_opt)
            {
                line[0]=':';
                chainID_list.push_back(line);
            }
            else chainID_list.push_back("");
        }
        else
        {
            FASTA_lines.back()[0]+=line;
            for (l=0;l<line.size();l++) mol_vec.back()+=
                ('a'<=line[l] && line[l]<='z')-('A'<=line[l] && line[l]<='Z');
        }
    }

    line.clear();
    if      (compress_type>=1) fin_gz.close();
    else if (compress_type==0) fin.close();
    return FASTA_lines.size();
}

int read_PDB(const vector<string> &PDB_lines, double **a, char *seq,
    vector<string> &resi_vec, const int read_resi)
{
    size_t i;
    for (i=0;i<PDB_lines.size();i++)
    {
        a[i][0] = atof(PDB_lines[i].substr(30, 8).c_str());
        a[i][1] = atof(PDB_lines[i].substr(38, 8).c_str());
        a[i][2] = atof(PDB_lines[i].substr(46, 8).c_str());
        seq[i]  = AAmap(PDB_lines[i].substr(17, 3));

        if (read_resi>=2) resi_vec.push_back(PDB_lines[i].substr(22,5)+
                                             PDB_lines[i][21]);
        if (read_resi==1) resi_vec.push_back(PDB_lines[i].substr(22,5));
    }
    seq[i]='\0'; 
    return i;
}

double dist(double x[3], double y[3])
{
    double d1=x[0]-y[0];
    double d2=x[1]-y[1];
    double d3=x[2]-y[2];
 
    return (d1*d1 + d2*d2 + d3*d3);
}

double dot(double *a, double *b)
{
    return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
}

void transform(double t[3], double u[3][3], double *x, double *x1)
{
    x1[0]=t[0]+dot(&u[0][0], x);
    x1[1]=t[1]+dot(&u[1][0], x);
    x1[2]=t[2]+dot(&u[2][0], x);
}

void do_rotation(double **x, double **x1, int len, double t[3], double u[3][3])
{
    for(int i=0; i<len; i++)
    {
        transform(t, u, &x[i][0], &x1[i][0]);
    }    
}

/* read user specified pairwise alignment from 'fname_lign' to 'sequence'.
 * This function should only be called by main function, as it will
 * terminate a program if wrong alignment is given */
void read_user_alignment(vector<string>&sequence, const string &fname_lign,
    const int i_opt)
{
    if (fname_lign == "")
        PrintErrorAndQuit("Please provide a file name for option -i!");
    // open alignment file
    int n_p = 0;// number of structures in alignment file
    string line;
    
    ifstream fileIn(fname_lign.c_str());
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
    else PrintErrorAndQuit("ERROR! Alignment file does not exist.");
    
    if (n_p < 2)
        PrintErrorAndQuit("ERROR: Fasta format is wrong, two proteins should be included.");
    if (sequence[0].size() != sequence[1].size())
        PrintErrorAndQuit("ERROR! FASTA file is wrong. The length in alignment should be equal for the two aligned proteins.");
    if (i_opt==3)
    {
        int aligned_resNum=0;
        for (size_t i=0;i<sequence[0].size();i++)
            aligned_resNum+=(sequence[0][i]!='-' && sequence[1][i]!='-');
        if (aligned_resNum<3)
            PrintErrorAndQuit("ERROR! Superposition is undefined for <3 aligned residues.");
    }
    line.clear();
    return;
}

/* read list of entries from 'name' to 'chain_list'.
 * dir_opt is the folder name (prefix).
 * suffix_opt is the file name extension (suffix_opt).
 * This function should only be called by main function, as it will
 * terminate a program if wrong alignment is given */
void file2chainlist(vector<string>&chain_list, const string &name,
    const string &dir_opt, const string &suffix_opt)
{
    ifstream fp(name.c_str());
    if (! fp.is_open())
        PrintErrorAndQuit(("Can not open file: "+name+'\n').c_str());
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

#endif
