/* File parsing and basic geometry operations */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <malloc.h>

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


void PrintErrorAndQuit(string sErrorString)
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
    if ('a'<=A && A<='z') return "  "+toupper(A);
    return "UNK";
}

char AAmap(const string &AA)
{
    if (AA.compare("ALA")==0) return 'A';
    if (AA.compare("ASX")==0) return 'B';
    if (AA.compare("CYS")==0) return 'C';
    if (AA.compare("ASP")==0) return 'D';
    if (AA.compare("GLU")==0) return 'E';
    if (AA.compare("PHE")==0) return 'F';
    if (AA.compare("GLY")==0) return 'G';
    if (AA.compare("HIS")==0) return 'H';
    if (AA.compare("ILE")==0) return 'I';
    if (AA.compare("LYS")==0) return 'K';
    if (AA.compare("LEU")==0) return 'L';
    if (AA.compare("MET")==0 || AA.compare("MSE")==0) return 'M';
    if (AA.compare("ASN")==0) return 'N';
    if (AA.compare("PYL")==0) return 'O';
    if (AA.compare("PRO")==0) return 'P';
    if (AA.compare("GLN")==0) return 'Q';
    if (AA.compare("ARG")==0) return 'R';
    if (AA.compare("SER")==0) return 'S';
    if (AA.compare("THR")==0) return 'T';
    if (AA.compare("SEC")==0) return 'U';
    if (AA.compare("VAL")==0) return 'V';
    if (AA.compare("TRP")==0) return 'W';    
    if (AA.compare("TYR")==0) return 'Y';
    if (AA.compare("GLX")==0) return 'Z';

    if (AA.compare(0,2," D")==0) return tolower(AA[2]);
    if (AA.compare(0,2,"  ")==0) return tolower(AA[2]);
    return 'X';
}

int get_PDB_lines(const string filename,
    vector<vector<string> >&PDB_lines, vector<string> &chainID_list,
    vector<int> &mol_vec, const int ter_opt=3, const int infmt_opt=0,
    const string atom_opt="auto", const int split_opt=0)
{
    int i=0; // resi
    string line;
    char chainID=0;
    string resi="";
    bool select_atom=false;
    int model_idx=0;
    vector<string> tmp_str_vec;
    
    int compress_type=0; // uncompressed file
    ifstream fin;
    redi::ipstream fin_gz; // if file is compressed
    if (filename.size()>=3 && 
        filename.substr(filename.size()-3,3)==".gz")
    {
        fin_gz.open("zcat "+filename);
        compress_type=1;
    }
    else if (filename.size()>=4 && 
        filename.substr(filename.size()-4,4)==".bz2")
    {
        fin_gz.open("bzcat "+filename);
        compress_type=2;
    }
    else fin.open(filename.c_str());

    if (infmt_opt==0) // PDB format
    {
        while (compress_type?fin_gz.good():fin.good())
        {
            if (compress_type) getline(fin_gz, line);
            else               getline(fin, line);
            if (i > 0)
            {
                if      (ter_opt>=1 && line.compare(0,3,"END")==0) break;
                else if (ter_opt>=3 && line.compare(0,3,"TER")==0) break;
            }
            if (split_opt && line.compare(0,3,"END")==0) chainID=0;
            if (line.compare(0, 6, "ATOM  ")==0 && line.size()>=54 &&
               (line[16]==' ' || line[16]=='A'))
            {
                if (atom_opt=="auto")
                {
                    if (line[17]==' ' && (line[18]=='D'||line[18]==' '))
                         select_atom=(line.compare(12,4," C3'")==0);
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
                                else i8_stream<<':'<<model_idx<<":_";
                            }
                            else
                            {
                                if (ter_opt>=1) i8_stream << ':' << chainID;
                                else i8_stream<<':'<<model_idx<<':'<<chainID;
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
                            else i8_stream<<':'<<model_idx<<":_";
                        }
                        else
                        {
                            if (ter_opt>=1) i8_stream << ':' << chainID;
                            else i8_stream<<':'<<model_idx<<':'<<chainID;
                        }
                        chainID_list.push_back(i8_stream.str());
                        PDB_lines.push_back(tmp_str_vec);
                        mol_vec.push_back(0);
                    }

                    if (resi==line.substr(22,5))
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
        int L=0;
        float x,y,z;
        stringstream i8_stream;
        while (compress_type?fin_gz.good():fin.good())
        {
            if (compress_type) fin_gz>>L>>x>>y>>z;
            else               fin   >>L>>x>>y>>z;
            if (compress_type) getline(fin_gz, line);
            else               getline(fin, line);
            if (!(compress_type?fin_gz.good():fin.good())) break;
            model_idx++;
            stringstream i8_stream;
            i8_stream << ':' << model_idx;
            chainID_list.push_back(i8_stream.str());
            PDB_lines.push_back(tmp_str_vec);
            mol_vec.push_back(0);
            for (i=0;i<L;i++)
            {
                if (compress_type) fin_gz>>x>>y>>z;
                else               fin   >>x>>y>>z;
                i8_stream<<"ATOM   "<<setw(4)<<i+1<<"  CA  UNK  "<<setw(4)
                    <<i+1<<"    "<<setiosflags(ios::fixed)<<setprecision(3)
                    <<setw(8)<<x<<setw(8)<<y<<setw(8)<<z;
                line=i8_stream.str();
                i8_stream.str(string());
                PDB_lines.back().push_back(line);
            }
            if (compress_type) getline(fin_gz, line);
            else               getline(fin, line);
        }
    }
    else if (infmt_opt==2) // xyz format
    {
        int L=0;
        char A;
        stringstream i8_stream;
        while (compress_type?fin_gz.good():fin.good())
        {
            if (compress_type) getline(fin_gz, line);
            else               getline(fin, line);
            L=atoi(line.c_str());
            if (compress_type) getline(fin_gz, line);
            else               getline(fin, line);
            for (i=0;i<line.size();i++)
                if (line[i]==' '||line[i]=='\t') break;
            if (!(compress_type?fin_gz.good():fin.good())) break;
            chainID_list.push_back(':'+line.substr(0,i));
            PDB_lines.push_back(tmp_str_vec);
            mol_vec.push_back(0);
            for (i=0;i<L;i++)
            {
                if (compress_type) getline(fin_gz, line);
                else               getline(fin, line);
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
    if (compress_type) fin_gz.close();
    else               fin.close();
    line.clear();
    if (!split_opt) chainID_list.push_back("");
    return PDB_lines.size();
}

/* extract pairwise sequence alignment from residue index vectors,
 * assuming that "sequence" contains two empty strings.
 * return length of alignment, including gap. */
int extract_aln_from_resi(vector<string> &sequence, char *seqx, char *seqy,
    const vector<string> resi_vec1, const vector<string> resi_vec2,
    const int byresi_opt)
{
    sequence.clear();
    sequence.push_back("");
    sequence.push_back("");

    int i1=0; // positions in resi_vec1
    int i2=0; // positions in resi_vec2
    int xlen=resi_vec1.size();
    int ylen=resi_vec2.size();
    map<char,int> chainID_map1;
    map<char,int> chainID_map2;
    if (byresi_opt==3)
    {
        vector<char> chainID_vec;
        char chainID;
        int i;
        for (i=0;i<xlen;i++)
        {
            chainID=resi_vec1[i][5];
            if (!chainID_vec.size()|| chainID_vec.back()!=chainID)
            {
                chainID_vec.push_back(chainID);
                chainID_map1[chainID]=chainID_vec.size();
            }
        }
        chainID_vec.clear();
        for (i=0;i<ylen;i++)
        {
            chainID=resi_vec2[i][5];
            if (!chainID_vec.size()|| chainID_vec.back()!=chainID)
            {
                chainID_vec.push_back(chainID);
                chainID_map2[chainID]=chainID_vec.size();
            }
        }
        chainID_vec.clear();
    }
    while(i1<xlen && i2<ylen)
    {
        if ((byresi_opt<=2 && resi_vec1[i1]==resi_vec2[i2]) || (byresi_opt==3
             && resi_vec1[i1].substr(0,5)==resi_vec2[i2].substr(0,5)
             && chainID_map1[resi_vec1[i1][5]]==chainID_map2[resi_vec2[i2][5]]))
        {
            sequence[0]+=seqx[i1++];
            sequence[1]+=seqy[i2++];
        }
        else if (atoi(resi_vec1[i1].substr(0,4).c_str())<=
                 atoi(resi_vec2[i2].substr(0,4).c_str()))
        {
            sequence[0]+=seqx[i1++];
            sequence[1]+='-';
        }
        else
        {
            sequence[0]+='-';
            sequence[1]+=seqy[i2++];
        }
    }
    chainID_map1.clear();
    chainID_map2.clear();
    return sequence[0].size();
}

int read_PDB(const vector<string> &PDB_lines, double **a, char *seq,
    vector<string> &resi_vec, const int byresi_opt)
{
    int i;
    for (i=0;i<PDB_lines.size();i++)
    {
        a[i][0] = atof(PDB_lines[i].substr(30, 8).c_str());
        a[i][1] = atof(PDB_lines[i].substr(38, 8).c_str());
        a[i][2] = atof(PDB_lines[i].substr(46, 8).c_str());
        seq[i]  = AAmap(PDB_lines[i].substr(17, 3));

        if (byresi_opt>=2) resi_vec.push_back(PDB_lines[i].substr(22,5)+
                                              PDB_lines[i][21]);
        if (byresi_opt==1) resi_vec.push_back(PDB_lines[i].substr(22,5));
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

/* read user specified pairwise alignment from 'fname_lign' to 'sequence'.
 * This function should only be called by main function, as it will
 * terminate a program if wrong alignment is given */
void read_user_alignment(vector<string>&sequence, const string &fname_lign,
    const bool I_opt)
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
        PrintErrorAndQuit("ERROR! FASTA file is wrong. The length in alignment should be equal respectively to the two aligned proteins.");
    if (I_opt)
    {
        int aligned_resNum=0;
        for (int i=0;i<sequence[0].size();i++) 
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
