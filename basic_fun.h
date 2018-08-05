/*
=============================================================
   Implementation of TM-align in C/C++   

   This program is written by Jianyi Yang at
   Yang Zhang lab
   And it is updated by Jianjie Wu at
   Yang Zhang lab
   Department of Computational Medicine and Bioinformatics 
   University of Michigan 
   100 Washtenaw Avenue, Ann Arbor, MI 48109-2218 
           
   Please report bugs and questions to zhng@umich.edu
=============================================================
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <malloc.h>

#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <iterator>
#include <algorithm>
#include <string>
#include <iomanip>
#include <map>

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
};

template <class A> void DeleteArray(A *** array, int Narray)
{
    for(int i=0; i<Narray; i++)
        if(*(*array+i)) delete [] *(*array+i);
    if(Narray) delete [] (*array);
    (*array)=NULL;
};


char AAmap(string AA)
{
    char A='X';
    if      (AA.compare("ALA")==0) A='A';
    else if (AA.compare("ASX")==0) A='B';
    else if (AA.compare("CYS")==0) A='C';
    else if (AA.compare("ASP")==0) A='D';
    else if (AA.compare("GLU")==0) A='E';
    else if (AA.compare("PHE")==0) A='F';
    else if (AA.compare("GLY")==0) A='G';
    else if (AA.compare("HIS")==0) A='H';
    else if (AA.compare("ILE")==0) A='I';
    else if (AA.compare("LYS")==0) A='K';
    else if (AA.compare("LEU")==0) A='L';
    else if (AA.compare("MET")==0 || AA.compare("MSE")==0) A='M';
    else if (AA.compare("ASN")==0) A='N';
    else if (AA.compare("PYL")==0) A='O';
    else if (AA.compare("PRO")==0) A='P';
    else if (AA.compare("GLN")==0) A='Q';
    else if (AA.compare("ARG")==0) A='R';
    else if (AA.compare("SER")==0) A='S';
    else if (AA.compare("THR")==0) A='T';
    else if (AA.compare("SEC")==0) A='U';
    else if (AA.compare("VAL")==0) A='V';
    else if (AA.compare("TRP")==0) A='W';    
    else if (AA.compare("TYR")==0) A='Y';
    else if (AA.compare("GLX")==0) A='Z';

    else if (AA.compare(0,2," D")==0) A=tolower(AA[2]);
    else if (AA.compare(0,2,"  ")==0) A=tolower(AA[2]);
    return A;
}

void get_xyz(string line, double *x, double *y, double *z, char *resname, int *no)
{
    char cstr[50];    
    
    strcpy(cstr, (line.substr(30, 8)).c_str());
    sscanf(cstr, "%lf", x);
    
    strcpy(cstr, (line.substr(38, 8)).c_str());
    sscanf(cstr, "%lf", y);
    
    strcpy(cstr, (line.substr(46, 8)).c_str());
    sscanf(cstr, "%lf", z);
    
    strcpy(cstr, (line.substr(17, 3)).c_str());
    *resname = AAmap(cstr);

    strcpy(cstr, (line.substr(22, 4)).c_str());
    sscanf(cstr, "%d", no);
}

int get_PDB_lines(const char *filename, vector<string> &PDB_lines, 
    vector<string> &resi_vec, const int byresi_opt, 
    const int ter_opt=3, const string atom_opt="auto")
{
    int i=0; // resi
    string line, str, i8;    
    char chainID=0;
    string resi="";
    bool select_atom=false;
    
    ifstream fin (filename);
    if (fin.is_open())
    {
        while (fin.good())
        {
            getline(fin, line);
            if (i > 0)
            {
                if      (ter_opt>=1 && line.compare(0,3,"END")==0) break;
                else if (ter_opt>=3 && line.compare(0,3,"TER")==0) break;
            }
            if (line.compare(0, 6, "ATOM  ")==0 && line.size()>=54 &&
               (line[16]==' ' || line[16]=='A'))
            {
                if (atom_opt=="auto")
                {
                    if (line[17]==' ' && (line[18]=='D'||line[18]==' '))
                        select_atom=(line.compare(12,4," C3'")==0);
                    else
                        select_atom=(line.compare(12,4," CA ")==0);
                }
                else    select_atom=(line.compare(12,4,atom_opt)==0);
                if (select_atom)
                {
                    if (!chainID) chainID=line[21];
                    else if (ter_opt>=2 && chainID!=line[21]) break;

                    if (resi==line.substr(22,5))
                        cerr<<"Warning! Duplicated residue "<<resi<<endl;
                    resi=line.substr(22,5); // including insertion code
                    if (byresi_opt==1) resi_vec.push_back(resi);
                    if (byresi_opt>=2) resi_vec.push_back(resi+line[21]);

                    /* change residue index in line */
                    stringstream i8_stream;
                    i8_stream << i;
                    i8=i8_stream.str();
                    if (i8.size()<4) i8=string(4-i8.size(), ' ')+i8;
                    line=line.substr(0,22)+i8+line.substr(26);

                    PDB_lines.push_back(line);
                    i++;
                }
            }
        }
        fin.close();
    }
    else return 0;
    line.clear();
    return i;
}

/* extract pairwise sequence alignment from residue index vectors,
 * assuming that "sequence" contains two empty strings.
 * return length of alignment, including gap. */
int extract_aln_from_resi(vector<string> &sequence, char *seqx, char *seqy,
    const vector<string> resi_vec1, const vector<string> resi_vec2,
    const int byresi_opt)
{
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

int read_PDB(const vector<string> &PDB_lines, double **a, char *seq, int *resno)
{
    int i;
    for (i=0;i<PDB_lines.size();i++)
        get_xyz(PDB_lines[i], &a[i][0], &a[i][1], &a[i][2], &seq[i], &resno[i]);
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
string Trim(string inputString)
{
    string result = inputString;
    int idxBegin = inputString.find_first_not_of(" \n\r\t");
    int idxEnd = inputString.find_last_not_of(" \n\r\t");
    if (idxBegin >= 0 && idxEnd >= 0)
        result = inputString.substr(idxBegin, idxEnd + 1 - idxBegin);
    return result;
}
