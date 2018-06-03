/*
===============================================================================
   Implementation of TM-align in C/C++   

   This program is written by Jianyi Yang at
   Yang Zhang lab
   And it is updated by Jianjie Wu at
   Yang Zhang lab
   Department of Computational Medicine and Bioinformatics 
   University of Michigan 
   100 Washtenaw Avenue, Ann Arbor, MI 48109-2218 
           
   Please report bugs and questions to jianjiew@umich.edu or zhng@umich.edu
===============================================================================
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

#include <map>

#include "basic_define.h"

using namespace std;


void PrintErrorAndQuit(string sErrorString)
{
	cout << sErrorString << endl;
	exit(1);
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
    else if (AA.compare("CYS")==0 || AA.compare("CYX")==0) A='C';
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

    if      (AA.compare(0,2," D")==0) A=tolower(AA[2]);
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

int get_PDB_len(char *filename)
{
	int i = 0;
	string line;
	string atom("ATOM ");

	ifstream fin(filename);
	if (fin.is_open())
	{
		while (fin.good())
		{
			getline(fin, line);
			if (line.compare(0, atom.length(), atom) == 0)
				i++;
		}
		fin.close();
	}
	else
	{
		char message[5000];
		sprintf(message, "Can not open file: %s\n", filename);
		PrintErrorAndQuit(message);
	}

	return i;
}

int read_PDB(char *filename, double **a, char *seq, int *resno, 
    const int ter_opt=3, const string atom_opt=" CA ")
{
    int i=0; // resi
    string line, str, i8;    
    char chainID=0;
	string resn="";
    
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
				if (line.compare(12, 4, atom_opt)==0)
				{
                    if (!chainID) chainID=line[21];
                    else if (ter_opt>=2 && chainID!=line[21]) break;

                    if (resn==line.substr(22,5))
                        cerr<<"Warning! Duplicated residue "<<resn<<endl;
                    resn=line.substr(22,5);

                    // change residue index in line
                    stringstream i8_stream;
                    i8_stream << i;
                    i8=i8_stream.str();
                    if (i8.size()<4)
                    {
                        i8=string(4-i8.size(), ' ')+i8;
                    }
                    line=line.substr(0,22)+i8+line.substr(26);

					get_xyz(line, &a[i][0], &a[i][1], &a[i][2], &seq[i], &resno[i]);
                    i++;
				}
			}
        }
        fin.close();
    }
    else
    {
		char message[5000];
		sprintf(message, "Can not open file: %s\n", filename);
        PrintErrorAndQuit(message);
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

void output_align1(int *invmap0, int len)
{
	for(int i=0; i<len; i++)
	{
		if(invmap0[i]>=0)
		{
			cout << invmap0[i]+1 << " ";
		}			
		else
			cout << invmap0[i] << " ";

	}	
	cout << endl << endl;	
}


int output_align(int *invmap0, int len)
{
	int n_ali=0;
	for(int i=0; i<len; i++)
	{		
		cout <<  invmap0[i] << " ";
		n_ali++;		
	}	
	cout << endl << endl;	

	return n_ali;
}
// output aligned residues to string
string output_align_to_string(int *invmap, int len)
{
	string result = "";
	int step = 10;
	int rest = len % step;
	int loop = len / step;
	int i,k;
	char temp[1000];
	for ( i = 0; i<loop; i++)
	{
		k = i * step;
		for (int j = 0; j < step; j++)
		{
			sprintf(temp, "%4d ", invmap[k + j]);
			result += string(temp);
		}
		result += "\n";
	}
	k = i * step;
	for (int j = 0; j < rest; j++)
	{
		sprintf(temp, "%4d ", invmap[k + j]);
		result += string(temp);
	}
	result += "\n";

	return result;
}

// output aligned coords to string
string output_coord_to_string(double **x, double **y, int len)
{
	string result = "";
	int i, k;
	char temp[1000];
	result += "xtm coords:\n";
	for (i = 0; i<len; i++)
	{
		sprintf(temp, "%8.3f %8.3f %8.3f\n", x[i][0], x[i][1], x[i][2]);
		result += string(temp);
	}
	result += "ytm coords:\n";
	for (i = 0; i<len; i++)
	{
		sprintf(temp, "%8.3f %8.3f %8.3f\n", y[i][0], y[i][1], y[i][2]);
		result += string(temp);
	}
	return result;
}
