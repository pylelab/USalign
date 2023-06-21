#include <fstream>
#include <map>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include "pstream.h"

using namespace std;

void print_help()
{
    cout <<
"Add chain ID to PDB format file.\n"
"\n"
"Usage: addChainID input.pdb output.pdb chainID\n"
    <<endl;
    exit(EXIT_SUCCESS);
}

void splitlines(const string &line, vector<string> &lines,
    const char delimiter='\n')
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
            lines.push_back("");
        }
        lines.back()+=line[pos];
    }
}

void addChainID(const string &infile,const string &outfile,
    const string &chainID)
{
    stringstream buf;
    if (infile=="-") buf<<cin.rdbuf();
#if defined(REDI_PSTREAM_H_SEEN)
    else if (infile.size()>3 && infile.substr(infile.size()-3)==".gz")
    {
        redi::ipstream fp_gz; // if file is compressed
        fp_gz.open("gunzip -c "+infile);
        buf<<fp_gz.rdbuf();
        fp_gz.close();
    }
#endif
    else
    {
        ifstream fp;
        fp.open(infile.c_str(),ios::in); //ifstream fp(filename,ios::in);
        buf<<fp.rdbuf();
        fp.close();
    }
    vector<string> lines;
    splitlines(buf.str(),lines);
    buf.str(string());
    size_t l;

    for (l=0;l<lines.size();l++)
    {
        if (lines[l].substr(0,6)=="ATOM  " ||
            lines[l].substr(0,6)=="HETATM")
        {
            if (lines[l].size()<22)
            {
                cerr<<"incomplete:"<<lines[l]<<endl;
                continue;
            }
            buf<<lines[l].substr(0,20)<<chainID<<lines[l].substr(22)<<endl;
        }
        else if (lines[l].size())
            buf<<lines[l]<<endl;
        lines[l].clear();
    }

    if (outfile=="-")
        cout<<buf.str();
    else
    {
        ofstream fout;
        fout.open(outfile.c_str(),ios::out);
        fout<<buf.str();
        fout.close();
    }
    buf.str(string());
    vector<string>().swap(lines);
    return;
}

int main(int argc, char *argv[])
{
    if (argc < 2) print_help();

    string infile ="";
    string outfile="";
    string chainID="";

    for (int i=1; i<argc; i++)
    {
        if      ( infile.size()==0) infile =argv[i];
        else if (outfile.size()==0) outfile=argv[i];
        else if (chainID.size()==0) chainID=argv[i];
        else
        {
            cerr<<"ERROR! no such option "<<argv[i]<<endl;
            exit(1);
        }
    }

    if (outfile.size()==0) outfile="-";
    if (chainID.size()==0) chainID="  ";
    else if (chainID.size()==1) chainID=" "+chainID;
    else if (chainID.size()>2) chainID=chainID.substr(chainID.size()-2,2);

    addChainID(infile,outfile,chainID);
    
    infile.clear();
    outfile.clear();
    return 0;
}
