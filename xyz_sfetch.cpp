#include <iostream>
#include <fstream>
#include "pstream.h"

using namespace std;

void print_help()
{
    cout <<
"Usage: xyz_sfetch ca.xyz\n"
"    List all entries in xyz file 'ca.xyz'\n"
"    them to 'subset.xyz'\n"
"\n"
"Usage: xyz_sfetch ca.xyz list > subset.xyz\n"
"    From xyz file 'ca.xyz', extract all entries listed by 'list'. Output\n"
"    them to 'subset.xyz'\n"
    <<endl;
    exit(EXIT_SUCCESS);
}

int main(int argc, char *argv[])
{
    if (argc < 2) print_help();


    /* get argument */
    string filename="";
    string list_opt="";

    for(int i=1; i<argc; i++)
    {
        if (filename.size()==0) filename=argv[i];
        else list_opt=argv[i];
    }

    if(filename.size()==0||filename=="-h") print_help();

    /* Check xyz file compress type */
    int compress_type=0;   // uncompressed file
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

    /* list all entries in xyz file */
    string line;
    int L,i;
    if (list_opt.size()==0)
    {
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
            cout<<line.substr(0,i)<<'\t'<<L<<endl;
            for (i=0;i<L;i++)
            {
                if (compress_type) getline(fin_gz, line);
                else               getline(fin, line);
            }
        }
        if (compress_type) fin_gz.close();
        else               fin.close();
        /* clean up */
        line.clear();
        filename.clear();
        list_opt.clear();
        return 0;
    }

    /* read entry list */
    vector<string> chain_list;
    ifstream fp(list_opt.c_str());
    while (fp.good())
    {
        getline(fp, line);
        for (i=0;i<line.size();i++)
            if (line[i]==' '||line[i]=='\t') break;
        if (line.size() && i) chain_list.push_back(line.substr(0,i));
    }
    fp.close();

    /* extract listed entries in xyz file */
    string txt;
    while (compress_type?fin_gz.good():fin.good())
    {
        if (compress_type) getline(fin_gz, line);
        else               getline(fin, line);
        L=atoi(line.c_str());
        txt+=line+'\n';
        if (compress_type) getline(fin_gz, line);
        else               getline(fin, line);
        for (i=0;i<line.size();i++)
            if (line[i]==' '||line[i]=='\t') break;
        if (!(compress_type?fin_gz.good():fin.good())) break;
        if (find(chain_list.begin(), chain_list.end(),
            line.substr(0,i))==chain_list.end())
        {
            txt.clear();
            continue;
        }
        txt+=line;
        for (i=0;i<L;i++)
        {
            if (compress_type) getline(fin_gz, line);
            else               getline(fin, line);
            txt+='\n'+line;
        }
        cout<<txt<<endl;
        txt.clear();
    }
    if (compress_type) fin_gz.close();
    else               fin.close();
    /* clean up */
    txt.clear();
    line.clear();
    filename.clear();
    list_opt.clear();
    chain_list.clear();
    return 0;
}
