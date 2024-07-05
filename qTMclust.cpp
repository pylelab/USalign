/* Different filters are used when different header files are included.
 * At least one of HwRMSD.h and TMalign.h should be included.
 * HwRMSD.h implement HwRMSD filter.
 * No filter will be used if only TMalign.h is included. */

#include "HwRMSD.h"
#include "TMalign.h"

using namespace std;

void print_extra_help()
{
    cout <<
"Additional options:\n"
"    -fast    Fast but slightly inaccurate final alignment\n"
"\n"
"    -atom    4-character atom name used to represent a residue.\n"
"             Default is \" C3'\" for RNA/DNA and \" CA \" for proteins\n"
"             (note the spaces before and after CA).\n"
"\n"
"    -mol     Molecule type: RNA or protein\n"
"             Default is detect molecule type automatically\n"
"\n"
"    -het     Whether to align residues marked as 'HETATM' in addition to 'ATOM  '\n"
"             0: (default) only align 'ATOM  ' residues\n"
"             1: align both 'ATOM  ' and 'HETATM' residues\n"
"\n"
"    -infmt   Input format\n"
"            -1: (default) automatically detect PDB or PDBx/mmCIF format\n"
"             0: PDB format\n"
"             1: SPICKER format\n"
"             2: xyz format\n"
"             3: PDBx/mmCIF format\n"
"    -chain   Chains to parse in structure_2. Use _ for a chain without chain ID.\n"
"             Multiple chains can be separated by commas, e.g.,\n"
"             USalign -chain1 C,D,E,F 5jdo.pdb -chain2 A,B,C,D 3wtg.pdb -ter 0\n"
"\n"
    <<endl;
}

void print_help(bool h_opt=false)
{
    cout << "\n"
"qTMclust: Structure Clustering by Sequence-Indepedent Structure Alignment\n"
"\n"
"Usage 1: (alignment within a folder of PDB files)\n"
"    qTMclust -dir chain_folder/ chain_list -TMcut 0.5 -o cluster.txt\n"
"\n"
"Usage 2: (alignment within chains or within models of a single PDB file)\n"
"    qTMclust -split 2 -ter 1 multichain.pdb -TMcut 0.5 -o cluster.txt\n"
"    qTMclust -split 1 -ter 0 multimodel.pdb -TMcut 0.5 -o cluster.txt\n"
"\n"
"Options:\n"
"    -TMcut   TM-score cutoff in the range of [0.45,1) for considering two\n"
"             structures being similar. Default is 0.5.\n"
"\n"
"    -s       Which TM-score to use when aligning structures with different lengths?\n"
"             1: the larger TM-score, i.e. normalized by shorter length\n"
"             2: (default) the smaller TM-score, i.e. normalized by longer length\n"
"             3: average of the two TM-scores\n"
"             4: harmonic average of the two TM-scores\n"
"             5: geometric average of the two TM-scores\n"
"             6: root mean square of the two TM-scores\n"
"\n"
"    -o       Output the cluster result to file.\n"
"             Default is print result to screen.\n"
"\n"
"    -dir     Perform all-against-all alignment among the list of PDB\n"
"             chains listed by 'chain_list' under 'chain_folder'. Note\n"
"             that the slash is necessary.\n"
"             $ qTMclust -dir chain_folder/ chain_list\n"
"\n"
"    -suffix  (Only when -dir is set, default is empty)\n"
"             add file name suffix to files listed by chain_list\n"
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
"    -init    tentative clustering\n"
"\n"
"    -h       Print the full help message, including additional options.\n"
"\n"
    <<endl;

    if (h_opt) print_extra_help();

    exit(EXIT_SUCCESS);
}

void filter_lower_bound(double &lb_HwRMSD, double &lb_TMfast, 
    const double TMcut, const int s_opt,const int mol_type)
{
    lb_HwRMSD=0.5*TMcut;
    lb_TMfast=0.9*TMcut;
    if (s_opt<=1)
    {
        if (mol_type>0) // RNA
        {
            lb_HwRMSD=0.02*TMcut;
            lb_TMfast=0.60*TMcut;
        }
        else // protein
        {
            lb_HwRMSD=0.25*TMcut;
            lb_TMfast=0.80*TMcut;
        }
    }
    return;
}

void read_init_cluster(const string&filename, 
    map<string, map<string,bool> > &init_cluster)
{
    ifstream fin;
    string line;
    vector<string> line_vec;
    map<string, bool> tmp_map;
    size_t i,j;
    fin.open(filename.c_str());
    while (fin.good())
    {
        getline(fin,line);
        split(line,line_vec,'\t');
        for (i=0;i<line_vec.size();i++)
        {
            for (j=0;j<line_vec.size();j++)
                if (i!=j) tmp_map[line_vec[j]]=1;
            init_cluster[line_vec[i]]=tmp_map;
            map<string, bool> ().swap(tmp_map);
        }
        for (i=0;i<line_vec.size();i++) line_vec[i].clear(); line_vec.clear();
    }
    fin.close();
    vector<string>().swap(line_vec);
}

int main(int argc, char *argv[])
{
    if (argc < 2) print_help();


    clock_t t1, t2;
    t1 = clock();

    /**********************/
    /*    get argument    */
    /**********************/
    string xname       = "";
    double TMcut       = 0.5;
    string fname_clust = ""; // file name for output cluster result
    string fname_init  = "";
    string fname_lign  = ""; // file name for user alignment
    vector<string> sequence; // get value from alignment file
    double Lnorm_ass, d0_scale;

    bool h_opt = false; // print full help message
    int  i_opt = 0;     // 3 for -I, stick to user given alignment
    int  a_opt = 0;     // flag for -a, do not normalized by average length
    int  s_opt = 2;     // flag for -s, normalized by longer length
    bool u_opt = false; // flag for -u, normalized by user specified length
    bool d_opt = false; // flag for -d, user specified d0

    int    infmt_opt =-1;    // PDB or PDBx/mmCIF format
    int    ter_opt   =3;     // TER, END, or different chainID
    int    split_opt =0;     // do not split chain
    bool   fast_opt  =false; // flags for -fast, fTM-align algorithm
    int    het_opt   =0;     // do not read HETATM residues
    string atom_opt  ="auto";// use C alpha atom for protein and C3' for RNA
    string mol_opt   ="auto";// auto-detect the molecule type as protein/RNA
    string suffix_opt="";    // set -suffix to empty
    string dir_opt   ="";    // set -dir to empty
    int    byresi_opt=0;     // set -byresi to 0
    vector<string> chain_list;
    vector<string> chain2parse;
    vector<string> model2parse;
    map<string, map<string,bool> > init_cluster;

    for(int i = 1; i < argc; i++)
    {
        if ( (!strcmp(argv[i],"-u")||!strcmp(argv[i],"-L")) && i < (argc-1) )
        {
            PrintErrorAndQuit("Sorry! -u has not been implemented yet");
            Lnorm_ass = atof(argv[i + 1]); u_opt = true; i++;
        }
        else if ( !strcmp(argv[i],"-d") && i < (argc-1) )
        {
            PrintErrorAndQuit("Sorry! -d has not been implemented yet");
            d0_scale = atof(argv[i + 1]); d_opt = true; i++;
        }
        else if (!strcmp(argv[i], "-I") && i < (argc-1) )
        {
            fname_lign = argv[i + 1];      i_opt = 3; i++;
        }
        else if ( !strcmp(argv[i],"-o") && i < (argc-1) )
        {
            fname_clust = argv[i + 1]; i++;
        }
        else if ( !strcmp(argv[i],"-a") && i < (argc-1))
        {
            PrintErrorAndQuit("Sorry! -a is not used for clustering");
        }
        else if ( !strcmp(argv[i],"-s") && i < (argc-1) )
        {
            s_opt=atoi(argv[i + 1]); i++;
            if (s_opt<1 || s_opt>6)
                PrintErrorAndQuit("-s must be within 1 to 6");
        }
        else if ( !strcmp(argv[i],"-h") )
        {
            h_opt = true;
        }
        else if (!strcmp(argv[i], "-fast"))
        {
            fast_opt = true;
        }
        else if ( !strcmp(argv[i],"-infmt") && i < (argc-1) )
        {
            infmt_opt=atoi(argv[i + 1]); i++;
        }
        else if ( !strcmp(argv[i],"-ter") && i < (argc-1) )
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
        else if ( !strcmp(argv[i],"-TMcut") && i < (argc-1) )
        {
            TMcut=atof(argv[i + 1]); i++;
            if (TMcut>1 or TMcut<0.45)
                PrintErrorAndQuit("TMcut must be in the range of [0.45,1)");
        }
        else if ( !strcmp(argv[i],"-byresi") && i < (argc-1) )
        {
            PrintErrorAndQuit("Sorry! -byresi has not been implemented yet");
            byresi_opt=atoi(argv[i + 1]); i++;
        }
        else if ( !strcmp(argv[i],"-het") && i < (argc-1) )
        {
            het_opt=atoi(argv[i + 1]); i++;
        }
        else if ( !strcmp(argv[i],"-init") && i < (argc-1) )
        {
            read_init_cluster(argv[i+1],init_cluster); i++;
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
        else if (xname.size() == 0) xname=argv[i];
        else PrintErrorAndQuit(string("ERROR! Undefined option ")+argv[i]);
    }

    if(xname.size()==0) print_help(h_opt);

    if (suffix_opt.size() && dir_opt.size()==0)
        PrintErrorAndQuit("-suffix is only valid if -dir, -dir1 or -dir2 is set");
    if (atom_opt.size()!=4)
        PrintErrorAndQuit("ERROR! Atom name must have 4 characters, including space.");
    if (mol_opt!="auto" && mol_opt!="protein" && mol_opt!="RNA")
        PrintErrorAndQuit("ERROR! Molecule type must be either RNA or protein.");
    else if (mol_opt=="protein" && atom_opt=="auto")
        atom_opt=" CA ";
    else if (mol_opt=="RNA" && atom_opt=="auto")
        atom_opt=" C3'";

    if (u_opt && Lnorm_ass<=0)
        PrintErrorAndQuit("Wrong value for option -u!  It should be >0");
    if (d_opt && d0_scale<=0)
        PrintErrorAndQuit("Wrong value for option -d!  It should be >0");
    if (split_opt==1 && ter_opt!=0)
        PrintErrorAndQuit("-split 1 should be used with -ter 0");
    else if (split_opt==2 && ter_opt!=0 && ter_opt!=1)
        PrintErrorAndQuit("-split 2 should be used with -ter 0 or 1");
    if (split_opt<0 || split_opt>2)
        PrintErrorAndQuit("-split can only be 0, 1 or 2");

    /* read initial alignment file from 'align.txt' */
    if (i_opt) read_user_alignment(sequence, fname_lign, i_opt);

    if (byresi_opt) i_opt=3;

    /* parse file list */
    if (dir_opt.size()==0) chain_list.push_back(xname);
    else file2chainlist(chain_list, xname, dir_opt, suffix_opt);

    /* declare previously global variables */
    vector<vector<string> >PDB_lines; // text of chain
    vector<int>    mol_vec;           // molecule type of chain1, RNA if >0
    vector<string> chainID_list;      // list of chainID
    size_t xchainnum=0;         // number of chains in a PDB file
    size_t i,j;                 // number of residues/chains in a PDB is
                                // usually quite limited. Yet, the number of
                                // files can be very large. size_t is safer
                                // than int for very long list of files
    int    xlen,ylen;           // chain length
    double **xa,**ya;           // xyz coordinate
    vector<string> resi_vec;    // residue index for chain, dummy variable
    vector<pair<int,size_t> >chainLen_list; // vector of (length,index) pair
    vector<vector<char> > seq_vec;
    vector<vector<char> > sec_vec;
    vector<vector<vector<float> > >xyz_vec;

    /* parse files */
    string chain_name;
    vector<char>  seq_tmp;
    vector<char>  sec_tmp;
    vector<float> flt_tmp(3,0);
    vector<vector<float> >xyz_tmp;
    int r; // residue index
    size_t newchainnum;
    double ub_HwRMSD=0.90*TMcut+0.10;
    double lb_HwRMSD=0.5*TMcut;
    double ub_TMfast=0.90*TMcut+0.10;
    double lb_TMfast=0.9*TMcut;
    if      (s_opt==2 || s_opt==4 || s_opt==5) a_opt=-2; // normalized by longer length, i.e. smaller TM
    else if (s_opt==1 || s_opt==5) a_opt=-1; // normalized by shorter length, i.e. larger TM
    else if (s_opt==3) a_opt= 1; // normalized by average length

#ifdef TMalign_HwRMSD_h
    /* These parameters controls HwRMSD filter. iter_opt typically should be
     * >=3. Many alignments converge within iter_opt=5. Occassionally
     * some alignments require iter_opt=10. Higher iter_opt takes more time,
     * even though HwRMSD iter_opt 10 still takes far less time than TMalign
     * -fast -TMcut 0.5.
     * After HwRMSD filter, at least min_repr_num and at most max_repr_num
     * are used for subsequent TMalign. The actual number of representatives
     * are decided by xlen */
    const int glocal    =0; // global alignment
    const int iter_opt  =10;
    const int min_repr_num=10;
    const int max_repr_num=50;
#endif

    for (i=0;i<chain_list.size();i++)
    {
        xname=chain_list[i];
        newchainnum=get_PDB_lines(xname, PDB_lines, chainID_list,
            mol_vec, ter_opt, infmt_opt, atom_opt, false, split_opt, het_opt,
            chain2parse, model2parse);
        if (!newchainnum)
        {
            cerr<<"Warning! Cannot parse file: "<<xname
                <<". Chain number 0."<<endl;
            continue;
        }
        chain_name=xname.substr(dir_opt.size(),
            xname.size()-dir_opt.size()-suffix_opt.size());
        for (j=0;j<newchainnum;j++)
        {
            chainID_list[j+xchainnum]=chain_name+chainID_list[j+xchainnum];
            xlen=PDB_lines[j].size();
            cout<<"Parsing "<<xname<<'\t'<<chainID_list[j+xchainnum]
                <<" ("<<xlen<<" residues)."<<endl;
            if (mol_opt=="RNA") mol_vec[j+xchainnum]=1;
            else if (mol_opt=="protein") mol_vec[j+xchainnum]=-1;

            NewArray(&xa, xlen, 3);
            seq_tmp.assign(xlen+1,'A');
            sec_tmp.assign(xlen+1,0);

            read_PDB(PDB_lines[j], xa, &seq_tmp[0], resi_vec, byresi_opt);

            if (mol_vec[j]<=0) make_sec(xa, xlen, &sec_tmp[0]);
            else make_sec(&seq_tmp[0],xa,xlen,&sec_tmp[0],atom_opt);

            xyz_tmp.assign(xlen,flt_tmp);
            for (r=0;r<xlen;r++)
            {
                xyz_tmp[r][0]=xa[r][0];
                xyz_tmp[r][1]=xa[r][1];
                xyz_tmp[r][2]=xa[r][2];
            }

            seq_vec.push_back(seq_tmp);
            sec_vec.push_back(sec_tmp);
            xyz_vec.push_back(xyz_tmp);

            chainLen_list.push_back(
                make_pair(PDB_lines[j].size(),j+xchainnum));

            seq_tmp.clear();
            sec_tmp.clear();
            xyz_tmp.clear();
            DeleteArray(&xa, xlen);
            PDB_lines[j].clear();
        }
        PDB_lines.clear();
        xchainnum+=newchainnum;
    }
    flt_tmp.clear();
    chain_list.clear();

    // swap completely destroy the vector and free up the memory capacity
    vector<vector<string> >().swap(PDB_lines);
    size_t Nstruct=chainLen_list.size();

    /* sort by chain length */
    stable_sort(chainLen_list.begin(),chainLen_list.end(),
        greater<pair<int,int> >());
    cout<<"Clustering "<<chainLen_list.size()
        <<" chains with TM-score cutoff >="<<TMcut<<'\n'
        <<"Longest chain "<<chainID_list[chainLen_list[0].second]<<'\t'
        <<chainLen_list[0].first<<" residues.\n"
        <<"Shortest chain "<<chainID_list[chainLen_list.back().second]<<'\t'
        <<chainLen_list.back().first<<" residues."<<endl;

    /* set the first cluster */
    vector<size_t> clust_mem_vec(Nstruct,-1); // cluster membership
    vector<size_t> clust_repr_vec; // the same as number of clusters
    size_t chain_i=chainLen_list[0].second;
    clust_repr_vec.push_back(chain_i);
    clust_mem_vec[chain_i]=0;
    map<size_t,size_t> clust_repr_map;

    /* perform alignment */
    size_t chain_j;
    const double fast_lb=50.;  // proteins shorter than fast_lb never use -fast
    const double fast_ub=1000.;// proteins longer than fast_ub always use -fast
    double Lave;               // average protein length for chain_i and chain_j
    size_t sizePROT;           // number of representatives for current chain
    vector<size_t> index_vec;  // index of cluster representatives for the chain
    bool found_clust;          // whether current chain hit previous cluster

    for (i=1;i<Nstruct;i++)
    {
        chain_i=chainLen_list[i].second;
        xlen=xyz_vec[chain_i].size();
        if (xlen<=5) // TMalign cannot handle L<=5
        {
            clust_mem_vec[chain_i]=clust_repr_vec.size();
            clust_repr_vec.push_back(clust_repr_vec.size());
            continue;
        }

        NewArray(&xa, xlen, 3);
        for (r=0;r<xlen;r++)
        {
            xa[r][0]=xyz_vec[chain_i][r][0];
            xa[r][1]=xyz_vec[chain_i][r][1];
            xa[r][2]=xyz_vec[chain_i][r][2];
        }

        // j-1 is index of old cluster. here, we starts from the latest
        // cluster because proteins with similar length are more likely
        // to be similar. we cannot use j as index because size_t j cannot
        // be negative at the end of this loop
        for (j=clust_repr_vec.size();j>0;j--)
        {
            chain_j=clust_repr_vec[j-1];
            ylen=xyz_vec[chain_j].size();
            if (mol_vec[chain_i]*mol_vec[chain_j]<0)    continue;
            else if (s_opt==2 && xlen<TMcut*ylen)       continue;
            else if (s_opt==3 && xlen<(2*TMcut-1)*ylen) continue;
            else if (s_opt==4 && xlen*(2/TMcut-1)<ylen) continue;
            else if (s_opt==5 && xlen<TMcut*TMcut*ylen) continue;
            else if (s_opt==6 && xlen*xlen<(2*TMcut*TMcut-1)*ylen*ylen) continue;
            index_vec.push_back(chain_j);
        }
        sizePROT=index_vec.size();

        string key=chainID_list[chain_i];
        cout<<'>'<<chainID_list[chain_i]<<'\t'<<xlen<<'\t'
            <<setiosflags(ios::fixed)<<setprecision(2)
            <<100.*i/Nstruct<<"%(#"<<i<<")\t"
            <<"#repr="<<sizePROT<<"/"<<clust_repr_vec.size()<<endl;

#ifdef TMalign_HwRMSD_h
        vector<pair<double,size_t> > HwRMSDscore_list;
        double TM;
        size_t init_count=0;
        for (j=0;j<sizePROT;j++)
        {
            chain_j=index_vec[j];
            string value=chainID_list[chain_j];
            if (init_cluster.count(key) && init_count>=2 && 
                HwRMSDscore_list.size()>=init_cluster[key].size() && !init_cluster[key].count(value))
                continue;
            ylen=xyz_vec[chain_j].size();
            if (mol_vec[chain_i]*mol_vec[chain_j]<0)    continue;
            else if (s_opt==2 && xlen<TMcut*ylen)       continue;
            else if (s_opt==3 && xlen<(2*TMcut-1)*ylen) continue;
            else if (s_opt==4 && xlen*(2/TMcut-1)<ylen) continue;
            else if (s_opt==5 && xlen<TMcut*TMcut*ylen) continue;
            else if (s_opt==6 && xlen*xlen<(2*TMcut*TMcut-1)*ylen*ylen) continue;

            if (s_opt<=1) filter_lower_bound(lb_HwRMSD, lb_TMfast, 
                TMcut, s_opt, mol_vec[chain_i]+mol_vec[chain_j]);
            
            //cout<<chainID_list[chain_i]<<" => "<<chainID_list[chain_j]<<endl;
            
            NewArray(&ya, ylen, 3);
            for (r=0;r<ylen;r++)
            {
                ya[r][0]=xyz_vec[chain_j][r][0];
                ya[r][1]=xyz_vec[chain_j][r][1];
                ya[r][2]=xyz_vec[chain_j][r][2];
            }

            /* declare variable specific to this pair of HwRMSD */
            double t0[3], u0[3][3];
            double TM1, TM2;
            double TM3, TM4, TM5;     // for s_opt, u_opt, d_opt
            double d0_0, TM_0;
            double d0A, d0B, d0u, d0a;
            double d0_out=5.0;
            string seqM, seqxA, seqyA;// for output alignment
            double rmsd0 = 0.0;
            int L_ali;                // Aligned length in standard_TMscore
            double Liden=0;
            double TM_ali, rmsd_ali;  // TMscore and rmsd in standard_TMscore
            int n_ali=0;
            int n_ali8=0;
            int *invmap = new int[ylen+1];

            /* entry function for structure alignment */
            HwRMSD_main(
                xa, ya, &seq_vec[chain_i][0], &seq_vec[chain_j][0],
                &sec_vec[chain_i][0], &sec_vec[chain_j][0], t0, u0,
                TM1, TM2, TM3, TM4, TM5,
                d0_0, TM_0, d0A, d0B, d0u,
                d0a, d0_out, seqM, seqxA, seqyA,
                rmsd0, L_ali, Liden, TM_ali,
                rmsd_ali, n_ali, n_ali8, xlen, ylen,
                sequence, Lnorm_ass,
                d0_scale, i_opt,
                a_opt, u_opt, d_opt, mol_vec[chain_i]+mol_vec[chain_j],
                invmap, glocal, iter_opt);

            TM=TM3; // average length
            if      (s_opt==1) TM=TM2; // shorter length
            else if (s_opt==2) TM=TM1; // longer length
            else if (s_opt==3) TM=(TM1+TM2)/2;     // average TM
            else if (s_opt==4) TM=2/(1/TM1+1/TM2); // harmonic average
            else if (s_opt==5) TM=sqrt(TM1*TM2);   // geometric average
            else if (s_opt==6) TM=sqrt((TM1*TM1+TM2*TM2)/2); // root mean square

            Lave=sqrt(xlen*ylen); // geometry average because O(L1*L2)
            if (TM>=lb_HwRMSD || Lave<=fast_lb)
            {
                if (init_cluster.count(key) && init_cluster[key].count(value))
                {
                    HwRMSDscore_list.push_back(make_pair(TM+1,index_vec[j]));
                    init_count++;
                    if (init_count==init_cluster[key].size()) break;
                }
                else
                    HwRMSDscore_list.push_back(make_pair(TM,index_vec[j]));
            }

            /* clean up after each HwRMSD */
            seqM.clear();
            seqxA.clear();
            seqyA.clear();
            DeleteArray(&ya, ylen);
            delete [] invmap;

            /* if a good hit is guaranteed to be found, stop the loop */
            if (TM>=ub_HwRMSD) break;
        }

        stable_sort(HwRMSDscore_list.begin(),HwRMSDscore_list.end(),
            greater<pair<double,size_t> >());

        int cur_repr_num_cutoff=min_repr_num;
        if (xlen<=fast_lb) cur_repr_num_cutoff=max_repr_num;
        else if (xlen>fast_lb && xlen<fast_ub) cur_repr_num_cutoff+=
            (fast_ub-xlen)/(fast_ub-fast_lb)*(max_repr_num-min_repr_num);
        //if (init_count>=2) cur_repr_num_cutoff=init_count;

        index_vec.clear();
        for (j=0;j<HwRMSDscore_list.size();j++)
        {
            TM=HwRMSDscore_list[j].first;
            chain_j=HwRMSDscore_list[j].second;
            ylen=xyz_vec[chain_j].size();
            Lave=sqrt(xlen*ylen); // geometry average because O(L1*L2)
            if (Lave>fast_lb && TM<TMcut*0.5 && 
                index_vec.size()>=cur_repr_num_cutoff) break;
            index_vec.push_back(chain_j);
            cout<<"#"<<chain_j<<"\t"<<chainID_list[chain_j]<<"\t"
                <<setiosflags(ios::fixed)<<setprecision(4)<<TM<<endl;
        }
        cout<<index_vec.size()<<" out of "
            <<HwRMSDscore_list.size()<<" entries"<<endl;
        HwRMSDscore_list.clear();
#endif

        found_clust=false;
        for (j=0;j<index_vec.size();j++)
        {
            chain_j=index_vec[j];
            ylen=xyz_vec[chain_j].size();
            if (mol_vec[chain_i]*mol_vec[chain_j]<0)    continue;
            else if (s_opt==2 && xlen<TMcut*ylen)       continue;
            else if (s_opt==3 && xlen<(2*TMcut-1)*ylen) continue;
            else if (s_opt==4 && xlen*(2/TMcut-1)<ylen) continue;
            else if (s_opt==5 && xlen<TMcut*TMcut*ylen) continue;
            else if (s_opt==6 && xlen*xlen<(2*TMcut*TMcut-1)*ylen*ylen) continue;
            if (s_opt<=1) filter_lower_bound(lb_HwRMSD, lb_TMfast,
                TMcut, s_opt, mol_vec[chain_i]+mol_vec[chain_j]);

            NewArray(&ya, ylen, 3);
            for (r=0;r<ylen;r++)
            {
                ya[r][0]=xyz_vec[chain_j][r][0];
                ya[r][1]=xyz_vec[chain_j][r][1];
                ya[r][2]=xyz_vec[chain_j][r][2];
            }

            Lave=sqrt(xlen*ylen); // geometry average because O(L1*L2)
            bool overwrite_fast_opt=(fast_opt==true || Lave>=fast_ub);
            
            /* declare variable specific to this pair of TMalign */
            double t0[3], u0[3][3];
            double TM1, TM2;
            double TM3, TM4, TM5;     // for s_opt, u_opt, d_opt
            double d0_0, TM_0;
            double d0A, d0B, d0u, d0a;
            double d0_out=5.0;
            string seqM, seqxA, seqyA;// for output alignment
            double rmsd0 = 0.0;
            int L_ali;                // Aligned length in standard_TMscore
            double Liden=0;
            double TM_ali, rmsd_ali;  // TMscore and rmsd in standard_TMscore
            int n_ali=0;
            int n_ali8=0;
            vector<double> do_vec;
            
            /* entry function for structure alignment */
            int status=TMalign_main(
                xa, ya, &seq_vec[chain_i][0], &seq_vec[chain_j][0],
                &sec_vec[chain_i][0], &sec_vec[chain_j][0],
                t0, u0, TM1, TM2, TM3, TM4, TM5,
                d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                seqM, seqxA, seqyA, do_vec,
                rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                xlen, ylen, sequence, Lnorm_ass, d0_scale,
                i_opt, a_opt, u_opt, d_opt, overwrite_fast_opt,
                mol_vec[chain_i]+mol_vec[chain_j],TMcut);

            cout<<status<<'\t'<<chainID_list[chain_j]<<'\t'
                <<setiosflags(ios::fixed)<<setprecision(4)
                <<TM2<<'\t'<<TM1<<'\t'<<overwrite_fast_opt<<endl;

            seqM.clear();
            seqxA.clear();
            seqyA.clear();
            do_vec.clear();

            double TM=TM3; // average length
            if      (s_opt==1) TM=TM2; // shorter length
            else if (s_opt==2) TM=TM1; // longer length
            else if (s_opt==3) TM=(TM1+TM2)/2;     // average TM
            else if (s_opt==4) TM=2/(1/TM1+1/TM2); // harmonic average
            else if (s_opt==5) TM=sqrt(TM1*TM2);   // geometric average
            else if (s_opt==6) TM=sqrt((TM1*TM1+TM2*TM2)/2); // root mean square

            if (TM<lb_TMfast || 
               (TM<TMcut && (fast_opt || overwrite_fast_opt==false)))
            {
                DeleteArray(&ya, ylen);
                continue;
            }

            if (TM>=ub_TMfast || 
               (TM>=TMcut && (fast_opt || overwrite_fast_opt==false)))
            {
                clust_mem_vec[chain_i]=clust_repr_map[chain_j];
                DeleteArray(&ya, ylen);
                found_clust=true;
                break;
            }

            if (TM<lb_TMfast && overwrite_fast_opt==false)
            {
                TMalign_main(
                    xa, ya, &seq_vec[chain_i][0], &seq_vec[chain_j][0],
                    &sec_vec[chain_i][0], &sec_vec[chain_j][0],
                    t0, u0, TM1, TM2, TM3, TM4, TM5,
                    d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                    seqM, seqxA, seqyA, do_vec,
                    rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                    xlen, ylen, sequence, Lnorm_ass, d0_scale,
                    i_opt, a_opt, u_opt, d_opt, false,
                    mol_vec[chain_i]+mol_vec[chain_j],TMcut);
                seqM.clear();
                seqxA.clear();
                seqyA.clear();
                do_vec.clear();
                DeleteArray(&ya, ylen);
                
                TM=TM3;                // average length
                if      (s_opt==1) TM=TM2; // shorter length
                else if (s_opt==2) TM=TM1; // longer length
                else if (s_opt==3) TM=(TM1+TM2)/2;     // average TM
                else if (s_opt==4) TM=2/(1/TM1+1/TM2); // harmonic average
                else if (s_opt==5) TM=sqrt(TM1*TM2);   // geometric average
                else if (s_opt==6) TM=sqrt((TM1*TM1+TM2*TM2)/2); // root mean square
                cout<<"*\t"<<chainID_list[chain_j]<<'\t'<<TM2<<'\t'<<TM1<<endl;
                if (TM>=TMcut)
                {
                    clust_mem_vec[chain_i]=clust_repr_map[chain_j];
                    found_clust=true;
                    break;
                }
            }
        }
        DeleteArray(&xa, xlen);
        index_vec.clear();

        if (!found_clust) // new cluster
        {
            clust_mem_vec[chain_i]=clust_repr_vec.size();
            clust_repr_map[chain_i]=clust_repr_vec.size();
            clust_repr_vec.push_back(chain_i);
        }
        else // member structures are not used further
        {
            vector<char> ().swap(seq_vec[chain_i]);
            vector<char> ().swap(sec_vec[chain_i]);
            vector<vector<float> > ().swap(xyz_vec[chain_i]);
        }
    }

    /* clean up */
    mol_vec.clear();
    xyz_vec.clear();
    seq_vec.clear();
    sec_vec.clear();

    /* print out cluster */
    stringstream txt;
    for (j=0;j<clust_repr_vec.size();j++)
    {
        chain_j=clust_repr_vec[j]; // cluster representative
        txt<<chainID_list[chain_j];
        for (chain_i=0;chain_i<clust_mem_vec.size();chain_i++)
        {
            if (chain_i!=chain_j && clust_mem_vec[chain_i]==j)
                txt<<'\t'<<chainID_list[chain_i];
        }
        txt<<'\n';
    }
    if (fname_clust.size() && fname_clust!="-")
    {
        ofstream fp(fname_clust.c_str());
        fp<<txt.str();
        fp.close();
    }
    else cout<<txt.str()<<endl;

    /* clean up */
    txt.str(string());
    clust_repr_vec.clear();
    clust_mem_vec.clear();
    chainID_list.clear();
    clust_repr_map.clear();
    vector<string>().swap(chain2parse);
    vector<string>().swap(model2parse);
    map<string, map<string,bool> >().swap(init_cluster);

    t2 = clock();
    float diff = ((float)t2 - (float)t1)/CLOCKS_PER_SEC;
    printf("#Total CPU time is %5.2f seconds\n", diff);
    return 0;
}
