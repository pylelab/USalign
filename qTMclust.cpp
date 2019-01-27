//#include "Contactlib.h"
#include "TMalign.h"

using namespace std;

void print_extra_help()
{
    cout <<
"Additional options:\n"
"    -fast    Fast but slightly inaccurate final alignment\n"
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
"    -infmt   Input format\n"
"            -1: (default) automatically detect PDB or PDBx/mmCIF format\n"
"             0: PDB format\n"
"             1: SPICKER format\n"
"             2: xyz format\n"
"             3: PDBx/mmCIF format\n"
    <<endl;
}

void print_help(bool h_opt=false)
{
    cout << "\n"
"qTMclust: Structure Clustering by Sequence-Indepedent Structure Alignment\n"
"\n"
"Usage 1: (alignment within a folder of PDB files)\n"
"    qTMclust -dir chain_folder/ chain_list TMcut -o cluster.txt\n"
"\n"
"Usage 2: (alignment within chains or within models of a single PDB file)\n"
"    qTMclust -split 2 -ter 1 multichain.pdb TMcut -o cluster.txt\n"
"    qTMclust -split 1 -ter 0 multimodel.pdb TMcut -o cluster.txt\n"
"\n"
"Options:\n"
"    TMcut    TM-score cutoff in the range of [0.5,1) for considering two\n"
"             structures being similar.\n"
"\n"
"    -a       How to normalize TM-score:\n"
"             0: (default) normalize TM-score by longer protein length\n"
"            -1: normalize TM-score by shorter protein length\n"
"             1: normalize TM-score by average protein length\n"
"\n"
"    -o       Output the cluster result to file.\n"
"             Default is print result to screen.\n"
"\n"
"    -dir     Perform all-against-all alignment among the list of PDB\n"
"             chains listed by 'chain_list' under 'chain_folder'. Note\n"
"             that the slash is necessary.\n"
"             $ qTMclust -dir chain_folder/ chain_list\n"
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
"    -h       Print the full help message, including additional options.\n"
"\n"
    <<endl;

    if (h_opt) print_extra_help();

    exit(EXIT_SUCCESS);
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
    string fname_lign  = ""; // file name for user alignment
    vector<string> sequence; // get value from alignment file
    double Lnorm_ass, d0_scale;

    bool h_opt = false; // print full help message
    bool I_opt = false; // flag for -I, stick to user given alignment
    int  a_opt = 0;     // flag for -a, normalized by longer length
    bool u_opt = false; // flag for -u, normalized by user specified length
    bool d_opt = false; // flag for -d, user specified d0

    int    infmt_opt =-1;    // PDB or PDBx/mmCIF format
    int    ter_opt   =3;     // TER, END, or different chainID
    int    split_opt =0;     // do not split chain
    bool   fast_opt  =false; // flags for -fast, fTM-align algorithm
    string atom_opt  ="auto";// use C alpha atom for protein and C3' for RNA
    string mol_opt   ="auto";// auto-detect the molecule type as protein/RNA
    string suffix_opt="";    // set -suffix to empty
    string dir_opt   ="";    // set -dir to empty
    int    byresi_opt=0;     // set -byresi to 0
    vector<string> chain_list;

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
            PrintErrorAndQuit("Sorry! -I has not been implemented yet");
            fname_lign = argv[i + 1];      I_opt = true; i++;
        }
        else if ( !strcmp(argv[i],"-o") && i < (argc-1) )
        {
            fname_clust = argv[i + 1]; i++;
        }
        else if ( !strcmp(argv[i],"-a") && i < (argc-1) )
        {
            a_opt=atoi(argv[i + 1]); i++;
            if (a_opt!=0 && a_opt!=-1 && a_opt!=1)
                PrintErrorAndQuit("-a must be -1, 0, or 1");
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
        else if ( !strcmp(argv[i],"-byresi") && i < (argc-1) )
        {
            PrintErrorAndQuit("Sorry! -byresi has not been implemented yet");
            byresi_opt=atoi(argv[i + 1]); i++;
        }
        else if (xname.size() == 0) xname=argv[i];
        else
        {
            TMcut=atof(argv[i]);
            if (TMcut>1 or TMcut<0.5)
                PrintErrorAndQuit("TMcut must be in the range of [0.5,1)");
        }
    }

    if(xname.size()==0) print_help(h_opt);

    if (suffix_opt.size() && dir_opt.size()==0)
        PrintErrorAndQuit("-suffix is only valid if -dir, -dir1 or -dir2 is set");
    if (atom_opt.size()!=4)
        PrintErrorAndQuit("ERROR! atom name must have 4 characters, including space.");
    if (mol_opt!="auto" && mol_opt!="protein" && mol_opt!="RNA")
        PrintErrorAndQuit("ERROR! molecule type must be either RNA or protein.");
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
    if (I_opt) read_user_alignment(sequence, fname_lign, I_opt);

    if (byresi_opt) I_opt=true;

    /* parse file list */
    if (dir_opt.size()==0) chain_list.push_back(xname);
    else file2chainlist(chain_list, xname, dir_opt, suffix_opt);

    /* declare previously global variables */
    vector<vector<string> >PDB_lines; // text of chain
    vector<int>    mol_vec;           // molecule type of chain1, RNA if >0
    vector<string> chainID_list;      // list of chainID
    int    xchainnum=0;         // number of chains in a PDB file
    int    i,j;                 // file index
    int    xlen,ylen;           // chain length
    double **xa,**ya;           // xyz coordinate
    vector<string> resi_vec;    // residue index for chain, dummy variable
    vector<pair<int,int> >chainLen_list; // vector of (length,index) pair
    vector<vector<char> >seq_vec;
    vector<vector<int> >sec_vec;
    vector<vector<vector<float> > >xyz_vec;

    /* parse files */
    string chain_name;
    vector<char> seq_tmp;
    vector<int>  sec_tmp;
    vector<float> flt_tmp(3,0);
    vector<vector<float> >xyz_tmp;
    int r; // residue index
    int newchainnum;

#ifdef TMalign_Contactlib_h
    vector<vector<vector<int> > >vindexLocal_vec;
    vector<vector<vector<int> > >vindexRemote_vec;
    vector<vector<int> >vindexLocal;
    vector<vector<int> >vindexRemote;
    vector<int> SeqNum; // residue index in int, dummy variable
    map<int, char> ssDSSP;

    /* These parameters controls Contactlib filter. building and searching
     * Contactlib database also takes some time. Therefore, if there are
     * only a handful of existing cluster representatives, or the new chain
     * is very short and thus unlikely to have much contact groups, the
     * filter is skipped. Note that Contactlib database could be very large.
     * In this case, it might be necessary to split the database. */
    const int repr_num_cutoff=100; // The number of remaining representatives
                                   // after filter
    const int sizePROT_cutoff=1000;// minimum number of representatives
                                   // above which Contactlib filter is used.
    const int xlen_cutoff=100;     // maximum number of residues below which
                                   // Contactlib filter is not used
#endif

    for (i=0;i<chain_list.size();i++)
    {
        xname=chain_list[i];
        newchainnum=get_PDB_lines(xname, PDB_lines, chainID_list,
            mol_vec, ter_opt, infmt_opt, atom_opt, split_opt);
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
            chainLen_list.push_back(
                make_pair(PDB_lines[j].size(),j+xchainnum));
            if (mol_opt=="RNA") mol_vec[j+xchainnum]=1;
            else if (mol_opt=="protein") mol_vec[j+xchainnum]=-1;

            NewArray(&xa, xlen, 3);
            seq_tmp.assign(xlen+1,'A');
            sec_tmp.assign(xlen,0);

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

#ifdef TMalign_Contactlib_h
            /* make contactlib*/
            if (xlen>xlen_cutoff)
            {
                SeqNum.assign(xlen,0);
                for (r=0;r<xlen;r++)
                {
                    SeqNum[r]=r;
                    ssDSSP[r]=sec_tmp[r];
                }
                ContactlibBuild_main(xlen, xa, SeqNum, ssDSSP, sizeLocal,
                    sizeRemote, mingapRemote, vindexLocal, vindexRemote);
                SeqNum.clear();
                ssDSSP.clear();
            }
            vindexLocal_vec.push_back(vindexLocal);
            vindexRemote_vec.push_back(vindexRemote);

            /* clean up */
            vindexLocal.clear();
            vindexRemote.clear();
#endif
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
    int Nstruct=chainLen_list.size();
#ifdef TMalign_Contactlib_h
    if (Nstruct<=sizePROT_cutoff)
    {
        vector<vector<vector<int> > >().swap(vindexLocal_vec);
        vector<vector<vector<int> > >().swap(vindexRemote_vec);
    }
#endif

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
    vector<int> clust_mem_vec(Nstruct,-1); // cluster membership
    vector<int> clust_repr_vec; // the same as number of clusters
    int chain_i=chainLen_list[0].second;
    clust_repr_vec.push_back(chain_i);
    clust_mem_vec[chain_i]=0;
    map<int,int> clust_repr_map;

    /* perform alignment */
    int chain_j;
    const double fast_lb=50.;   // proteins shorter than fast_lb never use -fast
    const double fast_ub=1000.; // proteins longer than fast_ub always use -fast
    double Lave;                // average protein length for chain_i and chain_j
    int sizePROT;          // number of representatives for current chain
    vector<int> index_vec; // index of cluster representatives for this chain

#ifdef TMalign_Contactlib_h
    vector<int> psizeLocal; // local contact group number of each protein
    vector<int> psizeRemote;// remote contact group number of each protein
    vector<pair<double,int> > contactlibScore_list;
    double scoreLocal1,  scoreLocal2; // 1 for target, 2 for database
    double scoreRemote1, scoreRemote2;
    double scoreLocal,   scoreRemote, scoreCombine;
    int   tidsizeLocal, tidsizeRemote;
#endif

    for (i=1;i<Nstruct;i++)
    {
        chain_i=chainLen_list[i].second;
        xlen=chainLen_list[i].first;
        if (clust_mem_vec[chain_i]>=0) continue;
        else if (xlen<=5) // TMalign cannot handle L<=5
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

        // j is index of old cluster. here, we starts from the latest
        // cluster because proteins with similar length are more likely
        // to be similar
        for (j=clust_repr_vec.size()-1;j>=0;j--)
        {
            chain_j=clust_repr_vec[j];
            ylen=xyz_vec[chain_j].size();
            if (a_opt==0 && xlen<TMcut*ylen) continue;
            else if (a_opt==1 && xlen<TMcut*(xlen+ylen)/2) continue;
            index_vec.push_back(chain_j);
        }
        sizePROT=index_vec.size();

        cout<<'>'<<chainID_list[chain_i]<<'\t'<<xlen<<'\t'
            <<setiosflags(ios::fixed)<<setprecision(2)
            <<100.*i/Nstruct<<"%(#"<<i<<")\t"<<"#repr="<<sizePROT<<endl;

#ifdef TMalign_Contactlib_h
        /* Contactlib filter by itself also take sometime. It is
         * only used with larger number of cluster representatives. */
        if (sizePROT>sizePROT_cutoff && xlen>xlen_cutoff)
        {
            int *thitsLocal;   // number of target hits
            int *thitsRemote;  // number of target hits
            int *dbhitsLocal;  // number of database hits
            int *dbhitsRemote; // number of database hits
            int *hitsLocal;    // number of hits between target and database
            int *hitsRemote;   // number of hits between target and database
            thitsLocal  =new int[sizePROT];
            thitsRemote =new int[sizePROT];
            dbhitsLocal =new int[sizePROT];
            dbhitsRemote=new int[sizePROT];
            hitsLocal   =new int[sizePROT];
            hitsRemote  =new int[sizePROT];
            for (j=0;j<sizePROT;j++)
                thitsLocal[j] =dbhitsLocal[j] =hitsLocal[j] =
                thitsRemote[j]=dbhitsRemote[j]=hitsRemote[j]=0;

            Contactlib_main(vindexLocal_vec[chain_i],
                vindexRemote_vec[chain_i],
                vindexLocal_vec, vindexRemote_vec, index_vec, 
                psizeLocal, psizeRemote, thitsLocal, thitsRemote,
                dbhitsLocal, dbhitsRemote, hitsLocal, hitsRemote);

            tidsizeLocal =vindexLocal_vec[chain_i].size();
            tidsizeRemote=vindexRemote_vec[chain_i].size();
            for (j=0;j<sizePROT;j++)
            {
                scoreLocal1  = 1.*thitsLocal[j]  / tidsizeLocal;
                scoreLocal2  = 1.*dbhitsLocal[j] / psizeLocal[j];
                scoreLocal   = sqrt(scoreLocal1  * scoreLocal2);

                scoreRemote1 = 1.*thitsRemote[j] / tidsizeRemote;
                scoreRemote2 = 1.*dbhitsRemote[j]/ psizeRemote[j];
                scoreRemote  = sqrt(scoreRemote1 * scoreRemote2);

                scoreCombine = scoreLocal*wLocal + scoreRemote*wRemote;
                contactlibScore_list.push_back(make_pair(scoreCombine,
                    index_vec[j]));
            }
            stable_sort(contactlibScore_list.begin(),contactlibScore_list.end(),
                greater<pair<double,int> >());

            index_vec.clear();
            for (j=0;j<sizePROT;j++)
            {
                if (j>=repr_num_cutoff) break;
                chain_j=contactlibScore_list[j].second;
                index_vec.push_back(chain_j);
                //cout<<"#"<<chain_j<<"\t"<<chainID_list[chain_j]<<"\t"
                    //<<setiosflags(ios::fixed)<<setprecision(4)
                    //<<contactlibScore_list[j].first<<endl;
            }

            contactlibScore_list.clear();
            psizeLocal.clear();
            psizeRemote.clear();
            delete[]thitsLocal;
            delete[]thitsRemote;
            delete[]dbhitsLocal;
            delete[]dbhitsRemote;
            delete[]hitsLocal;
            delete[]hitsRemote;
        }
#endif

        for (j=0;j<sizePROT;j++)
        {
            chain_j=index_vec[j];
            ylen=xyz_vec[chain_j].size();
            if (a_opt==0 && xlen<TMcut*ylen) continue;
            else if (a_opt==1 && xlen<TMcut*(xlen+ylen)/2) continue;

            NewArray(&ya, ylen, 3);
            for (r=0;r<ylen;r++)
            {
                ya[r][0]=xyz_vec[chain_j][r][0];
                ya[r][1]=xyz_vec[chain_j][r][1];
                ya[r][2]=xyz_vec[chain_j][r][2];
            }

            Lave=sqrt(xlen*ylen); // geometry average because O(L1*L2)
            bool overwrite_fast_opt=(fast_opt==true || Lave>=fast_lb);
            
            /* declare variable specific to this pair of TMalign */
            double t0[3], u0[3][3];
            double TM1, TM2;
            double TM3, TM4, TM5;     // for a_opt, u_opt, d_opt
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
            
            /* entry function for structure alignment */
            int status=TMalign_main(
                xa, ya, &seq_vec[chain_i][0], &seq_vec[chain_j][0],
                &sec_vec[chain_i][0], &sec_vec[chain_j][0],
                t0, u0, TM1, TM2, TM3, TM4, TM5,
                d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                seqM, seqxA, seqyA,
                rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                xlen, ylen, sequence, Lnorm_ass, d0_scale,
                false, I_opt, a_opt, u_opt, d_opt, overwrite_fast_opt,
                mol_vec[chain_i]+mol_vec[chain_j],TMcut);

            if (false)
            cout<<status<<'\t'<<chainID_list[chain_j]<<'\t'
                <<setiosflags(ios::fixed)<<setprecision(4)
                <<TM2<<'\t'<<TM1<<'\t'<<overwrite_fast_opt<<endl;

            seqM.clear();
            seqxA.clear();
            seqyA.clear();

            double TM=TM3; // average length
            if (a_opt==0) TM=TM1; // longer length
            else if (a_opt==-1) TM=TM2; // shorter length

            double dTM=0.1*(fast_ub-Lave)/(fast_ub-fast_lb);
            if (dTM<0) dTM=0;

            if (TM<TMcut-dTM || 
               (TM<TMcut && (fast_opt || overwrite_fast_opt==false)))
            {
                DeleteArray(&ya, ylen);
                continue;
            }

            if (TM>=TMcut+dTM || 
               (TM>=TMcut && (fast_opt || overwrite_fast_opt==false)))
            {
                clust_mem_vec[chain_i]=clust_repr_map[chain_j];
                DeleteArray(&ya, ylen);
                break;
            }

            if (Lave<fast_ub)
            {
                TMalign_main(
                    xa, ya, &seq_vec[chain_i][0], &seq_vec[chain_j][0],
                    &sec_vec[chain_i][0], &sec_vec[chain_j][0],
                    t0, u0, TM1, TM2, TM3, TM4, TM5,
                    d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                    seqM, seqxA, seqyA,
                    rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                    xlen, ylen, sequence, Lnorm_ass, d0_scale,
                    false, I_opt, a_opt, u_opt, d_opt, false,
                    mol_vec[chain_i]+mol_vec[chain_j],TMcut);
                seqM.clear();
                seqxA.clear();
                seqyA.clear();
                DeleteArray(&ya, ylen);
                
                TM=TM3;               // average length
                if (a_opt==0) TM=TM1; // longer length
                else if (a_opt==-1) TM=TM2; // shorter length
                cout<<"*\t"<<chainID_list[chain_j]<<'\t'<<TM2<<'\t'<<TM1<<endl;
                if (TM>=TMcut)
                {
                    clust_mem_vec[chain_i]=clust_repr_map[chain_j];
                    break;
                }
            }
        }
        DeleteArray(&xa, xlen);
        index_vec.clear();

        if (clust_mem_vec[chain_i]<0) // new cluster
        {
            clust_mem_vec[chain_i]=clust_repr_vec.size();
            clust_repr_map[chain_i]=clust_repr_vec.size();
            clust_repr_vec.push_back(chain_i);
        }
        else // member structures are not used further
        {
            vector<char> ().swap(seq_vec[chain_i]);
            vector<int> ().swap(sec_vec[chain_i]);
            vector<vector<float> > ().swap(xyz_vec[chain_i]);
#ifdef TMalign_Contactlib_h
            if (xlen>xlen_cutoff)
            {
                vector<vector<int> > ().swap(vindexLocal_vec[chain_i]);
                vector<vector<int> > ().swap(vindexRemote_vec[chain_i]);
            }
#endif
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

#ifdef TMalign_Contactlib_h
    if (Nstruct>sizePROT_cutoff)
    {
        vector<vector<vector<int> > >().swap(vindexLocal_vec);
        vector<vector<vector<int> > >().swap(vindexRemote_vec);
    }
#endif

    t2 = clock();
    float diff = ((float)t2 - (float)t1)/CLOCKS_PER_SEC;
    printf("Total CPU time is %5.2f seconds\n", diff);
    return 0;
}
