#include "se.h"

/* assign chain-chain correspondence */
double enhanced_greedy_search(double **TMave_mat,int *assign1_list,
    int *assign2_list, const int chain1_num, const int chain2_num)
{
    double total_score=0;
    double tmp_score=0;
    int i,j;
    int maxi,maxj;

    /* initialize parameters */
    for (i=0;i<chain1_num;i++) assign1_list[i]=-1;
    for (j=0;j<chain2_num;j++) assign2_list[j]=-1;

    /* greedy assignment: in each iteration, the highest chain pair is
     * assigned, until no assignable chain is left */
    while(1)
    {
        tmp_score=-1;
        for (i=0;i<chain1_num;i++)
        {
            if (assign1_list[i]>=0) continue;
            for (j=0;j<chain2_num;j++)
            {
                if (assign2_list[j]>=0 || TMave_mat[i][j]<=0) continue;
                if (TMave_mat[i][j]>tmp_score) 
                {
                    maxi=i;
                    maxj=j;
                    tmp_score=TMave_mat[i][j];
                }
            }
        }
        if (tmp_score<=0) break; // error: no assignable chain
        assign1_list[maxi]=maxj;
        assign2_list[maxj]=maxi;
        total_score+=tmp_score;
    }
    if (total_score<=0) return total_score; // error: no assignable chain

    /* iterative refinemnt */
    double delta_score;
    for (int iter=0;iter<getmin(chain1_num,chain2_num)*5;iter++)
    {
        delta_score=0;
        for (i=0;i<chain1_num;i++)
        {
            if (assign1_list[i]<0) continue;
            for (j=0;j<chain2_num;j++)
            {
                // attempt to swap (i,assign1_list[i]) with
                // (assign2_list[j],j)
                if (j==assign1_list[i] || TMave_mat[i][j]<=0) continue;

                delta_score=TMave_mat[i][j]-TMave_mat[i][assign1_list[i]];
                if (assign2_list[j]>=0) delta_score+=
                    TMave_mat[assign2_list[j]][assign1_list[i]]
                   -TMave_mat[assign2_list[j]][j];

                if (delta_score>0) // successful swap
                {
                    assign1_list[i]=j;
                    assign2_list[j]=i;
                    total_score+=delta_score;
                    break;
                }
            }
            if (delta_score>0) break;
        }
        if (delta_score<=0) break; // cannot swap any chain pair
    }
    return total_score;
}

void copy_chain_data(const vector<vector<double> >&a_vec_i,
    const vector<char>&seq_vec_i,const vector<char>&sec_vec_i,
    const int len,double **a,char *seq,char *sec)
{
    int r;
    for (r=0;r<len;r++)
    {
        a[r][0]=a_vec_i[r][0];
        a[r][1]=a_vec_i[r][1];
        a[r][2]=a_vec_i[r][2];
        seq[r]=seq_vec_i[r];
        sec[r]=sec_vec_i[r];
    }
    seq[len]=0;
    sec[len]=0;
}

void parse_chain_list(const vector<string>&chain_list,
    vector<vector<vector<double> > >&a_vec, vector<vector<char> >&seq_vec,
    vector<vector<char> >&sec_vec, vector<int>&mol_vec, vector<int>&len_vec,
    vector<string>&chainID_list, const int ter_opt, const int split_opt,
    const string mol_opt, const int infmt_opt, const string atom_opt,
    const int het_opt, int &len_aa, int &len_na)
{
    int i,chain_i,r;
    string name;
    int chainnum;
    double **xa;
    int len;
    char *seq,*sec;

    vector<vector<string> >PDB_lines;
    vector<double> tmp_atom_array(3,0);
    vector<vector<double> > tmp_chain_array;
    vector<char>tmp_seq_array;
    vector<char>tmp_sec_array;
    vector<string> resi_vec;

    for (i=0;i<chain_list.size();i++)
    {
        name=chain_list[i];
        chainnum=get_PDB_lines(name, PDB_lines, chainID_list,
            mol_vec, ter_opt, infmt_opt, atom_opt, split_opt, het_opt);
        if (!chainnum)
        {
            cerr<<"Warning! Cannot parse file: "<<name
                <<". Chain number 0."<<endl;
            continue;
        }
        for (chain_i=0;chain_i<chainnum;chain_i++)
        {
            len=PDB_lines[chain_i].size();
            if (!len)
            {
                cerr<<"Warning! Cannot parse file: "<<name
                    <<". Chain length 0."<<endl;
                continue;
            }
            else if (len<3)
            {
                cerr<<"Sequence is too short <3!: "<<name<<endl;
                continue;
            }
            NewArray(&xa, len, 3);
            seq = new char[len + 1];
            sec = new char[len + 1];
            len = read_PDB(PDB_lines[chain_i], xa, seq, resi_vec, 0);
            if (mol_vec[chain_i]>0 || mol_opt=="RNA")
                make_sec(seq, xa, len, sec,atom_opt);
            else make_sec(xa, len, sec); // secondary structure assignment
            
            /* store in vector */
            tmp_chain_array.assign(len,tmp_atom_array);
            vector<char>tmp_seq_array(len+1,0);
            vector<char>tmp_sec_array(len+1,0);
            for (r=0;r<len;r++)
            {
                tmp_chain_array[r][0]=xa[r][0];
                tmp_chain_array[r][1]=xa[r][1];
                tmp_chain_array[r][2]=xa[r][2];
                tmp_seq_array[r]=seq[r];
                tmp_sec_array[r]=sec[r];
            }
            a_vec.push_back(tmp_chain_array);
            seq_vec.push_back(tmp_seq_array);
            sec_vec.push_back(tmp_sec_array);
            len_vec.push_back(len);

            /* clean up */
            tmp_chain_array.clear();
            tmp_seq_array.clear();
            tmp_sec_array.clear();
            PDB_lines[chain_i].clear();
            DeleteArray(&xa, len);
            delete [] seq;
            delete [] sec;
        } // chain_i
        name.clear();
        PDB_lines.clear();
        mol_vec.clear();
    } // i
    tmp_atom_array.clear();

    if (mol_opt=="RNA") mol_vec.assign(a_vec.size(),1);
    else if (mol_opt=="protein") mol_vec.assign(a_vec.size(),-1);
    else
    {
        mol_vec.assign(a_vec.size(),0);
        for (i=0;i<a_vec.size();i++)
        {
            for (r=0;r<len_vec[i];r++)
            {
                if (seq_vec[i][r]>='a' && seq_vec[i][r]<='z') mol_vec[i]++;
                else mol_vec[i]--;
            }
        }
    }

    len_aa=0;
    len_na=0;
    for (i=0;i<a_vec.size();i++)
    {
        if (mol_vec[i]>0) len_na+=len_vec[i];
        else              len_aa+=len_vec[i];
    }
}

int copy_chain_pair_data(
    const vector<vector<vector<double> > >&xa_vec,
    const vector<vector<vector<double> > >&ya_vec,
    const vector<vector<char> >&seqx_vec, const vector<vector<char> >&seqy_vec,
    const vector<vector<char> >&secx_vec, const vector<vector<char> >&secy_vec,
    const vector<int> &mol_vec1, const vector<int> &mol_vec2,
    const vector<int> &xlen_vec, const vector<int> &ylen_vec,
    double **xa, double **ya, char *seqx, char *seqy, char *secx, char *secy,
    int chain1_num, int chain2_num,
    vector<vector<string> >&seqxA_mat, vector<vector<string> >&seqyA_mat,
    int *assign1_list, int *assign2_list, vector<string>&sequence)
{
    int i,j,r;
    sequence.clear();
    sequence.push_back("");
    sequence.push_back("");
    int mol_type=0;
    int xlen=0;
    int ylen=0;
    for (i=0;i<chain1_num;i++)
    {
        j=assign1_list[i];
        if (j<0) continue;
        for (r=0;r<xlen_vec[i];r++)
        {
            seqx[xlen]=seqx_vec[i][r];
            secx[xlen]=secx_vec[i][r];
            xa[xlen][0]= xa_vec[i][r][0];
            xa[xlen][1]= xa_vec[i][r][1];
            xa[xlen][2]= xa_vec[i][r][2];
            xlen++;
        }
        sequence[0]+=seqxA_mat[i][j];
        for (r=0;r<ylen_vec[j];r++)
        {
            seqy[ylen]=seqy_vec[j][r];
            secy[ylen]=secy_vec[j][r];
            ya[ylen][0]= ya_vec[j][r][0];
            ya[ylen][1]= ya_vec[j][r][1];
            ya[ylen][2]= ya_vec[j][r][2];
            ylen++;
        }
        sequence[1]+=seqyA_mat[i][j];
        mol_type+=mol_vec1[i]+mol_vec2[j];
    }
    seqx[xlen]=0;
    secx[xlen]=0;
    seqy[ylen]=0;
    secy[ylen]=0;
    return mol_type;
}

double MMalign_search(
    const vector<vector<vector<double> > >&xa_vec,
    const vector<vector<vector<double> > >&ya_vec,
    const vector<vector<char> >&seqx_vec, const vector<vector<char> >&seqy_vec,
    const vector<vector<char> >&secx_vec, const vector<vector<char> >&secy_vec,
    const vector<int> &mol_vec1, const vector<int> &mol_vec2,
    const vector<int> &xlen_vec, const vector<int> &ylen_vec,
    double **xa, double **ya, char *seqx, char *seqy, char *secx, char *secy,
    int len_aa, int len_na, int chain1_num, int chain2_num,
    double **TM1_mat, double **TM2_mat, double **TMave_mat,
    vector<vector<string> >&seqxA_mat, vector<vector<string> >&seqyA_mat,
    int *assign1_list, int *assign2_list, vector<string>&sequence,
    double d0_scale, bool fast_opt)
{
    double total_score=0;
    int i,j;
    int xlen=0;
    int ylen=0;
    for (i=0;i<chain1_num;i++)
    {
        if (assign1_list[i]<0) continue;
        xlen+=xlen_vec[i];
        ylen+=ylen_vec[assign1_list[i]];
    }
    if (xlen<=3 || ylen<=3) return total_score;

    seqx = new char[xlen+1];
    secx = new char[xlen+1];
    NewArray(&xa, xlen, 3);
    seqy = new char[ylen+1];
    secy = new char[ylen+1];
    NewArray(&ya, ylen, 3);

    int mol_type=copy_chain_pair_data(xa_vec, ya_vec, seqx_vec, seqy_vec,
        secx_vec, secy_vec, mol_vec1, mol_vec2, xlen_vec, ylen_vec,
        xa, ya, seqx, seqy, secx, secy, chain1_num, chain2_num,
        seqxA_mat, seqyA_mat, assign1_list, assign2_list, sequence);

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

    double Lnorm_ass=len_aa+len_na;

    /* entry function for structure alignment */
    TMalign_main(xa, ya, seqx, seqy, secx, secy,
        t0, u0, TM1, TM2, TM3, TM4, TM5,
        d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out, seqM, seqxA, seqyA,
        rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
        xlen, ylen, sequence, Lnorm_ass, d0_scale,
        3, false, true, false, fast_opt, mol_type, -1);

    /* clean up */
    delete [] seqx;
    delete [] seqy;
    delete [] secx;
    delete [] secy;
    DeleteArray(&xa,xlen);
    DeleteArray(&ya,ylen);

    /* re-compute chain level alignment */
    for (i=0;i<chain1_num;i++)
    {
        xlen=xlen_vec[i];
        if (xlen<3)
        {
            for (j=0;j<chain2_num;j++)
                TM1_mat[i][j]=TM2_mat[i][j]=TMave_mat[i][j]=-1;
            continue;
        }
        seqx = new char[xlen+1];
        secx = new char[xlen+1];
        NewArray(&xa, xlen, 3);
        copy_chain_data(xa_vec[i],seqx_vec[i],secx_vec[i],
            xlen,xa,seqx,secx);

        double **xt;
        NewArray(&xt, xlen, 3);
        do_rotation(xa, xt, xlen, t0, u0);

        for (j=0;j<chain2_num;j++)
        {
            if (mol_vec1[i]*mol_vec2[j]<0) //no protein-RNA alignment
            {
                TM1_mat[i][j]=TM2_mat[i][j]=TMave_mat[i][j]=-1;
                continue;
            }

            ylen=ylen_vec[j];
            if (ylen<3)
            {
                TM1_mat[i][j]=TM2_mat[i][j]=TMave_mat[i][j]=-1;
                continue;
            }
            seqy = new char[ylen+1];
            secy = new char[ylen+1];
            NewArray(&ya, ylen, 3);
            copy_chain_data(ya_vec[j],seqy_vec[j],secy_vec[j],
                ylen,ya,seqy,secy);

            /* declare variable specific to this pair of TMalign */
            d0_out=5.0;
            seqM.clear();
            seqxA.clear();
            seqyA.clear();
            rmsd0 = 0.0;
            Liden=0;
            int *invmap = new int[ylen+1];

            double Lnorm_ass=len_aa;
            if (mol_vec1[i]+mol_vec2[j]>0) Lnorm_ass=len_na;

            /* entry function for structure alignment */
            se_main(xt, ya, seqx, seqy, TM1, TM2, TM3, TM4, TM5,
                d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out, seqM, seqxA, seqyA,
                rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                xlen, ylen, sequence, Lnorm_ass, d0_scale,
                0, false, true, false,
                mol_vec1[i]+mol_vec2[j], 1, invmap);

            /* print result */
            TM1_mat[i][j]=TM2; // normalized by chain1
            TM2_mat[i][j]=TM1; // normalized by chain2
            seqxA_mat[i][j]=seqxA;
            seqyA_mat[i][j]=seqyA;

            TMave_mat[i][j]=TM4*Lnorm_ass;

            /* clean up */
            seqM.clear();
            seqxA.clear();
            seqyA.clear();

            delete[]seqy;
            delete[]secy;
            DeleteArray(&ya,ylen);
        }
        delete[]seqx;
        delete[]secx;
        DeleteArray(&xa,xlen);
        DeleteArray(&xt,xlen);
    }
    return total_score;
}

double MMalign_final(
    const string xname, const string yname,
    const vector<string> chainID_list1, const vector<string> chainID_list2,
    string fname_super, string fname_lign, string fname_matrix,
    const vector<vector<vector<double> > >&xa_vec,
    const vector<vector<vector<double> > >&ya_vec,
    const vector<vector<char> >&seqx_vec, const vector<vector<char> >&seqy_vec,
    const vector<vector<char> >&secx_vec, const vector<vector<char> >&secy_vec,
    const vector<int> &mol_vec1, const vector<int> &mol_vec2,
    const vector<int> &xlen_vec, const vector<int> &ylen_vec,
    double **xa, double **ya, char *seqx, char *seqy, char *secx, char *secy,
    int len_aa, int len_na, int chain1_num, int chain2_num,
    double **TM1_mat, double **TM2_mat, double **TMave_mat,
    vector<vector<string> >&seqxA_mat, 
    vector<vector<string> >&seqM_mat,
    vector<vector<string> >&seqyA_mat,
    int *assign1_list, int *assign2_list, vector<string>&sequence,
    double d0_scale, bool m_opt, bool o_opt, int outfmt_opt, int ter_opt,
    bool a_opt, bool d_opt, bool fast_opt, bool full_opt)
{
    double total_score=0;
    int i,j;
    int xlen=0;
    int ylen=0;
    for (i=0;i<chain1_num;i++)
    {
        if (assign1_list[i]<0) continue;
        xlen+=xlen_vec[i];
        ylen+=ylen_vec[assign1_list[i]];
    }
    if (xlen<=3 || ylen<=3) return total_score;

    seqx = new char[xlen+1];
    secx = new char[xlen+1];
    NewArray(&xa, xlen, 3);
    seqy = new char[ylen+1];
    secy = new char[ylen+1];
    NewArray(&ya, ylen, 3);

    int mol_type=copy_chain_pair_data(xa_vec, ya_vec, seqx_vec, seqy_vec,
        secx_vec, secy_vec, mol_vec1, mol_vec2, xlen_vec, ylen_vec,
        xa, ya, seqx, seqy, secx, secy, chain1_num, chain2_num,
        seqxA_mat, seqyA_mat, assign1_list, assign2_list, sequence);

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

    double Lnorm_ass=len_aa+len_na;

    /* entry function for structure alignment */
    TMalign_main(xa, ya, seqx, seqy, secx, secy,
        t0, u0, TM1, TM2, TM3, TM4, TM5,
        d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out, seqM, seqxA, seqyA,
        rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
        xlen, ylen, sequence, Lnorm_ass, d0_scale,
        3, a_opt, false, d_opt, fast_opt, mol_type, -1);

    /* prepare full complex alignment */
    string chainID1="";
    string chainID2="";
    sequence.clear();
    sequence.push_back(""); // seqxA
    sequence.push_back(""); // seqyA
    sequence.push_back(""); // seqM
    int aln_start=0;
    int aln_end=0;
    for (i=0;i<chain1_num;i++)
    {
        j=assign1_list[i];
        if (j<0) continue;
        chainID1+=chainID_list1[i];
        chainID2+=chainID_list2[j];
        sequence[0]+=seqxA_mat[i][j]+'*';
        sequence[1]+=seqyA_mat[i][j]+'*';

        aln_end+=seqxA_mat[i][j].size();
        seqM_mat[i][j]=seqM.substr(aln_start,aln_end-aln_start);
        sequence[2]+=seqM_mat[i][j]+'*';
        aln_start=aln_end;
    }

    output_results(xname, yname, chainID1.c_str(), chainID2.c_str(),
        xlen, ylen, t0, u0, TM1, TM2, TM3, TM4, TM5, rmsd0, d0_out,
        sequence[2].c_str(), sequence[0].c_str(), sequence[1].c_str(),
        Liden, n_ali8, L_ali, TM_ali, rmsd_ali,
        TM_0, d0_0, d0A, d0B, 0, d0_scale, d0a, d0u, 
        (m_opt?fname_matrix:"").c_str(), outfmt_opt, ter_opt, 
        (o_opt?fname_super:"").c_str(),
        false, a_opt, false, d_opt, false);

    /* clean up */
    seqM.clear();
    seqxA.clear();
    seqyA.clear();
    delete [] seqx;
    delete [] seqy;
    delete [] secx;
    delete [] secy;
    DeleteArray(&xa,xlen);
    DeleteArray(&ya,ylen);
    sequence[0].clear();
    sequence[1].clear();
    sequence[2].clear();

    if (!full_opt) return total_score;

    cout<<"# End of alignment for full complex. The following blocks list alignments for individual chains."<<endl;

    /* re-compute chain level alignment */
    for (i=0;i<chain1_num;i++)
    {
        j=assign1_list[i];
        if (j<0) continue;
        xlen=xlen_vec[i];
        seqx = new char[xlen+1];
        secx = new char[xlen+1];
        NewArray(&xa, xlen, 3);
        copy_chain_data(xa_vec[i],seqx_vec[i],secx_vec[i],
            xlen,xa,seqx,secx);

        double **xt;
        NewArray(&xt, xlen, 3);
        do_rotation(xa, xt, xlen, t0, u0);

        ylen=ylen_vec[j];
        if (ylen<3)
        {
            TM1_mat[i][j]=TM2_mat[i][j]=TMave_mat[i][j]=-1;
            continue;
        }
        seqy = new char[ylen+1];
        secy = new char[ylen+1];
        NewArray(&ya, ylen, 3);
        copy_chain_data(ya_vec[j],seqy_vec[j],secy_vec[j],
            ylen,ya,seqy,secy);

        /* declare variable specific to this pair of TMalign */
        d0_out=5.0;
        rmsd0 = 0.0;
        Liden=0;
        int *invmap = new int[ylen+1];
        seqM="";
        seqxA="";
        seqyA="";
        double Lnorm_ass=len_aa;
        if (mol_vec1[i]+mol_vec2[j]>0) Lnorm_ass=len_na;
        sequence[0]=seqxA_mat[i][j];
        sequence[1]=seqyA_mat[i][j];

        /* entry function for structure alignment */
        se_main(xt, ya, seqx, seqy, TM1, TM2, TM3, TM4, TM5,
            d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out, seqM, seqxA, seqyA,
            rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
            xlen, ylen, sequence, Lnorm_ass, d0_scale,
            1, a_opt, true, d_opt, mol_vec1[i]+mol_vec2[j], 1, invmap);

        //TM2=TM4*Lnorm_ass/xlen;
        //TM1=TM4*Lnorm_ass/ylen;
        //d0A=d0u;
        //d0B=d0u;

        /* print result */
        output_results(xname, yname,
            chainID_list1[i].c_str(), chainID_list2[j].c_str(),
            xlen, ylen, t0, u0, TM1, TM2, TM3, TM4, TM5, rmsd0, d0_out,
            seqM_mat[i][j].c_str(), seqxA_mat[i][j].c_str(),
            seqyA_mat[i][j].c_str(), Liden, n_ali8, L_ali, TM_ali, rmsd_ali,
            TM_0, d0_0, d0A, d0B, Lnorm_ass, d0_scale, d0a, d0u, 
            "", outfmt_opt, ter_opt, "", false, a_opt, false, d_opt, 0);

        /* clean up */
        seqxA.clear();
        seqM.clear();
        seqyA.clear();
        sequence[0].clear();
        sequence[1].clear();
        delete[]seqy;
        delete[]secy;
        DeleteArray(&ya,ylen);
        delete[]seqx;
        delete[]secx;
        DeleteArray(&xa,xlen);
        DeleteArray(&xt,xlen);
    }
    sequence.clear();
    return total_score;
}
