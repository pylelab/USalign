#include "se.h"

/* count the number of nucleic acid chains (na_chain_num) and
 * protein chains (aa_chain_num) in a complex */
int count_na_aa_chain_num(int &na_chain_num,int &aa_chain_num,
    const vector<int>&mol_vec)
{
    na_chain_num=0;
    aa_chain_num=0;
    for (size_t i=0;i<mol_vec.size();i++)
    {
        if (mol_vec[i]>0) na_chain_num++;
        else              aa_chain_num++;
    }
    return na_chain_num+aa_chain_num;
}

/* adjust chain assignment for dimer-dimer alignment 
 * return true if assignment is adjusted */
bool adjust_dimer_assignment(        
    const vector<vector<vector<double> > >&xa_vec,
    const vector<vector<vector<double> > >&ya_vec,
    const vector<int>&xlen_vec, const vector<int>&ylen_vec,
    const vector<int>&mol_vec1, const vector<int>&mol_vec2,
    int *assign1_list, int *assign2_list,
    const vector<vector<string> >&seqxA_mat,
    const vector<vector<string> >&seqyA_mat)
{
    /* check currently assigned chains */
    int i1,i2,j1,j2;
    i1=i2=j1=j2=-1;    
    int chain1_num=xa_vec.size();
    int i,j;
    for (i=0;i<chain1_num;i++)
    {
        if (assign1_list[i]>=0)
        {
            if (i1<0)
            {
                i1=i;
                j1=assign1_list[i1];
            }
            else
            {
                i2=i;
                j2=assign1_list[i2];
            }
        }
    }

    /* normalize d0 by L */
    int xlen=xlen_vec[i1]+xlen_vec[i2];
    int ylen=ylen_vec[j1]+ylen_vec[j2];
    int mol_type=mol_vec1[i1]+mol_vec1[i2]+
                 mol_vec2[j1]+mol_vec2[j2];
    double D0_MIN, d0, d0_search;
    double Lnorm=getmin(xlen,ylen);
    parameter_set4final(getmin(xlen,ylen), D0_MIN, Lnorm, d0, 
        d0_search, mol_type);

    double **xa,**ya, **xt;
    NewArray(&xa, xlen, 3);
    NewArray(&ya, ylen, 3);
    NewArray(&xt, xlen, 3);

    double RMSD = 0;
    double dd   = 0;
    double t[3];
    double u[3][3];
    size_t L_ali=0; // index of residue in aligned region
    size_t r=0;     // index of residue in full alignment

    /* total score using current assignment */
    L_ali=0;
    i=j=-1;
    for (r=0;r<seqxA_mat[i1][j1].size();r++)
    {
        i+=(seqxA_mat[i1][j1][r]!='-');
        j+=(seqyA_mat[i1][j1][r]!='-');
        if (seqxA_mat[i1][j1][r]=='-' || seqyA_mat[i1][j1][r]=='-') continue;
        xa[L_ali][0]=xa_vec[i1][i][0];
        xa[L_ali][1]=xa_vec[i1][i][1];
        xa[L_ali][2]=xa_vec[i1][i][2];
        ya[L_ali][0]=ya_vec[j1][j][0];
        ya[L_ali][1]=ya_vec[j1][j][1];
        ya[L_ali][2]=ya_vec[j1][j][2];
        L_ali++;
    }
    i=j=-1;
    for (r=0;r<seqxA_mat[i2][j2].size();r++)
    {
        i+=(seqxA_mat[i2][j2][r]!='-');
        j+=(seqyA_mat[i2][j2][r]!='-');
        if (seqxA_mat[i2][j2][r]=='-' || seqyA_mat[i2][j2][r]=='-') continue;
        xa[L_ali][0]=xa_vec[i2][i][0];
        xa[L_ali][1]=xa_vec[i2][i][1];
        xa[L_ali][2]=xa_vec[i2][i][2];
        ya[L_ali][0]=ya_vec[j2][j][0];
        ya[L_ali][1]=ya_vec[j2][j][1];
        ya[L_ali][2]=ya_vec[j2][j][2];
        L_ali++;
    }

    Kabsch(xa, ya, L_ali, 1, &RMSD, t, u);
    do_rotation(xa, xt, L_ali, t, u);

    double total_score1=0;
    for (r=0;r<L_ali;r++)
    {
        dd=dist(xt[r],ya[r]);
        total_score1+=1/(1+dd/d0*d0);
    }
    total_score1/=Lnorm;

    /* total score using reversed assignment */
    L_ali=0;
    i=j=-1;
    for (r=0;r<seqxA_mat[i1][j2].size();r++)
    {
        i+=(seqxA_mat[i1][j2][r]!='-');
        j+=(seqyA_mat[i1][j2][r]!='-');
        if (seqxA_mat[i1][j2][r]=='-' || seqyA_mat[i1][j2][r]=='-') continue;
        xa[L_ali][0]=xa_vec[i1][i][0];
        xa[L_ali][1]=xa_vec[i1][i][1];
        xa[L_ali][2]=xa_vec[i1][i][2];
        ya[L_ali][0]=ya_vec[j2][j][0];
        ya[L_ali][1]=ya_vec[j2][j][1];
        ya[L_ali][2]=ya_vec[j2][j][2];
        L_ali++;
    }
    i=j=-1;
    for (r=0;r<seqxA_mat[i2][j1].size();r++)
    {
        i+=(seqxA_mat[i2][j1][r]!='-');
        j+=(seqyA_mat[i2][j1][r]!='-');
        if (seqxA_mat[i2][j1][r]=='-' || seqyA_mat[i2][j1][r]=='-') continue;
        xa[L_ali][0]=xa_vec[i2][i][0];
        xa[L_ali][1]=xa_vec[i2][i][1];
        xa[L_ali][2]=xa_vec[i2][i][2];
        ya[L_ali][0]=ya_vec[j1][j][0];
        ya[L_ali][1]=ya_vec[j1][j][1];
        ya[L_ali][2]=ya_vec[j1][j][2];
        L_ali++;
    }

    Kabsch(xa, ya, L_ali, 1, &RMSD, t, u);
    do_rotation(xa, xt, L_ali, t, u);

    double total_score2=0;
    for (r=0;r<L_ali;r++)
    {
        dd=dist(xt[r],ya[r]);
        total_score2+=1/(1+dd/d0*d0);
    }
    total_score2/=Lnorm;

    /* swap chain assignment */
    if (total_score1<total_score2)
    {
        assign1_list[i1]=j2;
        assign1_list[i2]=j1;
        assign2_list[j1]=i2;
        assign2_list[j2]=i1;
    }

    /* clean up */
    DeleteArray(&xa, xlen);
    DeleteArray(&ya, ylen);
    DeleteArray(&xt, xlen);
    return total_score1<total_score2;
}

/* assign chain-chain correspondence */
double enhanced_greedy_search(double **TMave_mat,int *assign1_list,
    int *assign2_list, const int chain1_num, const int chain2_num)
{
    double total_score=0;
    double tmp_score=0;
    int i,j;
    int maxi=0;
    int maxj=0;

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
    //cout<<"assign1_list={";
    //for (i=0;i<chain1_num;i++) cout<<assign1_list[i]<<","; cout<<"}"<<endl;
    //cout<<"assign2_list={";
    //for (j=0;j<chain2_num;j++) cout<<assign2_list[j]<<","; cout<<"}"<<endl;

    /* iterative refinemnt */
    double delta_score;
    int *assign1_tmp=new int [chain1_num];
    int *assign2_tmp=new int [chain2_num];
    for (i=0;i<chain1_num;i++) assign1_tmp[i]=assign1_list[i];
    for (j=0;j<chain2_num;j++) assign2_tmp[j]=assign2_list[j];
    int old_i=-1;
    int old_j=-1;

    for (int iter=0;iter<getmin(chain1_num,chain2_num)*5;iter++)
    {
        delta_score=-1;
        for (i=0;i<chain1_num;i++)
        {
            old_j=assign1_list[i];
            for (j=0;j<chain2_num;j++)
            {
                // attempt to swap (i,old_j=assign1_list[i]) with (i,j)
                if (j==assign1_list[i] || TMave_mat[i][j]<=0) continue;
                old_i=assign2_list[j];

                assign1_tmp[i]=j;
                if (old_i>=0) assign1_tmp[old_i]=old_j;
                assign2_tmp[j]=i;
                if (old_j>=0) assign2_tmp[old_j]=old_i;

                delta_score=TMave_mat[i][j];
                if (old_j>=0) delta_score-=TMave_mat[i][old_j];
                if (old_i>=0) delta_score-=TMave_mat[old_i][j];
                if (old_i>=0 && old_j>=0) delta_score+=TMave_mat[old_i][old_j];

                if (delta_score>0) // successful swap
                {
                    assign1_list[i]=j;
                    if (old_i>=0) assign1_list[old_i]=old_j;
                    assign2_list[j]=i;
                    if (old_j>=0) assign2_list[old_j]=old_i;
                    total_score+=delta_score;
                    break;
                }
                else
                {
                    assign1_tmp[i]=assign1_list[i];
                    if (old_i>=0) assign1_tmp[old_i]=assign1_list[old_i];
                    assign2_tmp[j]=assign2_list[j];
                    if (old_j>=0) assign2_tmp[old_j]=assign2_list[old_j];
                }
            }
            if (delta_score>0) break;
        }
        if (delta_score<=0) break; // cannot swap any chain pair
    }

    /* clean up */
    delete[]assign1_tmp;
    delete[]assign2_tmp;
    return total_score;
}

double calculate_centroids(const vector<vector<vector<double> > >&a_vec,
    const int chain_num, double ** centroids)
{
    int L=0;
    int c,r; // index of chain and residue
    for (c=0; c<chain_num; c++)
    {
        centroids[c][0]=0;
        centroids[c][1]=0;
        centroids[c][2]=0;
        L=a_vec[c].size();
        for (r=0; r<L; r++)
        {
            centroids[c][0]+=a_vec[c][r][0];
            centroids[c][1]+=a_vec[c][r][1];
            centroids[c][2]+=a_vec[c][r][2];
        }
        centroids[c][0]/=L;
        centroids[c][1]/=L;
        centroids[c][2]/=L;
        //cout<<centroids[c][0]<<'\t'
            //<<centroids[c][1]<<'\t'
            //<<centroids[c][2]<<endl;
    }

    vector<double> d0_vec(chain_num,-1);
    int c2=0;
    double d0MM=0;
    for (c=0; c<chain_num; c++)
    {
        for (c2=0; c2<chain_num; c2++)
        {
            if (c2==c) continue;
            d0MM=sqrt(dist(centroids[c],centroids[c2]));
            if (d0_vec[c]<=0) d0_vec[c]=d0MM;
            else d0_vec[c]=getmin(d0_vec[c], d0MM);
        }
    }
    d0MM=0;
    for (c=0; c<chain_num; c++) d0MM+=d0_vec[c];
    d0MM/=chain_num;
    d0_vec.clear();
    //cout<<d0MM<<endl;
    return d0MM;
}

/* calculate MMscore of aligned chains
 * MMscore = sum(TMave_mat[i][j]) * sum(1/(1+dij^2/d0MM^2)) 
 *         / (L* getmin(chain1_num,chain2_num))
 * dij is the centroid distance between chain pair i and j
 * d0MM is scaling factor. TMave_mat[i][j] is the TM-score between
 * chain pair i and j multiple by getmin(Li*Lj) */
double calMMscore(double **TMave_mat,int *assign1_list,
    const int chain1_num, const int chain2_num, double **xcentroids,
    double **ycentroids, const double d0MM, double **r1, double **r2,
    double **xt, double t[3], double u[3][3], const int L)
{
    int Nali=0; // number of aligned chain
    int i,j;
    double MMscore=0;
    for (i=0;i<chain1_num;i++)
    {
        j=assign1_list[i];
        if (j<0) continue;

        r1[Nali][0]=xcentroids[i][0];
        r1[Nali][1]=xcentroids[i][1];
        r1[Nali][2]=xcentroids[i][2];

        r2[Nali][0]=ycentroids[j][0];
        r2[Nali][1]=ycentroids[j][1];
        r2[Nali][2]=ycentroids[j][2];

        Nali++;
        MMscore+=TMave_mat[i][j];
    }
    MMscore/=L;

    double RMSD = 0;
    double TMscore=0;
    if (Nali>=3)
    {
        /* Kabsch superposition */
        Kabsch(r1, r2, Nali, 1, &RMSD, t, u);
        do_rotation(r1, xt, Nali, t, u);

        /* calculate pseudo-TMscore */
        double dd=0;
        for (i=0;i<Nali;i++)
        {
            dd=dist(xt[i], r2[i]);
            TMscore+=1/(1+dd/(d0MM*d0MM));
        }
    }
    else if (Nali==2)
    {
        double dd=dist(r1[0],r2[0]);
        TMscore=1/(1+dd/(d0MM*d0MM));
    }
    else TMscore=1; // only one aligned chain.
    TMscore/=getmin(chain1_num,chain2_num);
    MMscore*=TMscore;
    return MMscore;
}

/* check if this is alignment of heterooligomer or homooligomer
 * return het_deg, which ranges from 0 to 1.
 * The larger the value, the more "hetero"; 
 * Tthe smaller the value, the more "homo" */
double check_heterooligomer(double **TMave_mat, const int chain1_num,
    const int chain2_num)
{
    double het_deg=0;
    double min_TM=-1;
    double max_TM=-1;
    int i,j;
    for (i=0;i<chain1_num;i++)
    {
        for (j=0;j<chain2_num;j++)
        {
            if (min_TM<0 || TMave_mat[i][j] <min_TM) min_TM=TMave_mat[i][j];
            if (max_TM<0 || TMave_mat[i][j]>=max_TM) max_TM=TMave_mat[i][j];
        }
    }
    het_deg=(max_TM-min_TM)/max_TM;
    //cout<<"min_TM="<<min_TM<<endl;
    //cout<<"max_TM="<<max_TM<<endl;
    return het_deg;
}

/* reassign chain-chain correspondence, specific for homooligomer */
double homo_refined_greedy_search(double **TMave_mat,int *assign1_list,
    int *assign2_list, const int chain1_num, const int chain2_num,
    double **xcentroids, double **ycentroids, const double d0MM,
    const int L, double **ut_mat)
{
    double MMscore_max=0;
    double MMscore=0;
    int i,j;
    int c1,c2;
    int max_i=-1; // the chain pair whose monomer u t yields highest MMscore
    int max_j=-1;

    int chain_num=getmin(chain1_num,chain2_num);
    int *assign1_tmp=new int [chain1_num];
    int *assign2_tmp=new int [chain2_num];
    double **xt;
    NewArray(&xt, chain1_num, 3);
    double t[3];
    double u[3][3];
    int ui,uj,ut_idx;
    double TMscore=0; // pseudo TM-score
    double TMsum  =0;
    double TMnow  =0;
    double TMmax  =0;
    double dd=0;

    size_t  total_pair=chain1_num*chain2_num; // total pair
    double *ut_tmc_mat=new double [total_pair]; // chain level TM-score
    vector<pair<double,int> > ut_tm_vec(total_pair,make_pair(0.0,0)); // product of both

    for (c1=0;c1<chain1_num;c1++)
    {
        for (c2=0;c2<chain2_num;c2++)
        {
            if (TMave_mat[c1][c2]<=0) continue;
            ut_idx=c1*chain2_num+c2;
            for (ui=0;ui<3;ui++)
                for (uj=0;uj<3;uj++) u[ui][uj]=ut_mat[ut_idx][ui*3+uj];
            for (uj=0;uj<3;uj++) t[uj]=ut_mat[ut_idx][9+uj];
            
            do_rotation(xcentroids, xt, chain1_num, t, u);

            for (i=0;i<chain1_num;i++) assign1_tmp[i]=-1;
            for (j=0;j<chain2_num;j++) assign2_tmp[j]=-1;


            for (i=0;i<chain1_num;i++)
            {
                for (j=0;j<chain2_num;j++)
                {
                    ut_idx=i*chain2_num+j;
                    ut_tmc_mat[ut_idx]=0;
                    ut_tm_vec[ut_idx].first=-1;
                    ut_tm_vec[ut_idx].second=ut_idx;
                    if (TMave_mat[i][j]<=0) continue;
                    dd=dist(xt[i],ycentroids[j]);
                    ut_tmc_mat[ut_idx]=1/(1+dd/(d0MM*d0MM));
                    ut_tm_vec[ut_idx].first=
                        ut_tmc_mat[ut_idx]*TMave_mat[i][j];
                    //cout<<"TM["<<ut_idx<<"]="<<ut_tm_vec[ut_idx].first<<endl;
                }
            }
            //cout<<"sorting "<<total_pair<<" chain pairs"<<endl;

            /* initial assignment */
            assign1_tmp[c1]=c2;
            assign2_tmp[c2]=c1;
            TMsum=TMave_mat[c1][c2];
            TMscore=ut_tmc_mat[c1*chain2_num+c2];

            /* further assignment */
            sort(ut_tm_vec.begin(), ut_tm_vec.end()); // sort in ascending order
            for (ut_idx=total_pair-1;ut_idx>=0;ut_idx--)
            {
                j=ut_tm_vec[ut_idx].second % chain2_num;
                i=int(ut_tm_vec[ut_idx].second / chain2_num);
                if (TMave_mat[i][j]<=0) break;
                if (assign1_tmp[i]>=0 || assign2_tmp[j]>=0) continue;
                assign1_tmp[i]=j;
                assign2_tmp[j]=i;
                TMsum+=TMave_mat[i][j];
                TMscore+=ut_tmc_mat[i*chain2_num+j];
                //cout<<"ut_idx="<<ut_tm_vec[ut_idx].second
                    //<<"\ti="<<i<<"\tj="<<j<<"\ttm="<<ut_tm_vec[ut_idx].first<<endl;
            }

            /* final MMscore */
            MMscore=(TMsum/L)*(TMscore/chain_num);
            if (max_i<0 || max_j<0 || MMscore>MMscore_max)
            {
                max_i=c1;
                max_j=c2;
                MMscore_max=MMscore;
                for (i=0;i<chain1_num;i++) assign1_list[i]=assign1_tmp[i];
                for (j=0;j<chain2_num;j++) assign2_list[j]=assign2_tmp[j];
                //cout<<"TMsum/L="<<TMsum/L<<endl;
                //cout<<"TMscore/chain_num="<<TMscore/chain_num<<endl;
                //cout<<"MMscore="<<MMscore<<endl;
                //cout<<"assign1_list={";
                //for (i=0;i<chain1_num;i++) 
                    //cout<<assign1_list[i]<<","; cout<<"}"<<endl;
                //cout<<"assign2_list={";
                //for (j=0;j<chain2_num;j++)
                    //cout<<assign2_list[j]<<","; cout<<"}"<<endl;
            }
        }
    }

    /* clean up */
    delete[]assign1_tmp;
    delete[]assign2_tmp;
    delete[]ut_tmc_mat;
    ut_tm_vec.clear();
    DeleteArray(&xt, chain1_num);
    return MMscore;
}

/* reassign chain-chain correspondence, specific for heterooligomer */
double hetero_refined_greedy_search(double **TMave_mat,int *assign1_list,
    int *assign2_list, const int chain1_num, const int chain2_num,
    double **xcentroids, double **ycentroids, const double d0MM, const int L)
{
    double MMscore_old=0;
    double MMscore=0;
    int i,j;

    double **r1;
    double **r2;
    double **xt;
    int chain_num=getmin(chain1_num,chain2_num);
    NewArray(&r1, chain_num, 3);
    NewArray(&r2, chain_num, 3);
    NewArray(&xt, chain_num, 3);
    double t[3];
    double u[3][3];

    /* calculate MMscore */
    MMscore=MMscore_old=calMMscore(TMave_mat, assign1_list, chain1_num,
        chain2_num, xcentroids, ycentroids, d0MM, r1, r2, xt, t, u, L);
    //cout<<"MMscore="<<MMscore<<endl;
    //cout<<"TMave_mat="<<endl;
    //for (i=0;i<chain1_num;i++)
    //{
        //for (j=0; j<chain2_num; j++)
        //{
            //if (j<chain2_num-1) cout<<TMave_mat[i][j]<<'\t';
            //else                cout<<TMave_mat[i][j]<<endl;
        //}
    //}

    /* iteratively refine chain assignment. in each iteration, attempt
     * to swap (i,old_j=assign1_list[i]) with (i,j) */
    double delta_score=-1;
    int *assign1_tmp=new int [chain1_num];
    int *assign2_tmp=new int [chain2_num];
    for (i=0;i<chain1_num;i++) assign1_tmp[i]=assign1_list[i];
    for (j=0;j<chain2_num;j++) assign2_tmp[j]=assign2_list[j];
    int old_i=-1;
    int old_j=-1;

    //cout<<"assign1_list={";
    //for (i=0;i<chain1_num;i++) cout<<assign1_list[i]<<","; cout<<"}"<<endl;
    //cout<<"assign2_list={";
    //for (j=0;j<chain2_num;j++) cout<<assign2_list[j]<<","; cout<<"}"<<endl;

    for (int iter=0;iter<chain1_num*chain2_num;iter++)
    {
        delta_score=-1;
        for (i=0;i<chain1_num;i++)
        {
            old_j=assign1_list[i];
            for (j=0;j<chain2_num;j++)
            {
                if (j==assign1_list[i] || TMave_mat[i][j]<=0) continue;
                old_i=assign2_list[j];

                assign1_tmp[i]=j;
                if (old_i>=0) assign1_tmp[old_i]=old_j;
                assign2_tmp[j]=i;
                if (old_j>=0) assign2_tmp[old_j]=old_i;
                
                MMscore=calMMscore(TMave_mat, assign1_tmp, chain1_num,
                    chain2_num, xcentroids, ycentroids, d0MM,
                    r1, r2, xt, t, u, L);

                //cout<<"(i,j,old_i,old_j,MMscore)=("<<i<<","<<j<<","
                    //<<old_i<<","<<old_j<<","<<MMscore<<")"<<endl;

                if (MMscore>MMscore_old) // successful swap
                {
                    assign1_list[i]=j;
                    if (old_i>=0) assign1_list[old_i]=old_j;
                    assign2_list[j]=i;
                    if (old_j>=0) assign2_list[old_j]=old_i;
                    delta_score=(MMscore-MMscore_old);
                    MMscore_old=MMscore;
                    //cout<<"MMscore="<<MMscore<<endl;
                    break;
                }
                else
                {
                    assign1_tmp[i]=assign1_list[i];
                    if (old_i>=0) assign1_tmp[old_i]=assign1_list[old_i];
                    assign2_tmp[j]=assign2_list[j];
                    if (old_j>=0) assign2_tmp[old_j]=assign2_list[old_j];
                }
            }
        }
        //cout<<"iter="<<iter<<endl;
        //cout<<"assign1_list={";
        //for (i=0;i<chain1_num;i++) cout<<assign1_list[i]<<","; cout<<"}"<<endl;
        //cout<<"assign2_list={";
        //for (j=0;j<chain2_num;j++) cout<<assign2_list[j]<<","; cout<<"}"<<endl;
        if (delta_score<=0) break; // cannot swap any chain pair
    }
    MMscore=MMscore_old;
    //cout<<"MMscore="<<MMscore<<endl;

    /* clean up */
    delete[]assign1_tmp;
    delete[]assign2_tmp;
    DeleteArray(&r1, chain_num);
    DeleteArray(&r2, chain_num);
    DeleteArray(&xt, chain_num);
    return MMscore;
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
    const int mirror_opt, const int het_opt, int &len_aa, int &len_na,  
    const int o_opt, vector<string>&resi_vec)
{
    size_t i;
    int chain_i,r;
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
    //vector<string> resi_vec;
    int read_resi=0;
    if (o_opt) read_resi=2;

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
            len = read_PDB(PDB_lines[chain_i], xa, seq, resi_vec, read_resi);
            if (mirror_opt) for (r=0;r<len;r++) xa[r][2]=-xa[r][2];
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

void MMalign_final(
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
    vector<vector<string> >&seqxA_mat, vector<vector<string> >&seqM_mat,
    vector<vector<string> >&seqyA_mat, int *assign1_list, int *assign2_list,
    vector<string>&sequence, const double d0_scale, const bool m_opt,
    const int o_opt, const int outfmt_opt, const int ter_opt,
    const int split_opt, const bool a_opt, const bool d_opt,
    const bool fast_opt, const bool full_opt, const int mirror_opt,
    const vector<string>&resi_vec1, const vector<string>&resi_vec2)
{
    int i,j;
    int xlen=0;
    int ylen=0;
    for (i=0;i<chain1_num;i++) xlen+=xlen_vec[i];
    for (j=0;j<chain2_num;j++) ylen+=ylen_vec[j];
    if (xlen<=3 || ylen<=3) return;

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

    /* prepare unaligned region */
    for (i=0;i<chain1_num;i++)
    {
        if (assign1_list[i]>=0) continue;
        chainID1+=chainID_list1[i];
        chainID2+=':';
        string s(seqx_vec[i].begin(),seqx_vec[i].end());
        sequence[0]+=s.substr(0,xlen_vec[i])+'*';
        sequence[1]+=string(xlen_vec[i],'-')+'*';
        s.clear();
        sequence[2]+=string(xlen_vec[i],' ')+'*';
    }
    for (j=0;j<chain2_num;j++)
    {
        if (assign2_list[j]>=0) continue;
        chainID1+=':';
        chainID2+=chainID_list2[j];
        string s(seqy_vec[j].begin(),seqy_vec[j].end());
        sequence[0]+=string(ylen_vec[j],'-')+'*';
        sequence[1]+=s.substr(0,ylen_vec[j])+'*';
        s.clear();
        sequence[2]+=string(ylen_vec[j],' ')+'*';
    }

    /* print alignment */
    output_results(xname, yname, chainID1.c_str(), chainID2.c_str(),
        xlen, ylen, t0, u0, TM1, TM2, TM3, TM4, TM5, rmsd0, d0_out,
        sequence[2].c_str(), sequence[0].c_str(), sequence[1].c_str(),
        Liden, n_ali8, L_ali, TM_ali, rmsd_ali,
        TM_0, d0_0, d0A, d0B, 0, d0_scale, d0a, d0u, 
        (m_opt?fname_matrix:"").c_str(), outfmt_opt, ter_opt, true,
        split_opt, o_opt, fname_super,
        false, a_opt, false, d_opt, mirror_opt, resi_vec1, resi_vec2);

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

    if (!full_opt) return;

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
            "", outfmt_opt, ter_opt, false, split_opt, 0,
            "", false, a_opt, false, d_opt, 0, resi_vec1, resi_vec2);

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
    return;
}
