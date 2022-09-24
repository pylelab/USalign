#include <cfloat>
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

/* count how many chains are paired */
int count_assign_pair(int *assign1_list,const int chain1_num)
{
    int pair_num=0;
    int i;
    for (i=0;i<chain1_num;i++) pair_num+=(assign1_list[i]>=0);
    return pair_num;
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

/* clear chains with L<3 */
void clear_full_PDB_lines(vector<vector<string> > PDB_lines,const string atom_opt)
{
    int chain_i;
    int Lch;
    int a;
    bool select_atom;
    string line;
    for (chain_i=0;chain_i<PDB_lines.size();chain_i++)
    {
        Lch=0;
        for (a=0;a<PDB_lines[chain_i].size();a++)
        {
            line=PDB_lines[chain_i][a];
            if (atom_opt=="auto")
            {
                if (line[17]==' ' && (line[18]=='D'||line[18]==' '))
                     select_atom=(line.compare(12,4," C3'")==0);
                else select_atom=(line.compare(12,4," CA ")==0);
            }
            else     select_atom=(line.compare(12,4,atom_opt)==0);
            Lch+=select_atom;
        }
        if (Lch<3)
        {
            for (a=0;a<PDB_lines[chain_i].size();a++)
                PDB_lines[chain_i][a].clear();
            PDB_lines[chain_i].clear();
        }
    }
    line.clear();
}

size_t get_full_PDB_lines(const string filename,
    vector<vector<string> >&PDB_lines, const int ter_opt,
    const int infmt_opt, const int split_opt, const int het_opt)
{
    size_t i=0; // resi i.e. atom index
    string line;
    char chainID=0;
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
        fin.open(filename.c_str());

    if (infmt_opt==0||infmt_opt==-1) // PDB format
    {
        while (compress_type?fin_gz.good():fin.good())
        {
            if (compress_type) getline(fin_gz, line);
            else               getline(fin, line);
            if (infmt_opt==-1 && line.compare(0,5,"loop_")==0) // PDBx/mmCIF
                return get_full_PDB_lines(filename,PDB_lines,
                    ter_opt, 3, split_opt,het_opt);
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
                if (!chainID)
                {
                    chainID=line[21];
                    PDB_lines.push_back(tmp_str_vec);
                }
                else if (ter_opt>=2 && chainID!=line[21]) break;
                if (split_opt==2 && chainID!=line[21])
                {
                    chainID=line[21];
                    PDB_lines.push_back(tmp_str_vec);
                } 
               
                PDB_lines.back().push_back(line);
                i++;
            }
        }
    }
    else if (infmt_opt==1) // SPICKER format
    {
        size_t L=0;
        float x,y,z;
        stringstream i8_stream;
        while (compress_type?fin_gz.good():fin.good())
        {
            if (compress_type) fin_gz>>L>>x>>y>>z;
            else               fin   >>L>>x>>y>>z;
            if (compress_type) getline(fin_gz, line);
            else               getline(fin, line);
            if (!(compress_type?fin_gz.good():fin.good())) break;
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
        size_t L=0;
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
            PDB_lines.push_back(tmp_str_vec);
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
        string resi="";
        string model_index=""; // the same as model_idx but type is string
        stringstream i8_stream;
        while (compress_type?fin_gz.good():fin.good())
        {
            if (compress_type) getline(fin_gz, line);
            else               getline(fin, line);
            if (line.size()==0) continue;
            if (loop_) loop_ = (line.size()>=2)?(line.compare(0,2,"# ")):(line.compare(0,1,"#"));
            if (!loop_)
            {
                if (line.compare(0,5,"loop_")) continue;
                while(1)
                {
                    if (compress_type)
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
                    if (compress_type) getline(fin_gz, line);
                    else               getline(fin, line);
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
                    prev_asym_id=asym_id;
                }
            }

            if (prev_asym_id!=asym_id)
            {
                if (prev_asym_id!="" && ter_opt>=2) break;
                if (split_opt>=2) PDB_lines.push_back(tmp_str_vec);
            }
            if (prev_asym_id!=asym_id) prev_asym_id=asym_id;

            if (_atom_site.count("auth_seq_id"))
                 resi=line_vec[_atom_site["auth_seq_id"]];
            else resi=line_vec[_atom_site["label_seq_id"]];
            if (_atom_site.count("pdbx_PDB_ins_code") && 
                line_vec[_atom_site["pdbx_PDB_ins_code"]]!="?")
                resi+=line_vec[_atom_site["pdbx_PDB_ins_code"]][0];
            else resi+=" ";

            i++;
            i8_stream<<"ATOM  "
                <<setw(5)<<i<<" "<<atom<<" "<<AA<<setw(2)<<asym_id.substr(0,2)
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

    if (compress_type) fin_gz.close();
    else               fin.close();
    line.clear();
    return PDB_lines.size();
}

void output_dock(const vector<string>&chain_list, const int ter_opt,
    const int split_opt, const int infmt_opt, const string atom_opt,
    const int mirror_opt, double **ut_mat, const string&fname_super)
{
    size_t i;
    int chain_i,a;
    string name;
    int chainnum;
    double x[3];  // before transform
    double x1[3]; // after transform
    string line;
    vector<vector<string> >PDB_lines;
    int m=0;
    double t[3];
    double u[3][3];
    int ui,uj;
    stringstream buf;
    string filename;
    int het_opt=1;
    for (i=0;i<chain_list.size();i++)
    {
        name=chain_list[i];
        chainnum=get_full_PDB_lines(name, PDB_lines,
            ter_opt, infmt_opt, split_opt, het_opt);
        if (!chainnum) continue;
        clear_full_PDB_lines(PDB_lines, atom_opt); // clear chains with <3 residue
        for (chain_i=0;chain_i<chainnum;chain_i++)
        {
            if (PDB_lines[chain_i].size()<3) continue;
            buf<<fname_super<<'.'<<m<<".pdb";
            filename=buf.str();
            buf.str(string());
            for (ui=0;ui<3;ui++) for (uj=0;uj<3;uj++) u[ui][uj]=ut_mat[m][ui*3+uj];
            for (uj=0;uj<3;uj++) t[uj]=ut_mat[m][9+uj];
            for (a=0;a<PDB_lines[chain_i].size();a++)
            {
                line=PDB_lines[chain_i][a];
                x[0]=atof(line.substr(30,8).c_str());
                x[1]=atof(line.substr(38,8).c_str());
                x[2]=atof(line.substr(46,8).c_str());
                if (mirror_opt) x[2]=-x[2];
                transform(t, u, x, x1);
                buf<<line.substr(0,30)<<setiosflags(ios::fixed)
                   <<setprecision(3)
                   <<setw(8)<<x1[0]<<setw(8)<<x1[1]<<setw(8)<<x1[2]
                   <<line.substr(54)<<'\n';
            }
            buf<<"TER"<<endl;
            ofstream fp;
            fp.open(filename.c_str());
            fp<<buf.str();
            fp.close();
            buf.str(string());
            PDB_lines[chain_i].clear();
            m++;
        } // chain_i
        name.clear();
        PDB_lines.clear();
    } // i
    vector<vector<string> >().swap(PDB_lines);
    line.clear();
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
    int read_resi=2;

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
    for (i=0;i<sequence.size();i++) sequence[i].clear();
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
    int len_aa, int len_na, int chain1_num, int chain2_num, double **TMave_mat,
    vector<vector<string> >&seqxA_mat, vector<vector<string> >&seqyA_mat,
    int *assign1_list, int *assign2_list, vector<string>&sequence,
    double d0_scale, bool fast_opt, const int i_opt=3)
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
        i_opt, false, true, false, fast_opt, mol_type, -1);

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
            for (j=0;j<chain2_num;j++) TMave_mat[i][j]=-1;
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
                TMave_mat[i][j]=-1;
                continue;
            }

            ylen=ylen_vec[j];
            if (ylen<3)
            {
                TMave_mat[i][j]=-1;
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
                0, false, 2, false, mol_vec1[i]+mol_vec2[j], 1, invmap);

            /* print result */
            seqxA_mat[i][j]=seqxA;
            seqyA_mat[i][j]=seqyA;

            TMave_mat[i][j]=TM4*Lnorm_ass;
            if (assign1_list[i]==j) total_score+=TMave_mat[i][j];

            /* clean up */
            seqM.clear();
            seqxA.clear();
            seqyA.clear();

            delete[]seqy;
            delete[]secy;
            DeleteArray(&ya,ylen);
            delete[]invmap;
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
    double **TMave_mat,
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
            TMave_mat[i][j]=-1;
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
            1, a_opt, 2, d_opt, mol_vec1[i]+mol_vec2[j], 1, invmap);

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
        delete[]invmap;
    }
    sequence.clear();
    return;
}

void copy_chain_assign_data(int chain1_num, int chain2_num, 
    vector<string> &sequence,
    vector<vector<string> >&seqxA_mat, vector<vector<string> >&seqyA_mat,
    int *assign1_list, int *assign2_list, double **TMave_mat,
    vector<vector<string> >&seqxA_tmp, vector<vector<string> >&seqyA_tmp,
    int *assign1_tmp,  int *assign2_tmp,  double **TMave_tmp)
{
    int i,j;
    for (i=0;i<sequence.size();i++) sequence[i].clear();
    sequence.clear();
    sequence.push_back("");
    sequence.push_back("");
    for (i=0;i<chain1_num;i++) assign1_tmp[i]=assign1_list[i];
    for (i=0;i<chain2_num;i++) assign2_tmp[i]=assign2_list[i];
    for (i=0;i<chain1_num;i++)
    {
        for (j=0;j<chain2_num;j++)
        {
            seqxA_tmp[i][j]=seqxA_mat[i][j];
            seqyA_tmp[i][j]=seqyA_mat[i][j];
            TMave_tmp[i][j]=TMave_mat[i][j];
            if (assign1_list[i]==j)
            {
                sequence[0]+=seqxA_mat[i][j];
                sequence[1]+=seqyA_mat[i][j];
            }
        }
    }
    return;
}

void MMalign_iter(double & max_total_score, const int max_iter,
    const vector<vector<vector<double> > >&xa_vec,
    const vector<vector<vector<double> > >&ya_vec,
    const vector<vector<char> >&seqx_vec, const vector<vector<char> >&seqy_vec,
    const vector<vector<char> >&secx_vec, const vector<vector<char> >&secy_vec,
    const vector<int> &mol_vec1, const vector<int> &mol_vec2,
    const vector<int> &xlen_vec, const vector<int> &ylen_vec,
    double **xa, double **ya, char *seqx, char *seqy, char *secx, char *secy,
    int len_aa, int len_na, int chain1_num, int chain2_num, double **TMave_mat,
    vector<vector<string> >&seqxA_mat, vector<vector<string> >&seqyA_mat,
    int *assign1_list, int *assign2_list, vector<string>&sequence,
    double d0_scale, bool fast_opt)
{
    /* tmp assignment */
    double total_score;
    int *assign1_tmp, *assign2_tmp;
    assign1_tmp=new int[chain1_num];
    assign2_tmp=new int[chain2_num];
    double **TMave_tmp;
    NewArray(&TMave_tmp,chain1_num,chain2_num);
    vector<string> tmp_str_vec(chain2_num,"");
    vector<vector<string> >seqxA_tmp(chain1_num,tmp_str_vec);
    vector<vector<string> >seqyA_tmp(chain1_num,tmp_str_vec);
    vector<string> sequence_tmp;
    copy_chain_assign_data(chain1_num, chain2_num, sequence_tmp,
        seqxA_mat, seqyA_mat, assign1_list, assign2_list, TMave_mat,
        seqxA_tmp, seqyA_tmp, assign1_tmp,  assign2_tmp,  TMave_tmp);
    
    for (int iter=0;iter<max_iter;iter++)
    {
        total_score=MMalign_search(xa_vec, ya_vec, seqx_vec, seqy_vec,
            secx_vec, secy_vec, mol_vec1, mol_vec2, xlen_vec, ylen_vec,
            xa, ya, seqx, seqy, secx, secy, len_aa, len_na,
            chain1_num, chain2_num, 
            TMave_tmp, seqxA_tmp, seqyA_tmp, assign1_tmp, assign2_tmp,
            sequence, d0_scale, fast_opt);
        total_score=enhanced_greedy_search(TMave_tmp, assign1_tmp,
            assign2_tmp, chain1_num, chain2_num);
        //if (total_score<=0) PrintErrorAndQuit("ERROR! No assignable chain");
        if (total_score<=max_total_score) break;
        max_total_score=total_score;
        copy_chain_assign_data(chain1_num, chain2_num, sequence,
            seqxA_tmp, seqyA_tmp, assign1_tmp,  assign2_tmp,  TMave_tmp,
            seqxA_mat, seqyA_mat, assign1_list, assign2_list, TMave_mat);
    }

    /* clean up everything */
    delete [] assign1_tmp;
    delete [] assign2_tmp;
    DeleteArray(&TMave_tmp,chain1_num);
    vector<string>().swap(tmp_str_vec);
    vector<vector<string> >().swap(seqxA_tmp);
    vector<vector<string> >().swap(seqyA_tmp);
}


/* Input: vectors x, y, rotation matrix t, u, scale factor d02, and gap_open
 * Output: j2i[1:len2] \in {1:len1} U {-1}
 * path[0:len1, 0:len2]=1,2,3, from diagonal, horizontal, vertical */
void NWDP_TM_dimer(bool **path, double **val, double **x, double **y,
    int len1, int len2, bool **mask,
    double t[3], double u[3][3], double d02, double gap_open, int j2i[])
{
    int i, j;
    double h, v, d;

    //initialization
    for(i=0; i<=len1; i++)
    {
        //val[i][0]=0;
        val[i][0]=i*gap_open;
        path[i][0]=false; //not from diagonal
    }

    for(j=0; j<=len2; j++)
    {
        //val[0][j]=0;
        val[0][j]=j*gap_open;
        path[0][j]=false; //not from diagonal
        j2i[j]=-1;    //all are not aligned, only use j2i[1:len2]
    }      
    double xx[3], dij;


    //decide matrix and path
    for(i=1; i<=len1; i++)
    {
        transform(t, u, &x[i-1][0], xx);
        for(j=1; j<=len2; j++)
        {
            d=FLT_MIN;
            if (mask[i][j])
            {
                dij=dist(xx, &y[j-1][0]);    
                d=val[i-1][j-1] +  1.0/(1+dij/d02);
            } 

            //symbol insertion in horizontal (= a gap in vertical)
            h=val[i-1][j];
            if(path[i-1][j]) h += gap_open; //aligned in last position

            //symbol insertion in vertical
            v=val[i][j-1];
            if(path[i][j-1]) v += gap_open; //aligned in last position


            if(d>=h && d>=v)
            {
                path[i][j]=true; //from diagonal
                val[i][j]=d;
            }
            else 
            {
                path[i][j]=false; //from horizontal
                if(v>=h) val[i][j]=v;
                else val[i][j]=h;
            }
        } //for i
    } //for j

    //trace back to extract the alignment
    i=len1;
    j=len2;
    while(i>0 && j>0)
    {
        if(path[i][j]) //from diagonal
        {
            j2i[j-1]=i-1;
            i--;
            j--;
        }
        else 
        {
            h=val[i-1][j];
            if(path[i-1][j]) h +=gap_open;

            v=val[i][j-1];
            if(path[i][j-1]) v +=gap_open;

            if(v>=h) j--;
            else i--;
        }
    }
}

/* +ss
 * Input: secondary structure secx, secy, and gap_open
 * Output: j2i[1:len2] \in {1:len1} U {-1}
 * path[0:len1, 0:len2]=1,2,3, from diagonal, horizontal, vertical */
void NWDP_TM_dimer(bool **path, double **val, const char *secx, const char *secy,
    const int len1, const int len2, bool **mask, const double gap_open, int j2i[])
{

    int i, j;
    double h, v, d;

    //initialization
    for(i=0; i<=len1; i++)
    {
        //val[i][0]=0;
        val[i][0]=i*gap_open;
        path[i][0]=false; //not from diagonal
    }

    for(j=0; j<=len2; j++)
    {
        //val[0][j]=0;
        val[0][j]=j*gap_open;
        path[0][j]=false; //not from diagonal
        j2i[j]=-1;    //all are not aligned, only use j2i[1:len2]
    }      

    //decide matrix and path
    for(i=1; i<=len1; i++)
    {
        for(j=1; j<=len2; j++)
        {
            d=FLT_MIN;
            if (mask[i][j])
                d=val[i-1][j-1] + 1.0*(secx[i-1]==secy[j-1]);

            //symbol insertion in horizontal (= a gap in vertical)
            h=val[i-1][j];
            if(path[i-1][j]) h += gap_open; //aligned in last position

            //symbol insertion in vertical
            v=val[i][j-1];
            if(path[i][j-1]) v += gap_open; //aligned in last position

            if(d>=h && d>=v)
            {
                path[i][j]=true; //from diagonal
                val[i][j]=d;
            }
            else 
            {
                path[i][j]=false; //from horizontal
                if(v>=h) val[i][j]=v;
                else val[i][j]=h;
            }
        } //for i
    } //for j

    //trace back to extract the alignment
    i=len1;
    j=len2;
    while(i>0 && j>0)
    {
        if(path[i][j]) //from diagonal
        {
            j2i[j-1]=i-1;
            i--;
            j--;
        }
        else 
        {
            h=val[i-1][j];
            if(path[i-1][j]) h +=gap_open;

            v=val[i][j-1];
            if(path[i][j-1]) v +=gap_open;

            if(v>=h) j--;
            else i--;
        }
    }
}

//heuristic run of dynamic programing iteratively to find the best alignment
//input: initial rotation matrix t, u
//       vectors x and y, d0
//output: best alignment that maximizes the TMscore, will be stored in invmap
double DP_iter_dimer(double **r1, double **r2, double **xtm, double **ytm,
    double **xt, bool **path, double **val, double **x, double **y,
    int xlen, int ylen, bool **mask, double t[3], double u[3][3], int invmap0[],
    int g1, int g2, int iteration_max, double local_d0_search,
    double D0_MIN, double Lnorm, double d0, double score_d8)
{
    double gap_open[2]={-0.6, 0};
    double rmsd; 
    int *invmap=new int[ylen+1];
    
    int iteration, i, j, k;
    double tmscore, tmscore_max, tmscore_old=0;    
    int score_sum_method=8, simplify_step=40;
    tmscore_max=-1;

    //double d01=d0+1.5;
    double d02=d0*d0;
    for(int g=g1; g<g2; g++)
    {
        for(iteration=0; iteration<iteration_max; iteration++)
        {           
            NWDP_TM_dimer(path, val, x, y, xlen, ylen, mask,
                t, u, d02, gap_open[g], invmap);
            
            k=0;
            for(j=0; j<ylen; j++) 
            {
                i=invmap[j];

                if(i>=0) //aligned
                {
                    xtm[k][0]=x[i][0];
                    xtm[k][1]=x[i][1];
                    xtm[k][2]=x[i][2];
                    
                    ytm[k][0]=y[j][0];
                    ytm[k][1]=y[j][1];
                    ytm[k][2]=y[j][2];
                    k++;
                }
            }

            tmscore = TMscore8_search(r1, r2, xtm, ytm, xt, k, t, u,
                simplify_step, score_sum_method, &rmsd, local_d0_search,
                Lnorm, score_d8, d0);

           
            if(tmscore>tmscore_max)
            {
                tmscore_max=tmscore;
                for(i=0; i<ylen; i++) invmap0[i]=invmap[i];
            }
    
            if(iteration>0)
            {
                if(fabs(tmscore_old-tmscore)<0.000001) break;       
            }
            tmscore_old=tmscore;
        }// for iteration           
        
    }//for gapopen
    
    
    delete []invmap;
    return tmscore_max;
}

void get_initial_ss_dimer(bool **path, double **val, const char *secx,
    const char *secy, int xlen, int ylen, bool **mask, int *y2x)
{
    double gap_open=-1.0;
    NWDP_TM_dimer(path, val, secx, secy, xlen, ylen, mask, gap_open, y2x);
}

bool get_initial5_dimer( double **r1, double **r2, double **xtm, double **ytm,
    bool **path, double **val, double **x, double **y, int xlen, int ylen,
    bool **mask, int *y2x,
    double d0, double d0_search, const bool fast_opt, const double D0_MIN)
{
    double GL, rmsd;
    double t[3];
    double u[3][3];

    double d01 = d0 + 1.5;
    if (d01 < D0_MIN) d01 = D0_MIN;
    double d02 = d01*d01;

    double GLmax = 0;
    int aL = getmin(xlen, ylen);
    int *invmap = new int[ylen + 1];

    // jump on sequence1-------------->
    int n_jump1 = 0;
    if (xlen > 250)
        n_jump1 = 45;
    else if (xlen > 200)
        n_jump1 = 35;
    else if (xlen > 150)
        n_jump1 = 25;
    else
        n_jump1 = 15;
    if (n_jump1 > (xlen / 3))
        n_jump1 = xlen / 3;

    // jump on sequence2-------------->
    int n_jump2 = 0;
    if (ylen > 250)
        n_jump2 = 45;
    else if (ylen > 200)
        n_jump2 = 35;
    else if (ylen > 150)
        n_jump2 = 25;
    else
        n_jump2 = 15;
    if (n_jump2 > (ylen / 3))
        n_jump2 = ylen / 3;

    // fragment to superimpose-------------->
    int n_frag[2] = { 20, 100 };
    if (n_frag[0] > (aL / 3))
        n_frag[0] = aL / 3;
    if (n_frag[1] > (aL / 2))
        n_frag[1] = aL / 2;

    // start superimpose search-------------->
    if (fast_opt)
    {
        n_jump1*=5;
        n_jump2*=5;
    }
    bool flag = false;
    for (int i_frag = 0; i_frag < 2; i_frag++)
    {
        int m1 = xlen - n_frag[i_frag] + 1;
        int m2 = ylen - n_frag[i_frag] + 1;

        for (int i = 0; i<m1; i = i + n_jump1) //index starts from 0, different from FORTRAN
        {
            for (int j = 0; j<m2; j = j + n_jump2)
            {
                for (int k = 0; k<n_frag[i_frag]; k++) //fragment in y
                {
                    r1[k][0] = x[k + i][0];
                    r1[k][1] = x[k + i][1];
                    r1[k][2] = x[k + i][2];

                    r2[k][0] = y[k + j][0];
                    r2[k][1] = y[k + j][1];
                    r2[k][2] = y[k + j][2];
                }

                // superpose the two structures and rotate it
                Kabsch(r1, r2, n_frag[i_frag], 1, &rmsd, t, u);

                double gap_open = 0.0;
                NWDP_TM_dimer(path, val, x, y, xlen, ylen, mask,
                    t, u, d02, gap_open, invmap);
                GL = get_score_fast(r1, r2, xtm, ytm, x, y, xlen, ylen,
                    invmap, d0, d0_search, t, u);
                if (GL>GLmax)
                {
                    GLmax = GL;
                    for (int ii = 0; ii<ylen; ii++) y2x[ii] = invmap[ii];
                    flag = true;
                }
            }
        }
    }

    delete[] invmap;
    return flag;
}

void get_initial_ssplus_dimer(double **r1, double **r2, double **score,
    bool **path, double **val, const char *secx, const char *secy,
    double **x, double **y, int xlen, int ylen, bool **mask,
    int *y2x0, int *y2x, const double D0_MIN, double d0)
{
    //create score matrix for DP
    score_matrix_rmsd_sec(r1, r2, score, secx, secy, x, y, xlen, ylen,
        y2x0, D0_MIN,d0);

    int i,j;
    for (i=0;i<xlen+1;i++) for (j=0;j<ylen+1;j++) score[i][j]=FLT_MIN;
    
    double gap_open=-1.0;
    NWDP_TM(score, path, val, xlen, ylen, gap_open, y2x);
}

/* Entry function for TM-align. Return TM-score calculation status:
 * 0   - full TM-score calculation 
 * 1   - terminated due to exception
 * 2-7 - pre-terminated due to low TM-score */
int TMalign_dimer_main(double **xa, double **ya,
    const char *seqx, const char *seqy, const char *secx, const char *secy,
    double t0[3], double u0[3][3],
    double &TM1, double &TM2, double &TM3, double &TM4, double &TM5,
    double &d0_0, double &TM_0,
    double &d0A, double &d0B, double &d0u, double &d0a, double &d0_out,
    string &seqM, string &seqxA, string &seqyA,
    double &rmsd0, int &L_ali, double &Liden,
    double &TM_ali, double &rmsd_ali, int &n_ali, int &n_ali8,
    const int xlen, const int ylen,
    bool **mask,
    const vector<string> sequence, const double Lnorm_ass,
    const double d0_scale, const int i_opt, const int a_opt,
    const bool u_opt, const bool d_opt, const bool fast_opt,
    const int mol_type, const double TMcut=-1)
{
    double D0_MIN;        //for d0
    double Lnorm;         //normalization length
    double score_d8,d0,d0_search,dcu0;//for TMscore search
    double t[3], u[3][3]; //Kabsch translation vector and rotation matrix
    double **score;       // Input score table for dynamic programming
    bool   **path;        // for dynamic programming  
    double **val;         // for dynamic programming  
    double **xtm, **ytm;  // for TMscore search engine
    double **xt;          //for saving the superposed version of r_1 or xtm
    double **r1, **r2;    // for Kabsch rotation

    /***********************/
    /* allocate memory     */
    /***********************/
    int minlen = min(xlen, ylen);
    NewArray(&score, xlen+1, ylen+1);
    NewArray(&path, xlen+1, ylen+1);
    NewArray(&val, xlen+1, ylen+1);
    NewArray(&xtm, minlen, 3);
    NewArray(&ytm, minlen, 3);
    NewArray(&xt, xlen, 3);
    NewArray(&r1, minlen, 3);
    NewArray(&r2, minlen, 3);

    /***********************/
    /*    parameter set    */
    /***********************/
    parameter_set4search(xlen, ylen, D0_MIN, Lnorm, 
        score_d8, d0, d0_search, dcu0);
    int simplify_step    = 40; //for simplified search engine
    int score_sum_method = 8;  //for scoring method, whether only sum over pairs with dis<score_d8

    int i;
    int *invmap0         = new int[ylen+1];
    int *invmap          = new int[ylen+1];
    double TM, TMmax=-1;
    for(i=0; i<ylen; i++) invmap0[i]=-1;

    double ddcc=0.4;
    if (Lnorm <= 40) ddcc=0.1;   //Lnorm was setted in parameter_set4search
    double local_d0_search = d0_search;

    //************************************************//
    //    get initial alignment from user's input:    //
    //    Stick to the initial alignment              //
    //************************************************//
    bool bAlignStick = false;
    if (i_opt==3)// if input has set parameter for "-I"
    {
        // In the original code, this loop starts from 1, which is
        // incorrect. Fortran starts from 1 but C++ should starts from 0.
        for (int j = 0; j < ylen; j++)// Set aligned position to be "-1"
            invmap[j] = -1;

        int i1 = -1;// in C version, index starts from zero, not from one
        int i2 = -1;
        int L1 = sequence[0].size();
        int L2 = sequence[1].size();
        int L = min(L1, L2);// Get positions for aligned residues
        for (int kk1 = 0; kk1 < L; kk1++)
        {
            if (sequence[0][kk1] != '-') i1++;
            if (sequence[1][kk1] != '-')
            {
                i2++;
                if (i2 >= ylen || i1 >= xlen) kk1 = L;
                else if (sequence[0][kk1] != '-') invmap[i2] = i1;
            }
        }

        //--------------- 2. Align proteins from original alignment
        double prevD0_MIN = D0_MIN;// stored for later use
        int prevLnorm = Lnorm;
        double prevd0 = d0;
        TM_ali = standard_TMscore(r1, r2, xtm, ytm, xt, xa, ya, xlen, ylen,
            invmap, L_ali, rmsd_ali, D0_MIN, Lnorm, d0, d0_search, score_d8,
            t, u, mol_type);
        D0_MIN = prevD0_MIN;
        Lnorm = prevLnorm;
        d0 = prevd0;
        TM = detailed_search_standard(r1, r2, xtm, ytm, xt, xa, ya, xlen, ylen,
            invmap, t, u, 40, 8, local_d0_search, true, Lnorm, score_d8, d0);
        if (TM > TMmax)
        {
            TMmax = TM;
            for (i = 0; i<ylen; i++) invmap0[i] = invmap[i];
        }
        bAlignStick = true;
    }

    /******************************************************/
    /*    get initial alignment with gapless threading    */
    /******************************************************/
    if (!bAlignStick)
    {
        get_initial(r1, r2, xtm, ytm, xa, ya, xlen, ylen, invmap0, d0,
            d0_search, fast_opt, t, u);
        TM = detailed_search(r1, r2, xtm, ytm, xt, xa, ya, xlen, ylen, invmap0,
            t, u, simplify_step, score_sum_method, local_d0_search, Lnorm,
            score_d8, d0);
        if (TM>TMmax) TMmax = TM;
        if (TMcut>0) copy_t_u(t, u, t0, u0);
        //run dynamic programing iteratively to find the best alignment
        TM = DP_iter_dimer(r1, r2, xtm, ytm, xt, path, val, xa, ya, xlen, ylen,
             mask, t, u, invmap, 0, 2, (fast_opt)?2:30,
             local_d0_search, D0_MIN, Lnorm, d0, score_d8);
        if (TM>TMmax)
        {
            TMmax = TM;
            for (int i = 0; i<ylen; i++) invmap0[i] = invmap[i];
            if (TMcut>0) copy_t_u(t, u, t0, u0);
        }

        if (TMcut>0) // pre-terminate if TM-score is too low
        {
            double TMtmp=approx_TM(xlen, ylen, a_opt,
                xa, ya, t0, u0, invmap0, mol_type);

            if (TMtmp<0.5*TMcut)
            {
                TM1=TM2=TM3=TM4=TM5=TMtmp;
                clean_up_after_approx_TM(invmap0, invmap, score, path, val,
                    xtm, ytm, xt, r1, r2, xlen, minlen);
                return 2;
            }
        }

        /************************************************************/
        /*    get initial alignment based on secondary structure    */
        /************************************************************/
        get_initial_ss_dimer(path, val, secx, secy, xlen, ylen, mask, invmap);
        TM = detailed_search(r1, r2, xtm, ytm, xt, xa, ya, xlen, ylen, invmap,
            t, u, simplify_step, score_sum_method, local_d0_search, Lnorm,
            score_d8, d0);
        if (TM>TMmax)
        {
            TMmax = TM;
            for (int i = 0; i<ylen; i++) invmap0[i] = invmap[i];
            if (TMcut>0) copy_t_u(t, u, t0, u0);
        }
        if (TM > TMmax*0.2)
        {
            TM = DP_iter_dimer(r1, r2, xtm, ytm, xt, path, val, xa, ya,
                xlen, ylen, mask, t, u, invmap, 0, 2,
                (fast_opt)?2:30, local_d0_search, D0_MIN, Lnorm, d0, score_d8);
            if (TM>TMmax)
            {
                TMmax = TM;
                for (int i = 0; i<ylen; i++) invmap0[i] = invmap[i];
                if (TMcut>0) copy_t_u(t, u, t0, u0);
            }
        }

        if (TMcut>0) // pre-terminate if TM-score is too low
        {
            double TMtmp=approx_TM(xlen, ylen, a_opt,
                xa, ya, t0, u0, invmap0, mol_type);

            if (TMtmp<0.52*TMcut)
            {
                TM1=TM2=TM3=TM4=TM5=TMtmp;
                clean_up_after_approx_TM(invmap0, invmap, score, path, val,
                    xtm, ytm, xt, r1, r2, xlen, minlen);
                return 3;
            }
        }

        /************************************************************/
        /*    get initial alignment based on local superposition    */
        /************************************************************/
        //=initial5 in original TM-align
        if (get_initial5_dimer( r1, r2, xtm, ytm, path, val, xa, ya,
            xlen, ylen, mask, invmap, d0, d0_search, fast_opt, D0_MIN))
        {
            TM = detailed_search(r1, r2, xtm, ytm, xt, xa, ya, xlen, ylen,
                invmap, t, u, simplify_step, score_sum_method,
                local_d0_search, Lnorm, score_d8, d0);
            if (TM>TMmax)
            {
                TMmax = TM;
                for (int i = 0; i<ylen; i++) invmap0[i] = invmap[i];
                if (TMcut>0) copy_t_u(t, u, t0, u0);
            }
            if (TM > TMmax*ddcc)
            {
                TM = DP_iter_dimer(r1, r2, xtm, ytm, xt, path, val, xa, ya,
                    xlen, ylen, mask, t, u, invmap, 0, 2, 2,
                    local_d0_search, D0_MIN, Lnorm, d0, score_d8);
                if (TM>TMmax)
                {
                    TMmax = TM;
                    for (int i = 0; i<ylen; i++) invmap0[i] = invmap[i];
                    if (TMcut>0) copy_t_u(t, u, t0, u0);
                }
            }
        }
        else
            cerr << "\n\nWarning: initial alignment from local superposition fail!\n\n" << endl;

        if (TMcut>0) // pre-terminate if TM-score is too low
        {
            double TMtmp=approx_TM(xlen, ylen, a_opt,
                xa, ya, t0, u0, invmap0, mol_type);

            if (TMtmp<0.54*TMcut)
            {
                TM1=TM2=TM3=TM4=TM5=TMtmp;
                clean_up_after_approx_TM(invmap0, invmap, score, path, val,
                    xtm, ytm, xt, r1, r2, xlen, minlen);
                return 4;
            }
        }

        /********************************************************************/
        /* get initial alignment by local superposition+secondary structure */
        /********************************************************************/
        //=initial3 in original TM-align
        get_initial_ssplus_dimer(r1, r2, score, path, val, secx, secy, xa, ya,
            xlen, ylen, mask, invmap0, invmap, D0_MIN, d0);
        TM = detailed_search(r1, r2, xtm, ytm, xt, xa, ya, xlen, ylen, invmap,
             t, u, simplify_step, score_sum_method, local_d0_search, Lnorm,
             score_d8, d0);
        if (TM>TMmax)
        {
            TMmax = TM;
            for (i = 0; i<ylen; i++) invmap0[i] = invmap[i];
            if (TMcut>0) copy_t_u(t, u, t0, u0);
        }
        if (TM > TMmax*ddcc)
        {
            TM = DP_iter_dimer(r1, r2, xtm, ytm, xt, path, val, xa, ya,
                xlen, ylen, mask, t, u, invmap, 0, 2,
                (fast_opt)?2:30, local_d0_search, D0_MIN, Lnorm, d0, score_d8);
            if (TM>TMmax)
            {
                TMmax = TM;
                for (i = 0; i<ylen; i++) invmap0[i] = invmap[i];
                if (TMcut>0) copy_t_u(t, u, t0, u0);
            }
        }

        if (TMcut>0) // pre-terminate if TM-score is too low
        {
            double TMtmp=approx_TM(xlen, ylen, a_opt,
                xa, ya, t0, u0, invmap0, mol_type);

            if (TMtmp<0.56*TMcut)
            {
                TM1=TM2=TM3=TM4=TM5=TMtmp;
                clean_up_after_approx_TM(invmap0, invmap, score, path, val,
                    xtm, ytm, xt, r1, r2, xlen, minlen);
                return 5;
            }
        }

        /*******************************************************************/
        /*    get initial alignment based on fragment gapless threading    */
        /*******************************************************************/
        //=initial4 in original TM-align
        get_initial_fgt(r1, r2, xtm, ytm, xa, ya, xlen, ylen,
            invmap, d0, d0_search, dcu0, fast_opt, t, u);
        TM = detailed_search(r1, r2, xtm, ytm, xt, xa, ya, xlen, ylen, invmap,
            t, u, simplify_step, score_sum_method, local_d0_search, Lnorm,
            score_d8, d0);
        if (TM>TMmax)
        {
            TMmax = TM;
            for (i = 0; i<ylen; i++) invmap0[i] = invmap[i];
            if (TMcut>0) copy_t_u(t, u, t0, u0);
        }
        if (TM > TMmax*ddcc)
        {
            TM = DP_iter_dimer(r1, r2, xtm, ytm, xt, path, val, xa, ya,
                xlen, ylen, mask, t, u, invmap, 1, 2, 2,
                local_d0_search, D0_MIN, Lnorm, d0, score_d8);
            if (TM>TMmax)
            {
                TMmax = TM;
                for (i = 0; i<ylen; i++) invmap0[i] = invmap[i];
                if (TMcut>0) copy_t_u(t, u, t0, u0);
            }
        }

        if (TMcut>0) // pre-terminate if TM-score is too low
        {
            double TMtmp=approx_TM(xlen, ylen, a_opt,
                xa, ya, t0, u0, invmap0, mol_type);

            if (TMtmp<0.58*TMcut)
            {
                TM1=TM2=TM3=TM4=TM5=TMtmp;
                clean_up_after_approx_TM(invmap0, invmap, score, path, val,
                    xtm, ytm, xt, r1, r2, xlen, minlen);
                return 6;
            }
        }

        //************************************************//
        //    get initial alignment from user's input:    //
        //************************************************//
        if (i_opt==1)// if input has set parameter for "-i"
        {
            for (int j = 0; j < ylen; j++)// Set aligned position to be "-1"
                invmap[j] = -1;

            int i1 = -1;// in C version, index starts from zero, not from one
            int i2 = -1;
            int L1 = sequence[0].size();
            int L2 = sequence[1].size();
            int L = min(L1, L2);// Get positions for aligned residues
            for (int kk1 = 0; kk1 < L; kk1++)
            {
                if (sequence[0][kk1] != '-')
                    i1++;
                if (sequence[1][kk1] != '-')
                {
                    i2++;
                    if (i2 >= ylen || i1 >= xlen) kk1 = L;
                    else if (sequence[0][kk1] != '-') invmap[i2] = i1;
                }
            }

            //--------------- 2. Align proteins from original alignment
            double prevD0_MIN = D0_MIN;// stored for later use
            int prevLnorm = Lnorm;
            double prevd0 = d0;
            TM_ali = standard_TMscore(r1, r2, xtm, ytm, xt, xa, ya,
                xlen, ylen, invmap, L_ali, rmsd_ali, D0_MIN, Lnorm, d0,
                d0_search, score_d8, t, u, mol_type);
            D0_MIN = prevD0_MIN;
            Lnorm = prevLnorm;
            d0 = prevd0;

            TM = detailed_search_standard(r1, r2, xtm, ytm, xt, xa, ya,
                xlen, ylen, invmap, t, u, 40, 8, local_d0_search, true, Lnorm,
                score_d8, d0);
            if (TM > TMmax)
            {
                TMmax = TM;
                for (i = 0; i<ylen; i++) invmap0[i] = invmap[i];
            }
            // Different from get_initial, get_initial_ss and get_initial_ssplus
            TM = DP_iter_dimer(r1, r2, xtm, ytm, xt, path, val, xa, ya,
                xlen, ylen, mask, t, u, invmap, 0, 2,
                (fast_opt)?2:30, local_d0_search, D0_MIN, Lnorm, d0, score_d8);
            if (TM>TMmax)
            {
                TMmax = TM;
                for (i = 0; i<ylen; i++) invmap0[i] = invmap[i];
            }
        }
    }



    //*******************************************************************//
    //    The alignment will not be changed any more in the following    //
    //*******************************************************************//
    //check if the initial alignment is generated appropriately
    bool flag=false;
    for(i=0; i<ylen; i++)
    {
        if(invmap0[i]>=0)
        {
            flag=true;
            break;
        }
    }
    if(!flag)
    {
        cout << "There is no alignment between the two structures! "
             << "Program stop with no result!" << endl;
        TM1=TM2=TM3=TM4=TM5=0;
        return 1;
    }

    /* last TM-score pre-termination */
    if (TMcut>0)
    {
        double TMtmp=approx_TM(xlen, ylen, a_opt,
            xa, ya, t0, u0, invmap0, mol_type);

        if (TMtmp<0.6*TMcut)
        {
            TM1=TM2=TM3=TM4=TM5=TMtmp;
            clean_up_after_approx_TM(invmap0, invmap, score, path, val,
                xtm, ytm, xt, r1, r2, xlen, minlen);
            return 7;
        }
    }

    //********************************************************************//
    //    Detailed TMscore search engine --> prepare for final TMscore    //
    //********************************************************************//
    //run detailed TMscore search engine for the best alignment, and
    //extract the best rotation matrix (t, u) for the best alignment
    simplify_step=1;
    if (fast_opt) simplify_step=40;
    score_sum_method=8;
    TM = detailed_search_standard(r1, r2, xtm, ytm, xt, xa, ya, xlen, ylen,
        invmap0, t, u, simplify_step, score_sum_method, local_d0_search,
        false, Lnorm, score_d8, d0);

    //select pairs with dis<d8 for final TMscore computation and output alignment
    int k=0;
    int *m1, *m2;
    double d;
    m1=new int[xlen]; //alignd index in x
    m2=new int[ylen]; //alignd index in y
    do_rotation(xa, xt, xlen, t, u);
    k=0;
    for(int j=0; j<ylen; j++)
    {
        i=invmap0[j];
        if(i>=0)//aligned
        {
            n_ali++;
            d=sqrt(dist(&xt[i][0], &ya[j][0]));
            if (d <= score_d8 || (i_opt == 3))
            {
                m1[k]=i;
                m2[k]=j;

                xtm[k][0]=xa[i][0];
                xtm[k][1]=xa[i][1];
                xtm[k][2]=xa[i][2];

                ytm[k][0]=ya[j][0];
                ytm[k][1]=ya[j][1];
                ytm[k][2]=ya[j][2];

                r1[k][0] = xt[i][0];
                r1[k][1] = xt[i][1];
                r1[k][2] = xt[i][2];
                r2[k][0] = ya[j][0];
                r2[k][1] = ya[j][1];
                r2[k][2] = ya[j][2];

                k++;
            }
        }
    }
    n_ali8=k;

    Kabsch(r1, r2, n_ali8, 0, &rmsd0, t, u);// rmsd0 is used for final output, only recalculate rmsd0, not t & u
    rmsd0 = sqrt(rmsd0 / n_ali8);


    //****************************************//
    //              Final TMscore             //
    //    Please set parameters for output    //
    //****************************************//
    double rmsd;
    simplify_step=1;
    score_sum_method=0;
    double Lnorm_0=ylen;


    //normalized by length of structure A
    parameter_set4final(Lnorm_0, D0_MIN, Lnorm, d0, d0_search, mol_type);
    d0A=d0;
    d0_0=d0A;
    local_d0_search = d0_search;
    TM1 = TMscore8_search(r1, r2, xtm, ytm, xt, n_ali8, t0, u0, simplify_step,
        score_sum_method, &rmsd, local_d0_search, Lnorm, score_d8, d0);
    TM_0 = TM1;

    //normalized by length of structure B
    parameter_set4final(xlen+0.0, D0_MIN, Lnorm, d0, d0_search, mol_type);
    d0B=d0;
    local_d0_search = d0_search;
    TM2 = TMscore8_search(r1, r2, xtm, ytm, xt, n_ali8, t, u, simplify_step,
        score_sum_method, &rmsd, local_d0_search, Lnorm, score_d8, d0);

    double Lnorm_d0;
    if (a_opt>0)
    {
        //normalized by average length of structures A, B
        Lnorm_0=(xlen+ylen)*0.5;
        parameter_set4final(Lnorm_0, D0_MIN, Lnorm, d0, d0_search, mol_type);
        d0a=d0;
        d0_0=d0a;
        local_d0_search = d0_search;

        TM3 = TMscore8_search(r1, r2, xtm, ytm, xt, n_ali8, t0, u0,
            simplify_step, score_sum_method, &rmsd, local_d0_search, Lnorm,
            score_d8, d0);
        TM_0=TM3;
    }
    if (u_opt)
    {
        //normalized by user assigned length
        parameter_set4final(Lnorm_ass, D0_MIN, Lnorm,
            d0, d0_search, mol_type);
        d0u=d0;
        d0_0=d0u;
        Lnorm_0=Lnorm_ass;
        local_d0_search = d0_search;
        TM4 = TMscore8_search(r1, r2, xtm, ytm, xt, n_ali8, t0, u0,
            simplify_step, score_sum_method, &rmsd, local_d0_search, Lnorm,
            score_d8, d0);
        TM_0=TM4;
    }
    if (d_opt)
    {
        //scaled by user assigned d0
        parameter_set4scale(ylen, d0_scale, Lnorm, d0, d0_search);
        d0_out=d0_scale;
        d0_0=d0_scale;
        //Lnorm_0=ylen;
        Lnorm_d0=Lnorm_0;
        local_d0_search = d0_search;
        TM5 = TMscore8_search(r1, r2, xtm, ytm, xt, n_ali8, t0, u0,
            simplify_step, score_sum_method, &rmsd, local_d0_search, Lnorm,
            score_d8, d0);
        TM_0=TM5;
    }

    /* derive alignment from superposition */
    int ali_len=xlen+ylen; //maximum length of alignment
    seqxA.assign(ali_len,'-');
    seqM.assign( ali_len,' ');
    seqyA.assign(ali_len,'-');
    
    //do_rotation(xa, xt, xlen, t, u);
    do_rotation(xa, xt, xlen, t0, u0);

    int kk=0, i_old=0, j_old=0;
    d=0;
    for(int k=0; k<n_ali8; k++)
    {
        for(int i=i_old; i<m1[k]; i++)
        {
            //align x to gap
            seqxA[kk]=seqx[i];
            seqyA[kk]='-';
            seqM[kk]=' ';                    
            kk++;
        }

        for(int j=j_old; j<m2[k]; j++)
        {
            //align y to gap
            seqxA[kk]='-';
            seqyA[kk]=seqy[j];
            seqM[kk]=' ';
            kk++;
        }

        seqxA[kk]=seqx[m1[k]];
        seqyA[kk]=seqy[m2[k]];
        Liden+=(seqxA[kk]==seqyA[kk]);
        d=sqrt(dist(&xt[m1[k]][0], &ya[m2[k]][0]));
        if(d<d0_out) seqM[kk]=':';
        else         seqM[kk]='.';
        kk++;  
        i_old=m1[k]+1;
        j_old=m2[k]+1;
    }

    //tail
    for(int i=i_old; i<xlen; i++)
    {
        //align x to gap
        seqxA[kk]=seqx[i];
        seqyA[kk]='-';
        seqM[kk]=' ';
        kk++;
    }    
    for(int j=j_old; j<ylen; j++)
    {
        //align y to gap
        seqxA[kk]='-';
        seqyA[kk]=seqy[j];
        seqM[kk]=' ';
        kk++;
    }
    seqxA=seqxA.substr(0,kk);
    seqyA=seqyA.substr(0,kk);
    seqM =seqM.substr(0,kk);

    /* free memory */
    clean_up_after_approx_TM(invmap0, invmap, score, path, val,
        xtm, ytm, xt, r1, r2, xlen, minlen);
    delete [] m1;
    delete [] m2;
    return 0; // zero for no exception
}

void MMalign_dimer(double & total_score, 
    const vector<vector<vector<double> > >&xa_vec,
    const vector<vector<vector<double> > >&ya_vec,
    const vector<vector<char> >&seqx_vec, const vector<vector<char> >&seqy_vec,
    const vector<vector<char> >&secx_vec, const vector<vector<char> >&secy_vec,
    const vector<int> &mol_vec1, const vector<int> &mol_vec2,
    const vector<int> &xlen_vec, const vector<int> &ylen_vec,
    double **xa, double **ya, char *seqx, char *seqy, char *secx, char *secy,
    int len_aa, int len_na, int chain1_num, int chain2_num, double **TMave_mat,
    vector<vector<string> >&seqxA_mat, vector<vector<string> >&seqyA_mat,
    int *assign1_list, int *assign2_list, vector<string>&sequence,
    double d0_scale, bool fast_opt)
{
    int i,j;
    int xlen=0;
    int ylen=0;
    vector<int> xlen_dimer;
    vector<int> ylen_dimer;
    for (i=0;i<chain1_num;i++)
    {
        j=assign1_list[i];
        if (j<0) continue;
        xlen+=xlen_vec[i];
        ylen+=ylen_vec[j];
        xlen_dimer.push_back(xlen_vec[i]);
        ylen_dimer.push_back(ylen_vec[j]);
    }
    if (xlen<=3 || ylen<=3) return;

    bool **mask; // mask out inter-chain region
    NewArray(&mask, xlen+1, ylen+1);
    for (i=0;i<xlen+1;i++) for (j=0;j<ylen+1;j++) mask[i][j]=false;
    for (i=0;i<xlen_dimer[0]+1;i++) mask[i][0]=true;
    for (j=0;j<ylen_dimer[0]+1;j++) mask[0][j]=true;
    int c,prev_xlen,prev_ylen;
    prev_xlen=1;
    prev_ylen=1;
    for (c=0;c<xlen_dimer.size();c++)
    {
        for (i=prev_xlen;i<prev_xlen+xlen_dimer[c];i++)
            for (j=prev_ylen;j<prev_ylen+ylen_dimer[c];j++) mask[i][j]=true;
        prev_xlen+=xlen_dimer[c];
        prev_ylen+=ylen_dimer[c];
    }
    vector<int>().swap(xlen_dimer);
    vector<int>().swap(ylen_dimer);

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

    TMalign_dimer_main(xa, ya, seqx, seqy, secx, secy,
        t0, u0, TM1, TM2, TM3, TM4, TM5,
        d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out, seqM, seqxA, seqyA,
        rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
        xlen, ylen, mask, sequence, Lnorm_ass, d0_scale,
        1, false, true, false, fast_opt, mol_type, -1);

    /* clean up TM-align */
    delete [] seqx;
    delete [] seqy;
    delete [] secx;
    delete [] secy;
    DeleteArray(&xa,xlen);
    DeleteArray(&ya,ylen);
    DeleteArray(&mask,xlen+1);

    /* re-compute chain level alignment */
    total_score=0;
    for (i=0;i<chain1_num;i++)
    {
        xlen=xlen_vec[i];
        if (xlen<3)
        {
            for (j=0;j<chain2_num;j++) TMave_mat[i][j]=-1;
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
                TMave_mat[i][j]=-1;
                continue;
            }

            ylen=ylen_vec[j];
            if (ylen<3)
            {
                TMave_mat[i][j]=-1;
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
                0, false, 2, false, mol_vec1[i]+mol_vec2[j], 1, invmap);

            /* print result */
            seqxA_mat[i][j]=seqxA;
            seqyA_mat[i][j]=seqyA;

            TMave_mat[i][j]=TM4*Lnorm_ass;
            if (assign1_list[i]==j)
            {
                if (TM4<=0) assign1_list[i]=assign2_list[j]=-1;
                else        total_score+=TMave_mat[i][j];
            }

            /* clean up */
            seqM.clear();
            seqxA.clear();
            seqyA.clear();

            delete[]seqy;
            delete[]secy;
            DeleteArray(&ya,ylen);
            delete[]invmap;
        }
        delete[]seqx;
        delete[]secx;
        DeleteArray(&xa,xlen);
        DeleteArray(&xt,xlen);
    }
    return;
}

void MMalign_cross(double & max_total_score, const int max_iter,
    const vector<vector<vector<double> > >&xa_vec,
    const vector<vector<vector<double> > >&ya_vec,
    const vector<vector<char> >&seqx_vec, const vector<vector<char> >&seqy_vec,
    const vector<vector<char> >&secx_vec, const vector<vector<char> >&secy_vec,
    const vector<int> &mol_vec1, const vector<int> &mol_vec2,
    const vector<int> &xlen_vec, const vector<int> &ylen_vec,
    double **xa, double **ya, char *seqx, char *seqy, char *secx, char *secy,
    int len_aa, int len_na, int chain1_num, int chain2_num, double **TMave_mat,
    vector<vector<string> >&seqxA_mat, vector<vector<string> >&seqyA_mat,
    int *assign1_list, int *assign2_list, vector<string>&sequence,
    double d0_scale, bool fast_opt)
{
    /* tmp assignment */
    int *assign1_tmp, *assign2_tmp;
    assign1_tmp=new int[chain1_num];
    assign2_tmp=new int[chain2_num];
    double **TMave_tmp;
    NewArray(&TMave_tmp,chain1_num,chain2_num);
    vector<string> tmp_str_vec(chain2_num,"");
    vector<vector<string> >seqxA_tmp(chain1_num,tmp_str_vec);
    vector<vector<string> >seqyA_tmp(chain1_num,tmp_str_vec);
    vector<string> sequence_tmp;
    copy_chain_assign_data(chain1_num, chain2_num, sequence_tmp,
        seqxA_mat, seqyA_mat, assign1_list, assign2_list, TMave_mat,
        seqxA_tmp, seqyA_tmp, assign1_tmp,  assign2_tmp,  TMave_tmp);

    double total_score=MMalign_search(xa_vec, ya_vec, seqx_vec, seqy_vec,
        secx_vec, secy_vec, mol_vec1, mol_vec2, xlen_vec, ylen_vec,
        xa, ya, seqx, seqy, secx, secy, len_aa, len_na, chain1_num, chain2_num,
        TMave_tmp, seqxA_tmp, seqyA_tmp, assign1_tmp, assign2_tmp, sequence_tmp,
        d0_scale, fast_opt, 1);
    if (total_score>max_total_score)
    {
        copy_chain_assign_data(chain1_num, chain2_num, sequence,
            seqxA_tmp, seqyA_tmp, assign1_tmp,  assign2_tmp,  TMave_tmp,
            seqxA_mat, seqyA_mat, assign1_list, assign2_list, TMave_mat);
        max_total_score=total_score;
    }

    if (max_iter) MMalign_iter(
        max_total_score, max_iter, xa_vec, ya_vec, seqx_vec, seqy_vec,
        secx_vec, secy_vec, mol_vec1, mol_vec2, xlen_vec, ylen_vec,
        xa, ya, seqx, seqy, secx, secy, len_aa, len_na, chain1_num, chain2_num,
        TMave_mat, seqxA_mat, seqyA_mat, assign1_list, assign2_list, sequence,
        d0_scale, fast_opt);

    /* clean up everything */
    delete [] assign1_tmp;
    delete [] assign2_tmp;
    DeleteArray(&TMave_tmp,chain1_num);
    vector<string>().swap(tmp_str_vec);
    vector<vector<string> >().swap(seqxA_tmp);
    vector<vector<string> >().swap(seqyA_tmp);
    vector<string>().swap(sequence_tmp);
    return;
}

/* return the number of chains that are trimmed */
int trimComplex(vector<vector<vector<double> > >&a_trim_vec,
    vector<vector<char> >&seq_trim_vec, vector<vector<char> >&sec_trim_vec,
    vector<int>&len_trim_vec,
    const vector<vector<vector<double> > >&a_vec,
    const vector<vector<char> >&seq_vec, const vector<vector<char> >&sec_vec,
    const vector<int> &len_vec, const vector<int> &mol_vec,
    const int Lchain_aa_max, const int Lchain_na_max)
{
    int trim_chain_count=0;
    int chain_num=a_vec.size();
    int i,j;
    int r1,r2;
    double dinter;
    double dinter_min;
    vector<pair<double,int> >dinter_vec;
    vector<bool> include_vec;
    vector<char> seq_empty;
    vector<vector<double> >  a_empty;
    vector<double> xcoor(3,0);
    vector<double> ycoor(3,0);
    int xlen,ylen;
    int Lchain_max;
    double expand=2;
    for (i=0;i<chain_num;i++)
    {
        xlen=len_vec[i];
        if (mol_vec[i]>0) Lchain_max=Lchain_na_max*expand;
        else              Lchain_max=Lchain_aa_max*expand;
        if (Lchain_max<3) Lchain_max=3;
        if (xlen<=Lchain_max || xlen<=3)
        {
            a_trim_vec.push_back(a_vec[i]);
            seq_trim_vec.push_back(seq_vec[i]);
            sec_trim_vec.push_back(sec_vec[i]);
            len_trim_vec.push_back(xlen);
            continue;
        }
        trim_chain_count++;
        for (r1=0;r1<xlen;r1++)
        {
            xcoor[0]=a_vec[i][r1][0];
            xcoor[1]=a_vec[i][r1][1];
            xcoor[2]=a_vec[i][r1][2];
            dinter_min=FLT_MAX;
            for (j=0;j<chain_num;j++)
            {
                if (i==j) continue;
                ylen=len_vec[j];
                for (r2=0;r2<ylen;r2++)
                {
                    ycoor[0]=a_vec[j][r2][0];
                    ycoor[1]=a_vec[j][r2][1];
                    ycoor[2]=a_vec[j][r2][2];
                    dinter=(xcoor[0]-ycoor[0])*(xcoor[0]-ycoor[0])+
                           (xcoor[1]-ycoor[1])*(xcoor[1]-ycoor[1])+
                           (xcoor[2]-ycoor[2])*(xcoor[2]-ycoor[2]);
                    if (dinter<dinter_min) dinter_min=dinter;
                }
            }
            dinter_vec.push_back(make_pair(dinter,r1));
        }
        sort(dinter_vec.begin(),dinter_vec.end());
        include_vec.assign(xlen,false);
        for (r1=0;r1<Lchain_max;r1++)
            include_vec[dinter_vec[r1].second]=true;
        dinter_vec.clear();

        a_trim_vec.push_back(a_empty);
        seq_trim_vec.push_back(seq_empty);
        sec_trim_vec.push_back(seq_empty);
        len_trim_vec.push_back(Lchain_max);
        for (r1=0;r1<xlen;r1++)
        {
            if (include_vec[r1]==false) continue;
            a_trim_vec[i].push_back(a_vec[i][r1]);
            seq_trim_vec[i].push_back(seq_vec[i][r1]);
            sec_trim_vec[i].push_back(sec_vec[i][r1]);
        }
        include_vec.clear();
    }
    vector<pair<double,int> >().swap(dinter_vec);
    vector<bool>().swap(include_vec);
    vector<double> ().swap(xcoor);
    vector<double> ().swap(ycoor);
    return trim_chain_count;
}

void writeTrimComplex(vector<vector<vector<double> > >&a_trim_vec,
    vector<vector<char> >&seq_trim_vec, vector<int>&len_trim_vec,
    vector<string>&chainID_list, vector<int>&mol_vec,
    const string &atom_opt, string filename)
{
    int c,r;
    int a=0;
    string chainID;
    string atom;
    ofstream fp(filename.c_str());
    for (c=0;c<chainID_list.size();c++)
    {
        chainID=chainID_list[c];
        if (chainID.size()==1) chainID=" "+chainID;
        else if (chainID.size()>2) chainID=chainID.substr(chainID.size()-2,2);
        if (chainID[0]==':') chainID=" "+chainID.substr(1);
        atom=atom_opt;
        if (atom_opt=="auto")
        {
            if (mol_vec[c]>0) atom=" C3'";
            else              atom=" CA ";
        }

        for (r=0;r<len_trim_vec[c];r++)
            fp<<"ATOM  "<<resetiosflags(ios::right)<<setw(5)<<++a<<' '
              <<atom<<' '<<AAmap(seq_trim_vec[c][r])<<chainID
              <<setw(4)<<r+1<<"    "
              <<setiosflags(ios::fixed)<<setprecision(3)
              <<setw(8)<<a_trim_vec[c][r][0]
              <<setw(8)<<a_trim_vec[c][r][1]
              <<setw(8)<<a_trim_vec[c][r][2]<<endl;
    }
    fp.close();
    atom.clear();
    chainID.clear();
    return;
}

void output_dock_rotation_matrix(const char* fname_matrix,
    const vector<string>&xname_vec, const vector<string>&yname_vec,
    double ** ut_mat, int *assign1_list)
{
    fstream fout;
    fout.open(fname_matrix, ios::out | ios::trunc);
    if (fout)// succeed
    {
        int i,k;
        for (i=0;i<xname_vec.size();i++)
        {
            if (assign1_list[i]<0) continue;
            fout << "------ The rotation matrix to rotate "
                 <<xname_vec[i]<<" to "<<yname_vec[i]<<" ------\n"
                 << "m               t[m]        u[m][0]        u[m][1]        u[m][2]\n";
            for (k = 0; k < 3; k++)
                fout<<k<<setiosflags(ios::fixed)<<setprecision(10)
                    <<' '<<setw(18)<<ut_mat[i][9+k]
                    <<' '<<setw(14)<<ut_mat[i][3*k+0]
                    <<' '<<setw(14)<<ut_mat[i][3*k+1]
                    <<' '<<setw(14)<<ut_mat[i][3*k+2]<<'\n';
        }
        fout << "\nCode for rotating Structure 1 from (x,y,z) to (X,Y,Z):\n"
                "for(i=0; i<L; i++)\n"
                "{\n"
                "   X[i] = t[0] + u[0][0]*x[i] + u[0][1]*y[i] + u[0][2]*z[i];\n"
                "   Y[i] = t[1] + u[1][0]*x[i] + u[1][1]*y[i] + u[1][2]*z[i];\n"
                "   Z[i] = t[2] + u[2][0]*x[i] + u[2][1]*y[i] + u[2][2]*z[i];\n"
                "}"<<endl;
        fout.close();
    }
    else
        cout << "Open file to output rotation matrix fail.\n";
}
