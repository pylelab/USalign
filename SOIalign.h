#ifndef SOIalign_h
#define SOIalign_h 1

#include "TMalign.h"

void print_invmap(int *invmap, const int ylen)
{
    int i,j;
    for (j=0;j<ylen;j++)
    {
        i=invmap[j];
        if (i>=0) cout<<" ("<<i<<","<<j<<")";
    }
    cout<<endl;
}

void assign_sec_bond(int **secx_bond, const char *secx, const int xlen)
{
    int i,j;
    int starti=-1;
    int endi=-1;
    char ss;
    char prev_ss=0;
    for (i=0; i<xlen; i++)
    {
        ss=secx[i];
        secx_bond[i][0]=secx_bond[i][1]=-1;
        if (ss!=prev_ss && !(ss=='C' && prev_ss=='T') 
                        && !(ss=='T' && prev_ss=='C'))
        {
            if (starti>=0) // previous SSE end
            {
                endi=i;
                for (j=starti;j<endi;j++)
                {
                    secx_bond[j][0]=starti;
                    secx_bond[j][1]=endi;
                }
            }
            if (ss=='H' || ss=='E' || ss=='<' || ss=='>') starti=i;
            else starti=-1;
        }
        prev_ss=secx[i];
    }
    if (starti>=0) // previous SSE end
    {
        endi=i;
        for (j=starti;j<endi;j++)
        {
            secx_bond[j][0]=starti;
            secx_bond[j][1]=endi;
        }
    }
    for (i=0;i<xlen;i++) if (secx_bond[i][1]-secx_bond[i][0]==1)
        secx_bond[i][0]=secx_bond[i][1]=-1;
}

void getCloseK(double **xa, const int xlen, const int closeK_opt, double **xk)
{
    double **score;
    NewArray(&score, xlen+1, xlen+1);
    vector<pair<double,int> > close_idx_vec(xlen, make_pair(0,0));
    int i,j,k;
    for (i=0;i<xlen;i++)
    {
        score[i+1][i+1]=0;
        for (j=i+1;j<xlen;j++) score[j+1][i+1]=score[i+1][j+1]=dist(xa[i], xa[j]);
    }
    for (i=0;i<xlen;i++)
    {
        for (j=0;j<xlen;j++)
        {
            close_idx_vec[j].first=score[i+1][j+1];
            close_idx_vec[j].second=j;
        }
        sort(close_idx_vec.begin(), close_idx_vec.end());
        for (k=0;k<closeK_opt;k++)
        {
            j=close_idx_vec[k % xlen].second;
            xk[i*closeK_opt+k][0]=xa[j][0];
            xk[i*closeK_opt+k][1]=xa[j][1];
            xk[i*closeK_opt+k][2]=xa[j][2];
        }
    }

    /* clean up */
    vector<pair<double,int> >().swap(close_idx_vec);
    DeleteArray(&score, xlen+1);
}

/* check if pairing i to j conform to sequantiality within the SSE */
inline bool sec2sq(const int i, const int j,
    int **secx_bond, int **secy_bond, int *fwdmap, int *invmap)
{
    if (i<0 || j<0) return true;
    int ii,jj;
    if (secx_bond[i][0]>=0)
    {
        for (ii=secx_bond[i][0];ii<secx_bond[i][1];ii++)
        {
            jj=fwdmap[ii];
            if (jj>=0 && (i-ii)*(j-jj)<=0) return false;
        }
    }
    if (secy_bond[j][0]>=0)
    {
        for (jj=secy_bond[j][0];jj<secy_bond[j][1];jj++)
        {
            ii=invmap[jj];
            if (ii>=0 && (i-ii)*(j-jj)<=0) return false;
        }
    }
    return true;
}

void soi_egs(double **score, const int xlen, const int ylen, int *invmap,
    int **secx_bond, int **secy_bond, const int mm_opt)
{
    int i,j;
    int *fwdmap=new int[xlen]; // j=fwdmap[i];
    for (i=0; i<xlen; i++) fwdmap[i]=-1;
    for (j=0; j<ylen; j++)
    {
        i=invmap[j];
        if (i>=0) fwdmap[i]=j;
    }

    /* stage 1 - make initial assignment, starting from the highest score pair */
    double max_score;
    int maxi,maxj;
    while(1)
    {
        max_score=0;
        maxi=maxj=-1;
        for (i=0;i<xlen;i++)
        {
            if (fwdmap[i]>=0) continue;
            for (j=0;j<ylen;j++)
            {
                if (invmap[j]>=0 || score[i+1][j+1]<=max_score) continue;
                if (mm_opt==6 && !sec2sq(i,j,secx_bond,secy_bond,
                    fwdmap,invmap)) continue;
                maxi=i;
                maxj=j;
                max_score=score[i+1][j+1];
            }
        }
        if (maxi<0) break; // no assignment;
        invmap[maxj]=maxi;
        fwdmap[maxi]=maxj;
    }

    double total_score=0;
    for (j=0;j<ylen;j++)
    {
        i=invmap[j];
        if (i>=0) total_score+=score[i+1][j+1];
    }

    /* stage 2 - swap assignment until total score cannot be improved */
    int iter;
    int oldi,oldj;
    double delta_score;
    for (iter=0; iter<getmin(xlen,ylen)*5; iter++)
    {
        //cout<<"total_score="<<total_score<<".iter="<<iter<<endl;
        //print_invmap(invmap,ylen);
        delta_score=-1;
        for (i=0;i<xlen;i++)
        {
            oldj=fwdmap[i];
            for (j=0;j<ylen;j++)
            {
                oldi=invmap[j];
                if (score[i+1][j+1]<=0 || oldi==i) continue;
                if (mm_opt==6 && (!sec2sq(i,j,secx_bond,secy_bond,fwdmap,invmap) ||
                            !sec2sq(oldi,oldj,secx_bond,secy_bond,fwdmap,invmap)))
                    continue;
                delta_score=score[i+1][j+1];
                if (oldi>=0 && oldj>=0) delta_score+=score[oldi+1][oldj+1];
                if (oldi>=0) delta_score-=score[oldi+1][j+1];
                if (oldj>=0) delta_score-=score[i+1][oldj+1];

                if (delta_score>0) // successful swap
                {
                    fwdmap[i]=j;
                    if (oldi>=0) fwdmap[oldi]=oldj;
                    invmap[j]=i;
                    if (oldj>=0) invmap[oldj]=oldi;
                    total_score+=delta_score;
                    break;
                }
            }
        }
        if (delta_score<=0) break; // cannot make further swap
    }

    /* clean up */
    delete[]fwdmap;
}

/* entry function for se
 * u_opt corresponds to option -L
 *       if u_opt==2, use d0 from Lnorm_ass for alignment
 * */
int soi_se_main(
    double **xa, double **ya, const char *seqx, const char *seqy,
    double &TM1, double &TM2, double &TM3, double &TM4, double &TM5,
    double &d0_0, double &TM_0,
    double &d0A, double &d0B, double &d0u, double &d0a, double &d0_out,
    string &seqM, string &seqxA, string &seqyA,
    double &rmsd0, int &L_ali, double &Liden,
    double &TM_ali, double &rmsd_ali, int &n_ali, int &n_ali8,
    const int xlen, const int ylen, 
    const double Lnorm_ass, const double d0_scale, const bool i_opt,
    const bool a_opt, const int u_opt, const bool d_opt, const int mol_type,
    const int outfmt_opt, int *invmap, double *dist_list,
    int **secx_bond, int **secy_bond, const int mm_opt)
{
    double D0_MIN;        //for d0
    double Lnorm;         //normalization length
    double score_d8,d0,d0_search,dcu0;//for TMscore search
    double **score;       // score for aligning a residue pair
    bool   **path;        // for dynamic programming  
    double **val;         // for dynamic programming  

    int *m1=NULL;
    int *m2=NULL;
    int i,j;
    double d;
    if (outfmt_opt<2)
    {
        m1=new int[xlen]; //alignd index in x
        m2=new int[ylen]; //alignd index in y
    }

    /***********************/
    /* allocate memory     */
    /***********************/
    NewArray(&score, xlen+1, ylen+1);
    NewArray(&path,  xlen+1, ylen+1);
    NewArray(&val,   xlen+1, ylen+1);
    //int *invmap          = new int[ylen+1];

    /* set d0 */
    parameter_set4search(xlen, ylen, D0_MIN, Lnorm,
        score_d8, d0, d0_search, dcu0); // set score_d8
    parameter_set4final(xlen, D0_MIN, Lnorm,
        d0B, d0_search, mol_type); // set d0B
    parameter_set4final(ylen, D0_MIN, Lnorm,
        d0A, d0_search, mol_type); // set d0A
    if (a_opt)
        parameter_set4final((xlen+ylen)*0.5, D0_MIN, Lnorm,
            d0a, d0_search, mol_type); // set d0a
    if (u_opt)
    {
        parameter_set4final(Lnorm_ass, D0_MIN, Lnorm,
            d0u, d0_search, mol_type); // set d0u
        if (u_opt==2)
        {
            parameter_set4search(Lnorm_ass, Lnorm_ass, D0_MIN, Lnorm,
                score_d8, d0, d0_search, dcu0); // set score_d8
        }
    }

    /* perform alignment */
    for(j=0; j<ylen; j++) invmap[j]=-1;
    double d02=d0*d0;
    double score_d82=score_d8*score_d8;
    double d2;
    for(i=0; i<xlen; i++)
    {
        for(j=0; j<ylen; j++)
        {
            d2=dist(xa[i], ya[j]);
            if (d2>score_d82) score[i+1][j+1]=0;
            else score[i+1][j+1]=1./(1+ d2/d02);
        }
    }
    if (mm_opt==6) NWDP_TM(score, path, val, xlen, ylen, -0.6, invmap);
    soi_egs(score, xlen, ylen, invmap, secx_bond, secy_bond, mm_opt);

    rmsd0=TM1=TM2=TM3=TM4=TM5=0;
    int k=0;
    n_ali=0;
    n_ali8=0;
    for(j=0; j<ylen; j++)
    {
        i=invmap[j];
        dist_list[j]=-1;
        if(i>=0)//aligned
        {
            n_ali++;
            d=sqrt(dist(&xa[i][0], &ya[j][0]));
            dist_list[j]=d;
            if (score[i+1][j+1]>0)
            {
                if (outfmt_opt<2)
                {
                    m1[k]=i;
                    m2[k]=j;
                }
                k++;
                TM2+=1/(1+(d/d0B)*(d/d0B)); // chain_1
                TM1+=1/(1+(d/d0A)*(d/d0A)); // chain_2
                if (a_opt) TM3+=1/(1+(d/d0a)*(d/d0a)); // -a
                if (u_opt) TM4+=1/(1+(d/d0u)*(d/d0u)); // -u
                if (d_opt) TM5+=1/(1+(d/d0_scale)*(d/d0_scale)); // -d
                rmsd0+=d*d;
            }
        }
    }
    n_ali8=k;
    TM2/=xlen;
    TM1/=ylen;
    TM3/=(xlen+ylen)*0.5;
    TM4/=Lnorm_ass;
    TM5/=ylen;
    if (n_ali8) rmsd0=sqrt(rmsd0/n_ali8);

    if (outfmt_opt>=2)
    {
        DeleteArray(&score, xlen+1);
        return 0;
    }

    /* extract aligned sequence */
    int ali_len=xlen+ylen;
    for (j=0;j<ylen;j++) ali_len-=(invmap[j]>=0);
    seqxA.assign(ali_len,'-');
    seqM.assign( ali_len,' ');
    seqyA.assign(ali_len,'-');

    int *fwdmap = new int [xlen+1];
    for (i=0;i<xlen;i++) fwdmap[i]=-1;
    
    for (j=0;j<ylen;j++)
    {
        seqyA[j]=seqy[j];
        i=invmap[j];
        if (i<0) continue;
        if (sqrt(dist(xa[i], ya[j]))<d0_out) seqM[j]=':';
        else seqM[j]='.';
        fwdmap[i]=j;
        seqxA[j]=seqx[i];
        Liden+=(seqxA[k]==seqyA[k]);
    }
    k=0;
    for (i=0;i<xlen;i++)
    {
        j=fwdmap[i];
        if (j>=0) continue;
        seqxA[ylen+k]=seqx[i];
        k++;
    }

    /* free memory */
    delete [] fwdmap;
    delete [] m1;
    delete [] m2;
    DeleteArray(&score, xlen+1);
    DeleteArray(&path, xlen+1);
    DeleteArray(&val, xlen+1);
    return 0; // zero for no exception
}

inline void SOI_super2score(double **xt, double **ya, const int xlen,
    const int ylen, double **score, double d0, double score_d8)
{
    int i,j;
    double d02=d0*d0;
    double score_d82=score_d8*score_d8;
    double d2;
    for (i=0; i<xlen; i++)
    {
        for(j=0; j<ylen; j++)
        {
            d2=dist(xt[i], ya[j]);
            if (d2>score_d82) score[i+1][j+1]=0;
            else score[i+1][j+1]=1./(1+ d2/d02);
        }
    }
}

//heuristic run of dynamic programing iteratively to find the best alignment
//input: initial rotation matrix t, u
//       vectors x and y, d0
//output: best alignment that maximizes the TMscore, will be stored in invmap
double SOI_iter(double **r1, double **r2, double **xtm, double **ytm,
    double **xt, double **score, bool **path, double **val, double **xa, double **ya,
    int xlen, int ylen, double t[3], double u[3][3], int *invmap0,
    int iteration_max, double local_d0_search,
    double Lnorm, double d0, double score_d8,
    int **secx_bond, int **secy_bond, const int mm_opt, const bool init_invmap=false)
{
    double rmsd; 
    int *invmap=new int[ylen+1];
    
    int iteration, i, j, k;
    double tmscore, tmscore_max, tmscore_old=0;    
    tmscore_max=-1;

    //double d01=d0+1.5;
    double d02=d0*d0;
    double score_d82=score_d8*score_d8;
    double d2;
    for (iteration=0; iteration<iteration_max; iteration++)
    {
        if (iteration==0 && init_invmap) 
            for (j=0;j<ylen;j++) invmap[j]=invmap0[j];
        else
        {
            for (j=0; j<ylen; j++) invmap[j]=-1;
            if (mm_opt==6) NWDP_TM(score, path, val, xlen, ylen, -0.6, invmap);
        }
        soi_egs(score, xlen, ylen, invmap, secx_bond, secy_bond, mm_opt);
    
        k=0;
        for (j=0; j<ylen; j++) 
        {
            i=invmap[j];
            if (i<0) continue;

            xtm[k][0]=xa[i][0];
            xtm[k][1]=xa[i][1];
            xtm[k][2]=xa[i][2];
            
            ytm[k][0]=ya[j][0];
            ytm[k][1]=ya[j][1];
            ytm[k][2]=ya[j][2];
            k++;
        }

        tmscore = TMscore8_search(r1, r2, xtm, ytm, xt, k, t, u,
            40, 8, &rmsd, local_d0_search, Lnorm, score_d8, d0);

        if (tmscore>tmscore_max)
        {
            tmscore_max=tmscore;
            for (j=0; j<ylen; j++) invmap0[j]=invmap[j];
        }
    
        if (iteration>0 && fabs(tmscore_old-tmscore)<0.000001) break;       
        tmscore_old=tmscore;
        do_rotation(xa, xt, xlen, t, u);
        SOI_super2score(xt, ya, xlen, ylen, score, d0, score_d8);
    }// for iteration
    
    delete []invmap;
    return tmscore_max;
}

void get_SOI_initial_assign(double **xk, double **yk, const int closeK_opt,
    double **score, bool **path, double **val, const int xlen, const int ylen,
    double t[3], double u[3][3], int invmap[], 
    double local_d0_search, double d0, double score_d8,
    int **secx_bond, int **secy_bond, const int mm_opt)
{
    int i,j,k;
    double **xfrag;
    double **xtran;
    double **yfrag;
    NewArray(&xfrag, closeK_opt, 3);
    NewArray(&xtran, closeK_opt, 3);
    NewArray(&yfrag, closeK_opt, 3);
    double rmsd;
    double d02=d0*d0;
    double score_d82=score_d8*score_d8;
    double d2;

    /* fill in score */
    for (i=0;i<xlen;i++)
    {
        for (k=0;k<closeK_opt;k++)
        {
            xfrag[k][0]=xk[i*closeK_opt+k][0];
            xfrag[k][1]=xk[i*closeK_opt+k][1];
            xfrag[k][2]=xk[i*closeK_opt+k][2];
        }

        for (j=0;j<ylen;j++)
        {
            for (k=0;k<closeK_opt;k++)
            {
                yfrag[k][0]=yk[j*closeK_opt+k][0];
                yfrag[k][1]=yk[j*closeK_opt+k][1];
                yfrag[k][2]=yk[j*closeK_opt+k][2];
            }
            Kabsch(xfrag, yfrag, closeK_opt, 1, &rmsd, t, u);
            do_rotation(xfrag, xtran, closeK_opt, t, u);
            
            //for (k=0; k<closeK_opt; k++)
            //{
                //d2=dist(xtran[k], yfrag[k]);
                //if (d2>score_d82) score[i+1][j+1]=0;
                //else score[i+1][j+1]=1./(1+d2/d02);
            //}
            k=closeK_opt-1;
            d2=dist(xtran[k], yfrag[k]);
            if (d2>score_d82) score[i+1][j+1]=0;
            else score[i+1][j+1]=1./(1+d2/d02);
        }
    }

    /* initial assignment */
    for (j=0;j<ylen;j++) invmap[j]=-1;
    if (mm_opt==6) NWDP_TM(score, path, val, xlen, ylen, -0.6, invmap);
    for (j=0; j<ylen;j++) i=invmap[j];
    soi_egs(score, xlen, ylen, invmap, secx_bond, secy_bond, mm_opt);

    /* clean up */
    DeleteArray(&xfrag, closeK_opt);
    DeleteArray(&xtran, closeK_opt);
    DeleteArray(&yfrag, closeK_opt);
}

void SOI_assign2super(double **r1, double **r2, double **xtm, double **ytm,
    double **xt, double **xa, double **ya,
    const int xlen, const int ylen, double t[3], double u[3][3], int invmap[], 
    double local_d0_search, double Lnorm, double d0, double score_d8)
{
    int i,j,k;
    double rmsd;
    double d02=d0*d0;
    double score_d82=score_d8*score_d8;
    double d2;

    k=0;
    for (j=0; j<ylen; j++)
    {
        i=invmap[j];
        if (i<0) continue;
        xtm[k][0]=xa[i][0];
        xtm[k][1]=xa[i][1];
        xtm[k][2]=xa[i][2];

        ytm[k][0]=ya[j][0];
        ytm[k][1]=ya[j][1];
        ytm[k][2]=ya[j][2];
        k++;
    }
    TMscore8_search(r1, r2, xtm, ytm, xt, k, t, u,
        40, 8, &rmsd, local_d0_search, Lnorm, score_d8, d0);
    do_rotation(xa, xt, xlen, t, u);
}

/* entry function for TM-align with circular permutation
 * i_opt, a_opt, u_opt, d_opt, TMcut are not implemented yet */
int SOIalign_main(double **xa, double **ya,
    double **xk, double **yk, const int closeK_opt,
    const char *seqx, const char *seqy, const char *secx, const char *secy,
    double t0[3], double u0[3][3],
    double &TM1, double &TM2, double &TM3, double &TM4, double &TM5,
    double &d0_0, double &TM_0,
    double &d0A, double &d0B, double &d0u, double &d0a, double &d0_out,
    string &seqM, string &seqxA, string &seqyA, int *invmap,
    double &rmsd0, int &L_ali, double &Liden,
    double &TM_ali, double &rmsd_ali, int &n_ali, int &n_ali8,
    const int xlen, const int ylen,
    const vector<string> sequence, const double Lnorm_ass,
    const double d0_scale, const int i_opt, const int a_opt,
    const bool u_opt, const bool d_opt, const bool fast_opt,
    const int mol_type, double *dist_list, 
    int **secx_bond, int **secy_bond, const int mm_opt)
{
    double D0_MIN;        //for d0
    double Lnorm;         //normalization length
    double score_d8,d0,d0_search,dcu0;//for TMscore search
    double t[3], u[3][3]; //Kabsch translation vector and rotation matrix
    double **score;       // Input score table for enhanced greedy search
    double **scoret;      // Transposed score table for enhanced greedy search
    bool   **path;        // for dynamic programming  
    double **val;         // for dynamic programming  
    double **xtm, **ytm;  // for TMscore search engine
    double **xt;          //for saving the superposed version of r_1 or xtm
    double **yt;          //for saving the superposed version of r_2 or ytm
    double **r1, **r2;    // for Kabsch rotation

    /***********************/
    /* allocate memory     */
    /***********************/
    int minlen = min(xlen, ylen);
    int maxlen = (xlen>ylen)?xlen:ylen;
    NewArray(&score,  xlen+1, ylen+1);
    NewArray(&scoret, ylen+1, xlen+1);
    NewArray(&path, maxlen+1, maxlen+1);
    NewArray(&val,  maxlen+1, maxlen+1);
    NewArray(&xtm, minlen, 3);
    NewArray(&ytm, minlen, 3);
    NewArray(&xt, xlen, 3);
    NewArray(&yt, ylen, 3);
    NewArray(&r1, minlen, 3);
    NewArray(&r2, minlen, 3);

    /***********************/
    /*    parameter set    */
    /***********************/
    parameter_set4search(xlen, ylen, D0_MIN, Lnorm, 
        score_d8, d0, d0_search, dcu0);
    int simplify_step    = 40; //for simplified search engine
    int score_sum_method = 8;  //for scoring method, whether only sum over pairs with dis<score_d8

    int i,j;
    int *fwdmap0         = new int[xlen+1];
    int *invmap0         = new int[ylen+1];
    
    double TMmax=-1, TM=-1;
    for(i=0; i<xlen; i++) fwdmap0[i]=-1;
    for(j=0; j<ylen; j++) invmap0[j]=-1;
    double local_d0_search = d0_search;
    int iteration_max=(fast_opt)?2:30;
    //if (mm_opt==6) iteration_max=1;

    /*************************************************************/
    /* initial alignment with sequence order dependent alignment */
    /*************************************************************/
    CPalign_main(
        xa, ya, seqx, seqy, secx, secy,
        t0, u0, TM1, TM2, TM3, TM4, TM5,
        d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out, seqM, seqxA, seqyA,
        rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
        xlen, ylen, sequence, Lnorm_ass, d0_scale,
        i_opt, a_opt, u_opt, d_opt, fast_opt,
        mol_type,-1);
    if (mm_opt==6)
    {
        i=0;
        j=0;
        for (int r=0;r<seqxA.size();r++)
        {
            if (seqxA[r]=='*') // circular permutation point
            {
                for (int jj=0;jj<j;jj++) if (invmap0[jj]>=0)
                    invmap0[jj]+=xlen - i;
                i=0;
                continue;
            }
            if (seqyA[r]!='-')
            {
                if (seqxA[r]!='-') invmap0[j]=i;
                j++;
            }
            if (seqxA[r]!='-') i++;
        }
        for (j=0;j<ylen;j++)
        {
            i=invmap0[j];
            if (i>=0) fwdmap0[i]=j;
        }
    }
    do_rotation(xa, xt, xlen, t0, u0);
    SOI_super2score(xt, ya, xlen, ylen, score, d0, score_d8);
    for (i=0;i<xlen;i++) for (j=0;j<ylen;j++) scoret[j+1][i+1]=score[i+1][j+1];
    TMmax=SOI_iter(r1, r2, xtm, ytm, xt, score, path, val, xa, ya,
        xlen, ylen, t0, u0, invmap0, iteration_max,
        local_d0_search, Lnorm, d0, score_d8, secx_bond, secy_bond, mm_opt, true);
    TM   =SOI_iter(r2, r1, ytm, xtm, yt,scoret, path, val, ya, xa,
        ylen, xlen, t0, u0, fwdmap0, iteration_max,
        local_d0_search, Lnorm, d0, score_d8, secy_bond, secx_bond, mm_opt, true);
    //cout<<"TM2="<<TM2<<"\tTM1="<<TM1<<"\tTMmax="<<TMmax<<"\tTM="<<TM<<endl;
    if (TM>TMmax)
    {
        TMmax = TM;
        for (j=0; j<ylen; j++) invmap0[j]=-1;
        for (i=0; i<xlen; i++) 
        {
            j=fwdmap0[i];
            if (j>=0) invmap0[j]=i;
        }
    }
    
    /***************************************************************/
    /* initial alignment with sequence order independent alignment */
    /***************************************************************/
    if (closeK_opt>=3)
    {
        get_SOI_initial_assign(xk, yk, closeK_opt, score, path, val,
            xlen, ylen, t, u, invmap, local_d0_search, d0, score_d8,
            secx_bond, secy_bond, mm_opt);
        for (i=0;i<xlen;i++) for (j=0;j<ylen;j++) scoret[j+1][i+1]=score[i+1][j+1];

        SOI_assign2super(r1, r2, xtm, ytm, xt, xa, ya,
            xlen, ylen, t, u, invmap, local_d0_search, Lnorm, d0, score_d8);
        TM=SOI_iter(r1, r2, xtm, ytm, xt, score, path, val, xa, ya,
            xlen, ylen, t, u, invmap, iteration_max,
            local_d0_search, Lnorm, d0, score_d8, secx_bond, secy_bond, mm_opt);
        if (TM>TMmax)
        {
            TMmax = TM;
            for (j = 0; j<ylen; j++) invmap0[j] = invmap[j];
        }

        for (i=0;i<xlen;i++) fwdmap0[i]=-1;
        if (mm_opt==6) NWDP_TM(scoret, path, val, ylen, xlen, -0.6, fwdmap0);
        soi_egs(scoret, ylen, xlen, fwdmap0, secy_bond, secx_bond, mm_opt);
        SOI_assign2super(r2, r1, ytm, xtm, yt, ya, xa,
            ylen, xlen, t, u, fwdmap0, local_d0_search, Lnorm, d0, score_d8);
        TM=SOI_iter(r2, r1, ytm, xtm, yt, scoret, path, val, ya, xa, ylen, xlen, t, u,
            fwdmap0, iteration_max, local_d0_search, Lnorm, d0, score_d8,secy_bond, secx_bond, mm_opt);
        if (TM>TMmax)
        {
            TMmax = TM;
            for (j=0; j<ylen; j++) invmap0[j]=-1;
            for (i=0; i<xlen; i++) 
            {
                j=fwdmap0[i];
                if (j>=0) invmap0[j]=i;
            }
        }
    }

    //*******************************************************************//
    //    The alignment will not be changed any more in the following    //
    //*******************************************************************//
    //check if the initial alignment is generated appropriately
    bool flag=false;
    for (i=0; i<xlen; i++) fwdmap0[i]=-1;
    for (j=0; j<ylen; j++)
    {
        i=invmap0[j];
        invmap[j]=i;
        if (i>=0)
        {
            fwdmap0[i]=j;
            flag=true;
        }
    }
    if(!flag)
    {
        cout << "There is no alignment between the two structures! "
             << "Program stop with no result!" << endl;
        TM1=TM2=TM3=TM4=TM5=0;
        return 1;
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
    
    double rmsd;
    simplify_step=1;
    score_sum_method=0;
    double Lnorm_0=ylen;

    //select pairs with dis<d8 for final TMscore computation and output alignment
    int k=0;
    int *m1, *m2;
    double d;
    m1=new int[xlen]; //alignd index in x
    m2=new int[ylen]; //alignd index in y
    copy_t_u(t, u, t0, u0);
    
    //****************************************//
    //              Final TMscore 1           //
    //****************************************//

    do_rotation(xa, xt, xlen, t, u);
    k=0;
    n_ali=0;
    for (i=0; i<xlen; i++)
    {
        j=fwdmap0[i];
        if(j>=0)//aligned
        {
            n_ali++;
            d=sqrt(dist(&xt[i][0], &ya[j][0]));
            if (d <= score_d8)
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
            else fwdmap0[i]=-1;
        }
    }
    n_ali8=k;

    Kabsch(r1, r2, n_ali8, 0, &rmsd0, t, u);// rmsd0 is used for final output, only recalculate rmsd0, not t & u
    rmsd0 = sqrt(rmsd0 / n_ali8);
    
    //normalized by length of structure A
    parameter_set4final(xlen+0.0, D0_MIN, Lnorm, d0, d0_search, mol_type);
    d0B=d0;
    local_d0_search = d0_search;
    TM2 = TMscore8_search(r1, r2, xtm, ytm, xt, n_ali8, t, u, simplify_step,
        score_sum_method, &rmsd, local_d0_search, Lnorm, score_d8, d0);

    //****************************************//
    //              Final TMscore 2           //
    //****************************************//
    
    do_rotation(xa, xt, xlen, t0, u0);
    k=0;
    for (j=0; j<ylen; j++)
    {
        i=invmap0[j];
        if(i>=0)//aligned
        {
            d=sqrt(dist(&xt[i][0], &ya[j][0]));
            if (d <= score_d8)
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
            else invmap[j]=invmap0[j]=-1;
        }
    }

    //normalized by length of structure B
    parameter_set4final(Lnorm_0, D0_MIN, Lnorm, d0, d0_search, mol_type);
    d0A=d0;
    d0_0=d0A;
    local_d0_search = d0_search;
    TM1 = TMscore8_search(r1, r2, xtm, ytm, xt, n_ali8, t0, u0, simplify_step,
        score_sum_method, &rmsd, local_d0_search, Lnorm, score_d8, d0);
    TM_0 = TM1;

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
        local_d0_search = d0_search;
        TM5 = TMscore8_search(r1, r2, xtm, ytm, xt, n_ali8, t0, u0,
            simplify_step, score_sum_method, &rmsd, local_d0_search, Lnorm,
            score_d8, d0);
        TM_0=TM5;
    }

    /* derive alignment from superposition */
    int ali_len=xlen+ylen;
    for (j=0;j<ylen;j++) ali_len-=(invmap0[j]>=0);
    seqxA.assign(ali_len,'-');
    seqM.assign( ali_len,' ');
    seqyA.assign(ali_len,'-');
    
    //do_rotation(xa, xt, xlen, t, u);
    do_rotation(xa, xt, xlen, t0, u0);

    Liden=0;
    //double SO=0;
    for (j=0;j<ylen;j++)
    {
        seqyA[j]=seqy[j];
        i=invmap0[j];
        dist_list[j]=-1;
        if (i<0) continue;
        d=sqrt(dist(xt[i], ya[j]));
        if (d<d0_out) seqM[j]=':';
        else seqM[j]='.';
        dist_list[j]=d;
        //SO+=(d<3.5);
        seqxA[j]=seqx[i];
        Liden+=(seqx[i]==seqy[j]);
    }
    //SO/=getmin(xlen,ylen);
    k=0;
    for (i=0;i<xlen;i++)
    {
        j=fwdmap0[i];
        if (j>=0) continue;
        seqxA[ylen+k]=seqx[i];
        k++;
    }
    //cout<<n_ali8<<'\t'
        //<<rmsd0<<'\t'
        //<<100.*SO<<endl;


    /* clean up */
    DeleteArray(&score, xlen+1);
    DeleteArray(&scoret,ylen+1);
    DeleteArray(&path,maxlen+1);
    DeleteArray(&val, maxlen+1);
    DeleteArray(&xtm, minlen);
    DeleteArray(&ytm, minlen);
    DeleteArray(&xt,xlen);
    DeleteArray(&yt,ylen);
    DeleteArray(&r1, minlen);
    DeleteArray(&r2, minlen);
    delete[]invmap0;
    delete[]fwdmap0;
    delete[]m1;
    delete[]m2;
    return 0;
}
#endif
