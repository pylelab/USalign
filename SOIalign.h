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

void soi_egs(double **score, const int xlen, const int ylen, int *invmap)
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
                if (invmap[j]==-1 && score[i][j]>max_score)
                {
                    maxi=i;
                    maxj=j;
                    max_score=score[i][j];
                }
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
        if (i>=0) total_score+=score[i][j];
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
                if (score[i][j]<=0 || oldi==i) continue;
                delta_score=score[i][j];
                if (oldi>=0 && oldj>=0) delta_score+=score[oldi][oldj];
                if (oldi>=0) delta_score-=score[oldi][j];
                if (oldj>=0) delta_score-=score[i][oldj];

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
    const int outfmt_opt, int *invmap)
{
    double D0_MIN;        //for d0
    double Lnorm;         //normalization length
    double score_d8,d0,d0_search,dcu0;//for TMscore search
    double **score;       // score for aligning a residue pair

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
    NewArray(&score, xlen, ylen);
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
    for(i=0; i<xlen; i++)
    {
        for(j=0; j<ylen; j++)
        {
            d=sqrt(dist(xa[i], ya[j]));
            if (d>score_d8) score[i][j]=0;
            else score[i][j]=1./(1+ d/d02);
        }
    }
    soi_egs(score, xlen, ylen, invmap);

    rmsd0=TM1=TM2=TM3=TM4=TM5=0;
    int k=0;
    n_ali=0;
    n_ali8=0;
    for(j=0; j<ylen; j++)
    {
        i=invmap[j];
        if(i>=0)//aligned
        {
            n_ali++;
            d=sqrt(dist(&xa[i][0], &ya[j][0]));
            if (score[i][j]>0)
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
        DeleteArray(&score, xlen);
        return 0;
    }

    /* extract aligned sequence */
    seqxA.assign(n_ali8,'-');
    seqM.assign( n_ali8,' ');
    seqyA.assign(n_ali8,'-');
    
    k=0;
    for (j=0;j<ylen;j++)
    {
        i=invmap[j];
        if (i<0) continue;
        if (sqrt(dist(xa[i], ya[j]))<d0_out) seqM[k]=':';
        seqxA[k]=seqx[i];
        seqyA[k]=seqy[j];
        Liden+=(seqxA[k]==seqyA[k]);
        k++;
    }


    /* free memory */
    //delete [] invmap;
    delete [] m1;
    delete [] m2;
    DeleteArray(&score, xlen);
    return 0; // zero for no exception
}

/* entry function for TM-align with circular permutation
 * i_opt, a_opt, u_opt, d_opt, TMcut are not implemented yet */
int SOIalign_main(double **xa, double **ya,
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
    const int mol_type)
{
    CPalign_main(xa, ya, seqx, seqy, secx, secy,
        t0, u0, TM1, TM2, TM3, TM4, TM5,
        d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
        seqM, seqxA, seqyA,
        rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
        xlen, ylen, sequence, Lnorm_ass, d0_scale,
        i_opt, a_opt, u_opt, d_opt, fast_opt,
        mol_type,-1);

    double **xa_cp;   // coordinates
    NewArray(&xa_cp, xlen, 3);
    for (int r=0; r<xlen; r++) transform(t0, u0, xa[r], xa_cp[r]);
    
    /* TODO: iterate between EGS and superimposition */

    seqxA.clear();
    seqM.clear();
    seqyA.clear();
    soi_se_main(
        xa_cp, ya, seqx, seqy, TM1, TM2, TM3, TM4, TM5,
        d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
        seqM, seqxA, seqyA,
        rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
        xlen, ylen, Lnorm_ass, d0_scale,
        i_opt, a_opt, u_opt, d_opt,
        mol_type, 0, invmap);

    /* clean up */
    DeleteArray(&xa_cp,xlen);
    return 0;
}
#endif
