/* Functions for the core TMalign algorithm, including the entry function
 * flexalign_main */
#ifndef flexalign_h
#define flexalign_h 1

#include "TMalign.h"

void t_u2tu(double t0[3],double u0[3][3], vector<double> &tu_tmp)
{
    int i,j,k;
    for (i=0;i<3;i++) tu_tmp[i]=t0[i];
    k=3;
    for (i=0;i<3;i++) for (j=0;j<3;j++)
    {
        tu_tmp[k]=u0[i][j];
        k++;
    }
}

void tu2t_u(vector<double> tu_tmp, double t0[3],double u0[3][3])
{
    int i,j,k;
    for (i=0;i<3;i++) t0[i]=tu_tmp[i];
    k=3;
    for (i=0;i<3;i++) for (j=0;j<3;j++)
    {
        u0[i][j]=tu_tmp[k];
        k++;
    }
}

void aln2invmap(const string &seqxA, const string &seqyA, int *invmap)
{
    int i,j,r;
    int ylen=0;
    for (r=0;r<seqyA.size();r++) ylen+=seqyA[r]!='-';
    for(j=0; j<ylen; j++) invmap[j]=-1;

    i=j=-1;
    for (r=0;r<seqxA.size();r++)
    {
        i+=seqxA[r]!='-';
        j+=seqyA[r]!='-';
        if (seqxA[r]!='-' && seqyA[r]!='-') invmap[j]=i;
    }
}

/* Entry function for TM-align. Return TM-score calculation status:
 * 0   - full TM-score calculation 
 * 1   - terminated due to exception
 * 2-7 - pre-terminated due to low TM-score */
int flexalign_main(double **xa, double **ya,
    const char *seqx, const char *seqy, const char *secx, const char *secy,
    double t0[3], double u0[3][3], vector<vector<double> >&tu_vec,
    double &TM1, double &TM2, double &TM3, double &TM4, double &TM5,
    double &d0_0, double &TM_0,
    double &d0A, double &d0B, double &d0u, double &d0a, double &d0_out,
    string &seqM, string &seqxA, string &seqyA,
    double &rmsd0, int &L_ali, double &Liden,
    double &TM_ali, double &rmsd_ali, int &n_ali, int &n_ali8,
    const int xlen, const int ylen,
    const vector<string> sequence, const double Lnorm_ass,
    const double d0_scale, const int i_opt, const int a_opt,
    const bool u_opt, const bool d_opt, const bool fast_opt,
    const int mol_type, const int hinge_opt)
{
    TMalign_main(xa, ya, seqx, seqy, secx, secy, t0, u0,
        TM1, TM2, TM3, TM4, TM5, d0_0, TM_0,
        d0A, d0B, d0u, d0a, d0_out, seqM, seqxA, seqyA,
        rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
        xlen, ylen, sequence, Lnorm_ass,
        d0_scale, i_opt, a_opt, u_opt, d_opt, fast_opt, mol_type);
    
    vector<double> tu_tmp(12,0);
    t_u2tu(t0,u0,tu_tmp);
    tu_vec.push_back(tu_tmp);
    
    int* invmap=new int[ylen+1];
    aln2invmap(seqxA, seqyA, invmap);
    double **xt;
    NewArray(&xt, xlen, 3);
    for (int r=0;r<seqM.size();r++) if (seqM[r]!=' ') seqM[r]='0';

    int minlen = min(xlen, ylen);
    int hinge;
    for (hinge=0;hinge<hinge_opt;hinge++)
    {
        if (minlen-n_ali8<5) break;
        int xlen_h=xlen - n_ali8;
        int ylen_h=ylen - n_ali8;
        char *seqx_h = new char[xlen_h + 1];
        char *seqy_h = new char[ylen_h + 1];
        char *secx_h = new char[xlen_h + 1];
        char *secy_h = new char[ylen_h + 1];
        seqx_h[xlen_h]=seqy_h[ylen_h]=0;
        secx_h[xlen_h]=secy_h[ylen_h]=0;
        double **xa_h, **ya_h;
        NewArray(&xa_h, xlen_h, 3);
        NewArray(&ya_h, ylen_h, 3);
        vector<int> r1toi(xlen_h,0);
        vector<int> r2toj(ylen_h,0);

        int i,j,r;
        int r1,r2;
        i=j=-1;
        r1=r2=0;
        for (r=0;r<seqxA.size();r++)
        {
            i+=(seqxA[r]!='-');
            j+=(seqyA[r]!='-');
            if (seqyA[r]=='-')
            {
                seqx_h[r1]=seqx[i];
                secx_h[r1]=secx[i];
                xa_h[r1][0]=xa[i][0];
                xa_h[r1][1]=xa[i][1];
                xa_h[r1][2]=xa[i][2];
                r1toi[r1]=i;
                r1++;
            }
            if (seqxA[r]=='-')
            {
                seqy_h[r2]=seqx[j];
                secy_h[r2]=secx[j];
                ya_h[r2][0]=ya[j][0];
                ya_h[r2][1]=ya[j][1];
                ya_h[r2][2]=ya[j][2];
                r2toj[r2]=j;
                r2++;
            }
        }
        
        double TM1_h, TM2_h;
        double TM3_h, TM4_h, TM5_h;     // for a_opt, u_opt, d_opt
        double d0_0_h, TM_0_h;
        double d0A_h, d0B_h, d0u_h, d0a_h;
        double d0_out_h=5.0;
        string seqM_h, seqxA_h, seqyA_h;// for output alignment
        double rmsd0_h = 0.0;
        int L_ali_h=0;                // Aligned length in standard_TMscore
        double Liden_h=0;
        double TM_ali_h, rmsd_ali_h;  // TMscore and rmsd in standard_TMscore
        int n_ali_h=0;
        int n_ali8_h=0;

        TMalign_main(xa_h, ya_h, seqx_h, seqy_h, secx_h, secy_h, t0, u0,
            TM1_h, TM2_h, TM3_h, TM4_h, TM5_h, d0_0_h, TM_0_h,
            d0A_h, d0B_h, d0u_h, d0a_h, d0_out_h, seqM_h, seqxA_h, seqyA_h,
            rmsd0_h, L_ali_h, Liden_h, TM_ali_h, rmsd_ali_h, n_ali_h, n_ali8_h,
            xlen_h, ylen_h, sequence, Lnorm_ass,
            d0_scale, i_opt, a_opt, u_opt, d_opt, fast_opt, mol_type);
        
        t_u2tu(t0,u0,tu_tmp);
        tu_vec.push_back(tu_tmp);
        do_rotation(xa, xt, xlen, t0, u0);
        
        TM1_h=TM1;
        TM2_h=TM2;
        TM3_h=TM3;
        TM4_h=TM4;
        TM5_h=TM5;
        seqM_h=seqM;
        seqxA_h=seqxA;
        seqyA_h=seqyA;
        rmsd0_h=rmsd0;
        n_ali_h=n_ali;
        n_ali8_h=n_ali8;
        int* invmap_h=new int[ylen+1];
        for (j=0;j<ylen+1;j++) invmap_h[j]=invmap[j];
        se_main(xt, ya, seqx, seqy, TM1_h, TM2_h, TM3_h, TM4_h, TM5_h, d0_0, TM_0,
            d0A, d0B, d0u, d0a, d0_out, seqM_h, seqxA_h, seqyA_h,
            rmsd0_h, L_ali, Liden, TM_ali, rmsd_ali, n_ali_h, n_ali8_h,
            xlen, ylen, sequence, Lnorm_ass, d0_scale, i_opt,
            a_opt, u_opt, d_opt, mol_type, 0, invmap, hinge+1);
        int new_ali=0;
        for (r=0;r<seqM_h.size();r++) new_ali+=(seqM_h[r]==hinge+'1');
        if (new_ali>=5)
        {
            TM1=TM1_h;
            TM2=TM2_h;
            TM3=TM3_h;
            TM4=TM4_h;
            TM5=TM5_h;
            seqM=seqM_h;
            seqxA=seqxA_h;
            seqyA=seqyA_h;
            rmsd0=rmsd0_h;
            n_ali=n_ali_h;
            n_ali8=n_ali8_h;
        }
        
        /* clean up */
        delete [] invmap_h;
        DeleteArray(&xa_h, xlen_h);
        DeleteArray(&ya_h, ylen_h);
        r1toi.clear();
        r2toj.clear();
        seqM_h.clear();
        seqxA_h.clear();
        seqyA_h.clear();
        delete [] seqx_h;
        delete [] secx_h;
        delete [] seqy_h;
        delete [] secy_h;
        if (new_ali<5) break;
    }
    
    /* re-derive alignment based on tu_vec */
    vector<char> seqM_char(ylen,' ');
    vector<double> di_vec(ylen,-1);
    int i,j,r;
    double d;
    for (hinge=tu_vec.size()-1;hinge>=0;hinge--)
    {
        tu2t_u(tu_vec[hinge],t0,u0);
        //output_rotation_matrix("-", t0, u0); 
        do_rotation(xa, xt, xlen, t0, u0);
        for (j=0;j<ylen;j++)
        {
            i=invmap[j];
            if (i<0) continue;
            d=sqrt(dist(xt[i], ya[j]));
            if (di_vec[j]<0 || d<=di_vec[j])
            {
                di_vec[j]=d;
                seqM_char[j]=hinge+'0';
            }
        }
    }
    j=-1;
    for (r=0;r<seqM.size();r++)
    {
        if (seqyA[r]=='-') continue;
        j++;
        seqM[r]=seqM_char[j];
    }
    
    rmsd0=TM1=TM2=TM3=TM4=TM5=0;
    for(j=0; j<ylen; j++)
    {
        i=invmap[j];
        if(i<0) continue;
        {
            d=di_vec[j];
            TM2+=1/(1+(d/d0B)*(d/d0B)); // chain_1
            TM1+=1/(1+(d/d0A)*(d/d0A)); // chain_2
            if (a_opt) TM3+=1/(1+(d/d0a)*(d/d0a)); // -a
            if (u_opt) TM4+=1/(1+(d/d0u)*(d/d0u)); // -u
            if (d_opt) TM5+=1/(1+(d/d0_scale)*(d/d0_scale)); // -d
            rmsd0+=d*d;
        }
    }
    TM2/=xlen;
    TM1/=ylen;
    TM3/=(xlen+ylen)*0.5;
    TM4/=Lnorm_ass;
    TM5/=ylen;
    if (n_ali8) rmsd0=sqrt(rmsd0/n_ali8);


    /* clean up */
    seqM_char.clear();
    di_vec.clear();
    DeleteArray(&xt, xlen);
    delete[] invmap;
    return 0; // zero for no exception
}

#endif
