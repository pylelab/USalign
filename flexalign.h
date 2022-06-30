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
    vector<double> tu_tmp(12,0);
    int round2=tu_vec.size();
    if (round2==0)
    {
        TMalign_main(xa, ya, seqx, seqy, secx, secy, t0, u0,
            TM1, TM2, TM3, TM4, TM5, d0_0, TM_0,
            d0A, d0B, d0u, d0a, d0_out, seqM, seqxA, seqyA,
            rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
            xlen, ylen, sequence, Lnorm_ass,
            d0_scale, i_opt, a_opt, u_opt, d_opt, fast_opt, mol_type);
    
        t_u2tu(t0,u0,tu_tmp);
        tu_vec.push_back(tu_tmp);
    }
    
    int i,j,r;
    int* invmap=new int[ylen+1];
    for (j=0;j<ylen+1;j++) invmap[j]=-1;
    double **xt;
    NewArray(&xt, xlen, 3);
    do_rotation(xa, xt, xlen, t0, u0);

    TM1= TM2= TM3= TM4= TM5=rmsd0=0;
    seqM="";
    seqxA="";
    seqyA="";
    n_ali=n_ali8=0;
    se_main(xt, ya, seqx, seqy, TM1, TM2, TM3, TM4, TM5, d0_0, TM_0,
        d0A, d0B, d0u, d0a, d0_out, seqM, seqxA, seqyA,
        rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
        xlen, ylen, sequence, Lnorm_ass, d0_scale, i_opt,
        a_opt, u_opt, d_opt, mol_type, 0, invmap, 1);
    if (round2)
    {
        /* aligned structure A vs unaligned structure B */
        int xlen_h=n_ali8;
        int ylen_h=ylen - n_ali8;
        char *seqx_h = new char[xlen + 1];
        char *seqy_h = new char[ylen + 1];
        char *secx_h = new char[xlen + 1];
        char *secy_h = new char[ylen + 1];
        seqx_h[xlen]=seqy_h[ylen]=0;
        secx_h[xlen]=secy_h[ylen]=0;
        double **xa_h, **ya_h;
        NewArray(&xa_h, xlen, 3);
        NewArray(&ya_h, ylen, 3);

        int r1,r2;
        i=j=-1;
        r1=r2=0;
        for (r=0;r<seqxA.size();r++)
        {
            i+=(seqxA[r]!='-');
            j+=(seqyA[r]!='-');
            if (seqxA[r]!='-' && seqyA[r]!='-')
            {
                seqx_h[r1]=seqx[i];
                secx_h[r1]=secx[i];
                xa_h[r1][0]=xa[i][0];
                xa_h[r1][1]=xa[i][1];
                xa_h[r1][2]=xa[i][2];
                r1++;
            }
            if (seqxA[r]=='-')
            {
                seqy_h[r2]=seqx[j];
                secy_h[r2]=secx[j];
                ya_h[r2][0]=ya[j][0];
                ya_h[r2][1]=ya[j][1];
                ya_h[r2][2]=ya[j][2];
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
        
        do_rotation(xa, xt, xlen, t0, u0);
        t_u2tu(t0,u0,tu_vec[0]);
        
        int* invmap_h=new int[ylen+1];
        for (j=0;j<ylen+1;j++) invmap_h[j]=-1;
        TM1_h= TM2_h= TM3_h= TM4_h= TM5_h=rmsd0_h=0;
        seqM_h="";
        seqxA_h="";
        seqyA_h="";
        n_ali_h=n_ali8_h=0;
        se_main(xt, ya, seqx, seqy, TM1_h, TM2_h, TM3_h, TM4_h, TM5_h, d0_0,
            TM_0, d0A, d0B, d0u, d0a, d0_out, seqM_h, seqxA_h, seqyA_h,
            rmsd0_h, L_ali, Liden, TM_ali, rmsd_ali, n_ali_h, n_ali8_h,
            xlen, ylen, sequence, Lnorm_ass, d0_scale, i_opt,
            a_opt, u_opt, d_opt, mol_type, 0, invmap_h, 1);
        
        /* unaligned structure A vs aligned structure B */
        xlen_h=xlen - n_ali8;
        ylen_h=n_ali8;

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
                r1++;
            }
            if (seqxA[r]!='-' && seqyA[r]!='-')
            {
                seqy_h[r2]=seqx[j];
                secy_h[r2]=secx[j];
                ya_h[r2][0]=ya[j][0];
                ya_h[r2][1]=ya[j][1];
                ya_h[r2][2]=ya[j][2];
                r2++;
            }
        }
        
        d0_out_h=5.0;
        L_ali_h=Liden_h=0;
        TM1= TM2= TM3= TM4= TM5=rmsd0=0;
        seqM="";
        seqxA="";
        seqyA="";
        n_ali=n_ali8=0;

        TMalign_main(xa_h, ya_h, seqx_h, seqy_h, secx_h, secy_h, t0, u0,
            TM1, TM2, TM3, TM4, TM5, d0_0_h, TM_0_h,
            d0A_h, d0B_h, d0u_h, d0a_h, d0_out_h, seqM, seqxA, seqyA,
            rmsd0, L_ali_h, Liden_h, TM_ali_h, rmsd_ali_h, n_ali, n_ali8,
            xlen_h, ylen_h, sequence, Lnorm_ass,
            d0_scale, i_opt, a_opt, u_opt, d_opt, fast_opt, mol_type);
        
        do_rotation(xa, xt, xlen, t0, u0);
        
        for (j=0;j<ylen+1;j++) invmap[j]=-1;
        TM1= TM2= TM3= TM4= TM5=rmsd0=0;
        seqM="";
        seqxA="";
        seqyA="";
        n_ali=n_ali8=0;
        se_main(xt, ya, seqx, seqy, TM1, TM2, TM3, TM4, TM5, d0_0,
            TM_0, d0A, d0B, d0u, d0a, d0_out, seqM, seqxA, seqyA,
            rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
            xlen, ylen, sequence, Lnorm_ass, d0_scale, i_opt,
            a_opt, u_opt, d_opt, mol_type, 0, invmap, 1);

        double TM_h=(TM1_h>TM2_h)?TM1_h:TM2_h;
        double TM  =(TM1  >TM2  )?TM1  :TM2  ;
        if (TM_h>TM)
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
            for (j=0;j<ylen+1;j++) invmap[j]=invmap_h[j];
        }
        else t_u2tu(t0,u0,tu_vec[0]);
        
        /* clean up */
        delete [] invmap_h;
        DeleteArray(&xa_h, xlen);
        DeleteArray(&ya_h, ylen);
        seqM_h.clear();
        seqxA_h.clear();
        seqyA_h.clear();
        delete [] seqx_h;
        delete [] secx_h;
        delete [] seqy_h;
        delete [] secy_h;
    }
    for (r=0;r<seqM.size();r++) if (seqM[r]=='1') seqM[r]='0';

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
            a_opt, u_opt, d_opt, mol_type, 0, invmap_h, hinge+1);
        int new_ali=0;
        for (r=0;r<seqM_h.size();r++) new_ali+=(seqM_h[r]==hinge+'1');
        if (n_ali8_h - n_ali8<5) new_ali=0;
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
            t_u2tu(t0,u0,tu_tmp);
            tu_vec.push_back(tu_tmp);
            for (j=0;j<ylen+1;j++) invmap[j]=invmap_h[j];
            //cout<<">hinge="<<hinge<<'\n'
                //<<seqxA<<'\n'<<seqM<<'\n'<<seqyA<<endl;
            //for (j=0;j<ylen;j++) if ((i=invmap[j])>=0) cout<<"("<<i<<","<<j<<")";
            //cout<<endl;
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

    if (tu_vec.size()<=1)
    {
        DeleteArray(&xt, xlen);
        delete[] invmap;
        return tu_vec.size();
    }
    
    /* re-derive alignment based on tu_vec */
    vector<char> seqM_char(ylen,' ');
    vector<double> di_vec(ylen,-1);
    double d;
    for (hinge=tu_vec.size()-1;hinge>=0;hinge--)
    {
        tu2t_u(tu_vec[hinge],t0,u0);
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

    /* smooth out AFP assignment: remove singleton insert */
    for (hinge=tu_vec.size()-1;hinge>=0;hinge--)
    {
        j=-1;
        for (r=0;r<seqM.size();r++)
        {
            if (seqyA[r]=='-') continue;
            j++;
            if (seqM_char[j]!=hinge+'0') continue;
            if (r<seqM.size()-1 && (seqM[r+1]==hinge+'0' || seqM[r+1]==' '))
                continue;
            if (r>0 && (seqM[r-1]==hinge+'0' || seqM[r-1]==' ')) continue;
            if (r<seqM.size()-1 && r>0 && seqM[r-1]!=seqM[r+1]) continue;
            if (r>0) seqM[r]=seqM_char[j]=seqM[r-1];
            else     seqM[r]=seqM_char[j]=seqM[r+1];
        }
    }
    /* smooth out AFP assignment: remove singleton at the end of fragment */
    char left_hinge=' ';
    char right_hinge=' ';
    for (hinge=tu_vec.size()-1;hinge>=0;hinge--)
    {
        j=-1;
        for (r=0;r<seqM.size();r++)
        {
            if (seqyA[r]=='-') continue;
            j++;
            if (seqM[r]!=hinge+'0') continue;
            if (r>0 && seqM[r-1]==' ' && r<seqM.size()-1 && seqM[r+1]==' ')
                continue;
            
            left_hinge=' ';
            for (i=r-1;i>=0;i--)
            {
                if (seqM[i]==' ') continue;
                left_hinge=seqM[i];
                break;
            }
            if (left_hinge==hinge+'0') continue;
            
            right_hinge=' ';
            for (i=r+1;i<seqM.size();i++)
            {
                if (seqM[i]==' ') continue;
                right_hinge=seqM[i];
                break;
            }
            if (right_hinge==hinge+'0') continue;
            if (left_hinge!=right_hinge && left_hinge!=' ' && right_hinge!=' ')
                continue;
            
            if     (right_hinge!=' ') seqM[r]=seqM_char[j]=right_hinge;
            else if (left_hinge!=' ') seqM[r]=seqM_char[j]=left_hinge;
        }
    }
    /* smooth out AFP assignment: remove dimer insert */
    for (hinge=tu_vec.size()-1;hinge>=0;hinge--)
    {
        j=-1;
        for (r=0;r<seqM.size()-1;r++)
        {
            if (seqyA[r]=='-') continue;
            j++;
            if (seqM[r]  !=hinge+'0'|| seqM[r+1]!=hinge+'0') continue;
            
            if (r<seqM.size()-2 && (seqM[r+2]==' ' || seqM[r+2]==hinge+'0'))
                continue;
            if (r>0 && (seqM[r-1]==' ' || seqM[r-1]==hinge+'0')) continue;
            if (r<seqM.size()-2 && r>0 && seqM[r-1]!=seqM[r+2]) continue;

            if (r>0) seqM[r]=seqM_char[j]=seqM[r+1]=seqM_char[j+1]=seqM[r-1];
            else     seqM[r]=seqM_char[j]=seqM[r+1]=seqM_char[j+1]=seqM[r+2];
        }
    }
    /* smooth out AFP assignment: remove disconnected singleton */
    int i1,i2;
    for (hinge=tu_vec.size()-1;hinge>=0;hinge--)
    {
        j=-1;
        for (r=0;r<seqM.size();r++)
        {
            if (seqyA[r]=='-') continue;
            j++;
            if (seqM[r]!=hinge+'0') continue;
            
            left_hinge=' ';
            for (i=r-1;i>=0;i--)
            {
                if (seqM[i]==' ') continue;
                left_hinge=seqM[i];
                i1=(r-i);
                break;
            }
            if (left_hinge==hinge+'0') continue;
            
            right_hinge=' ';
            for (i=r+1;i<seqM.size();i++)
            {
                if (seqM[i]==' ') continue;
                right_hinge=seqM[i];
                i2=(i-r);
                break;
            }
            if (right_hinge==hinge+'0') continue;
            
            if (right_hinge==' ') seqM[r]=seqM_char[j]=left_hinge;
            else if (left_hinge==' ') seqM[r]=seqM_char[j]=right_hinge;
            else
            {
                if (i1<i2) seqM[r]=seqM_char[j]=left_hinge;
                else       seqM[r]=seqM_char[j]=right_hinge;
            }
        }
    }
    
    /* recalculate all scores */
    for (hinge=tu_vec.size()-1;hinge>=0;hinge--)
    {
        tu2t_u(tu_vec[hinge],t0,u0);
        do_rotation(xa, xt, xlen, t0, u0);
        for (j=0;j<ylen;j++)
        {
            i=invmap[j];
            if (i<0) continue;
            if (seqM_char[j]!=hinge+'0') continue;
            d=sqrt(dist(xt[i], ya[j]));
            if (di_vec[j]<0 || d<=di_vec[j])
            {
                di_vec[j]=d;
                seqM_char[j]=hinge+'0';
            }
        }
    }
    rmsd0=TM1=TM2=TM3=TM4=TM5=0;
    Liden=0;
    for (r=0;r<seqM.size();r++) if (seqM[r]!=' ') Liden+=seqxA[r]==seqyA[r];
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
    for (hinge=tu_vec.size()-1;hinge>0;hinge--)
    {
        int afp_len=0;
        for (r=0;r<seqM.size();r++) afp_len+=seqM[r]==hinge+'0';
        if (afp_len) break;
        tu_vec.pop_back(); // remove unnecessary afp
    }

    /* clean up */
    seqM_char.clear();
    di_vec.clear();
    DeleteArray(&xt, xlen);
    delete[] invmap;
    return tu_vec.size();
}

/* extract rotation matrix based on TMscore8 */
void output_flexalign_rotation_matrix(const char* fname_matrix,
    const vector<vector<double> >&tu_vec, double t[3], double u[3][3])
{
    stringstream ss;
    char dest[1000];
    for (int hinge=0;hinge<tu_vec.size();hinge++)
    {
        tu2t_u(tu_vec[hinge],t,u);
        ss << "------ The rotation matrix to rotate Structure_1 to Structure_2 ------\n";
        sprintf(dest, "m %18s %14s %14s %14s\n", "t[m]", "u[m][0]", "u[m][1]", "u[m][2]");
        ss << string(dest);
        for (int k = 0; k < 3; k++)
        {
            sprintf(dest, "%d %18.10f %14.10f %14.10f %14.10f\n", k, t[k], u[k][0], u[k][1], u[k][2]);
            ss << string(dest);
        }
    }
    ss << "\nCode for rotating Structure 1 from (x,y,z) to (X,Y,Z):\n"
            "for(i=0; i<L; i++)\n"
            "{\n"
            "   X[i] = t[0] + u[0][0]*x[i] + u[0][1]*y[i] + u[0][2]*z[i];\n"
            "   Y[i] = t[1] + u[1][0]*x[i] + u[1][1]*y[i] + u[1][2]*z[i];\n"
            "   Z[i] = t[2] + u[2][0]*x[i] + u[2][1]*y[i] + u[2][2]*z[i];\n"
            "}\n";
    if (strcmp(fname_matrix,(char *)("-"))==0)
       cout<<ss.str();
    else
    {
        fstream fout;
        fout.open(fname_matrix, ios::out | ios::trunc);
        if (fout)
        {
            fout<<ss.str();
            fout.close();
        }
        else cout << "Open file to output rotation matrix fail.\n";
    }
    ss.str(string());
}

void output_flexalign_rasmol(const string xname, const string yname,
    const string fname_super,const vector<vector<double> >&tu_vec,
    double t[3], double u[3][3], const int ter_opt,
    const int mm_opt, const int split_opt, const int mirror_opt,
    const char *seqM, const char *seqxA, const char *seqyA,
    const vector<string>&resi_vec1, const vector<string>&resi_vec2,
    const string chainID1, const string chainID2,
    const int xlen, const int ylen, const double d0A, const int n_ali8,
    const double rmsd, const double TM1, const double Liden)
{
    stringstream buf;
    stringstream buf_all;
    stringstream buf_atm;
    stringstream buf_all_atm;
    stringstream buf_all_atm_lig;
    //stringstream buf_pdb;
    stringstream buf_tm;
    string line;
    double x[3];  // before transform
    double x1[3]; // after transform
    bool after_ter; // true if passed the "TER" line in PDB
    string asym_id; // chain ID
    
    map<string,int> resi2hinge_dict;
    int r,i,j;
    j=-1;
    char hinge_char=0;
    int ali_len=strlen(seqM);
    for (r=0;r<strlen(seqxA);r++)
    {
        if (seqxA[r]=='-') continue;
        j++;
        hinge_char=seqM[r];
        if (hinge_char==' ')
        {
            for (i=1;i<ali_len;i++)
            {
                if (r-i>=0 && seqM[r-i]!=' ')
                    hinge_char=seqM[r-i];
                else if (r+i<xlen && seqM[r+i]!=' ')
                    hinge_char=seqM[r+i];
                if (hinge_char!=' ') break;
            }
        }
        resi2hinge_dict[resi_vec1[j]]=hinge_char-'0';
    }
    string resi=resi_vec1[0];
    int read_resi=resi.size()-4;

    buf_tm<<"REMARK US-align"
        <<"\nREMARK Structure 1:"<<setw(11)<<left<<xname+chainID1<<" Size= "<<xlen
        <<"\nREMARK Structure 2:"<<setw(11)<<yname+chainID2<<right<<" Size= "<<ylen
        <<" (TM-score is normalized by "<<setw(4)<<ylen<<", d0="
        <<setiosflags(ios::fixed)<<setprecision(2)<<setw(6)<<d0A<<")"
        <<"\nREMARK Aligned length="<<setw(4)<<n_ali8<<", RMSD="
        <<setw(6)<<setiosflags(ios::fixed)<<setprecision(2)<<rmsd
        <<", TM-score="<<setw(7)<<setiosflags(ios::fixed)<<setprecision(5)<<TM1
        <<", ID="<<setw(5)<<setiosflags(ios::fixed)<<setprecision(3)
        <<((n_ali8>0)?Liden/n_ali8:0)<<endl;
    string rasmol_CA_header="load inline\nselect *A\nwireframe .45\nselect *B\nwireframe .20\nselect all\ncolor white\n";
    string rasmol_cartoon_header="load inline\nselect all\ncartoon\nselect *A\ncolor blue\nselect *B\ncolor red\nselect ligand\nwireframe 0.25\nselect solvent\nspacefill 0.25\nselect all\nexit\n"+buf_tm.str();
    if (!mm_opt) buf<<rasmol_CA_header;
    buf_all<<rasmol_CA_header;
    if (!mm_opt) buf_atm<<rasmol_cartoon_header;
    buf_all_atm<<rasmol_cartoon_header;
    buf_all_atm_lig<<rasmol_cartoon_header;

    /* selecting chains for -mol */
    string chain1_sele;
    string chain2_sele;
    if (!mm_opt)
    {
        if (split_opt==2 && ter_opt>=1) // align one chain from model 1
        {
            chain1_sele=chainID1.substr(1);
            chain2_sele=chainID2.substr(1);
        }
        else if (split_opt==2 && ter_opt==0) // align one chain from each model
        {
            for (i=1;i<chainID1.size();i++) if (chainID1[i]==',') break;
            chain1_sele=chainID1.substr(i+1);
            for (i=1;i<chainID2.size();i++) if (chainID2[i]==',') break;
            chain2_sele=chainID2.substr(i+1);
        }
    }


    /* for PDBx/mmCIF only */
    map<string,int> _atom_site;
    int atom_site_pos;
    vector<string> line_vec;
    string atom; // 4-character atom name
    string AA;   // 3-character residue name
    string inscode; // 1-character insertion code
    string model_index; // model index
    bool is_mmcif=false;

    /* used for CONECT record of chain1 */
    int ca_idx1=0; // all CA atoms
    int lig_idx1=0; // all atoms
    vector <int> idx_vec;

    /* used for CONECT record of chain2 */
    int ca_idx2=0; // all CA atoms
    int lig_idx2=0; // all atoms

    /* extract aligned region */
    vector<string> resi_aln1;
    vector<string> resi_aln2;
    int i1=-1;
    int i2=-1;
    if (!mm_opt)
    {
        for (i=0;i<strlen(seqM);i++)
        {
            i1+=(seqxA[i]!='-');
            i2+=(seqyA[i]!='-');
            if (seqM[i]==' ') continue;
            resi_aln1.push_back(resi_vec1[i1].substr(0,4));
            resi_aln2.push_back(resi_vec2[i2].substr(0,4));
            if (seqM[i]!=':') continue;
            buf    <<"select "<<resi_aln1.back()<<":A,"
                   <<resi_aln2.back()<<":B\ncolor red\n";
            buf_all<<"select "<<resi_aln1.back()<<":A,"
                   <<resi_aln2.back()<<":B\ncolor red\n";
        }
        buf<<"select all\nexit\n"<<buf_tm.str();
    }
    buf_all<<"select all\nexit\n"<<buf_tm.str();

    ifstream fin;
    /* read first file */
    after_ter=false;
    asym_id="";
    fin.open(xname.c_str());
    int hinge=0;
    while (fin.good())
    {
        getline(fin, line);
        if (ter_opt>=3 && line.compare(0,3,"TER")==0) after_ter=true;
        if (is_mmcif==false && line.size()>=54 &&
           (line.compare(0, 6, "ATOM  ")==0 ||
            line.compare(0, 6, "HETATM")==0)) // PDB format
        {
            if (line[16]!='A' && line[16]!=' ') continue;
            x[0]=atof(line.substr(30,8).c_str());
            x[1]=atof(line.substr(38,8).c_str());
            x[2]=atof(line.substr(46,8).c_str());
            if (mirror_opt) x[2]=-x[2];
            if (read_resi==1) resi=line.substr(22,5);
            else resi=line.substr(22,5)+line[21];
            hinge=0;
            if (resi2hinge_dict.count(resi)) hinge=resi2hinge_dict[resi];
            tu2t_u(tu_vec[hinge],t,u);
            transform(t, u, x, x1);
            //buf_pdb<<line.substr(0,30)<<setiosflags(ios::fixed)
                //<<setprecision(3)
                //<<setw(8)<<x1[0] <<setw(8)<<x1[1] <<setw(8)<<x1[2]
                //<<line.substr(54)<<'\n';

            if (after_ter && line.compare(0,6,"ATOM  ")==0) continue;
            lig_idx1++;
            buf_all_atm_lig<<line.substr(0,6)<<setw(5)<<lig_idx1
                <<line.substr(11,9)<<" A"<<line.substr(22,8)
                <<setiosflags(ios::fixed)<<setprecision(3)
                <<setw(8)<<x1[0]<<setw(8)<<x1[1] <<setw(8)<<x1[2]<<'\n';
            if (chain1_sele.size() && line[21]!=chain1_sele[0]) continue;
            if (after_ter || line.compare(0,6,"ATOM  ")) continue;
            if (ter_opt>=2)
            {
                if (ca_idx1 && asym_id.size() && asym_id!=line.substr(21,1)) 
                {
                    after_ter=true;
                    continue;
                }
                asym_id=line[21];
            }
            buf_all_atm<<"ATOM  "<<setw(5)<<lig_idx1
                <<line.substr(11,9)<<" A"<<line.substr(22,8)
                <<setiosflags(ios::fixed)<<setprecision(3)
                <<setw(8)<<x1[0]<<setw(8)<<x1[1] <<setw(8)<<x1[2]<<'\n';
            if (!mm_opt && find(resi_aln1.begin(),resi_aln1.end(),
                line.substr(22,4))!=resi_aln1.end())
            {
                buf_atm<<"ATOM  "<<setw(5)<<lig_idx1
                    <<line.substr(11,9)<<" A"<<line.substr(22,8)
                    <<setiosflags(ios::fixed)<<setprecision(3)
                    <<setw(8)<<x1[0]<<setw(8)<<x1[1] <<setw(8)<<x1[2]<<'\n';
            }
            if (line.substr(12,4)!=" CA " && line.substr(12,4)!=" C3'") continue;
            ca_idx1++;
            buf_all<<"ATOM  "<<setw(5)<<ca_idx1<<' '
                <<line.substr(12,4)<<' '<<line.substr(17,3)<<" A"<<line.substr(22,8)
                <<setiosflags(ios::fixed)<<setprecision(3)
                <<setw(8)<<x1[0]<<setw(8)<<x1[1]<<setw(8)<<x1[2]<<'\n';
            if (find(resi_aln1.begin(),resi_aln1.end(),
                line.substr(22,4))==resi_aln1.end()) continue;
            if (!mm_opt) buf<<"ATOM  "<<setw(5)<<ca_idx1<<' '
                <<line.substr(12,4)<<' '<<line.substr(17,3)<<" A"<<line.substr(22,8)
                <<setiosflags(ios::fixed)<<setprecision(3)
                <<setw(8)<<x1[0]<<setw(8)<<x1[1]<<setw(8)<<x1[2]<<'\n';
            idx_vec.push_back(ca_idx1);
        }
        else if (line.compare(0,5,"loop_")==0) // PDBx/mmCIF
        {
            while(1)
            {
                if (fin.good()) getline(fin, line);
                else PrintErrorAndQuit("ERROR! Unexpected end of "+xname);
                if (line.size()) break;
            }
            if (line.compare(0,11,"_atom_site.")) continue;
            _atom_site.clear();
            atom_site_pos=0;
            _atom_site[line.substr(11,line.size()-12)]=atom_site_pos;
            while(1)
            {
                if (fin.good()) getline(fin, line);
                else PrintErrorAndQuit("ERROR! Unexpected end of "+xname);
                if (line.size()==0) continue;
                if (line.compare(0,11,"_atom_site.")) break;
                _atom_site[line.substr(11,line.size()-12)]=++atom_site_pos;
            }

            if (is_mmcif==false)
            {
                //buf_pdb.str(string());
                is_mmcif=true;
            }

            while(1)
            {
                line_vec.clear();
                split(line,line_vec);
                if (line_vec[_atom_site["group_PDB"]]!="ATOM" &&
                    line_vec[_atom_site["group_PDB"]]!="HETATM") break;
                if (_atom_site.count("pdbx_PDB_model_num"))
                {
                    if (model_index.size() && model_index!=
                        line_vec[_atom_site["pdbx_PDB_model_num"]])
                        break;
                    model_index=line_vec[_atom_site["pdbx_PDB_model_num"]];
                }

                x[0]=atof(line_vec[_atom_site["Cartn_x"]].c_str());
                x[1]=atof(line_vec[_atom_site["Cartn_y"]].c_str());
                x[2]=atof(line_vec[_atom_site["Cartn_z"]].c_str());
                if (mirror_opt) x[2]=-x[2];


                if (_atom_site.count("auth_seq_id"))
                    resi=line_vec[_atom_site["auth_seq_id"]];
                else resi=line_vec[_atom_site["label_seq_id"]];
                if (_atom_site.count("pdbx_PDB_ins_code") && 
                    line_vec[_atom_site["pdbx_PDB_ins_code"]]!="?")
                    resi+=line_vec[_atom_site["pdbx_PDB_ins_code"]][0];
                else resi+=" ";
                if (read_resi>=2)
                {
                    if (_atom_site.count("auth_asym_id"))
                        asym_id=line_vec[_atom_site["auth_asym_id"]];
                    else asym_id=line_vec[_atom_site["label_asym_id"]];
                    if (asym_id==".") asym_id=" ";
                    resi+=asym_id[0];
                }
                hinge=0;
                if (resi2hinge_dict.count(resi)) hinge=resi2hinge_dict[resi];
                tu2t_u(tu_vec[hinge],t,u);
                transform(t, u, x, x1);

                if (_atom_site.count("label_alt_id")==0 || 
                    line_vec[_atom_site["label_alt_id"]]=="." ||
                    line_vec[_atom_site["label_alt_id"]]=="A")
                {
                    atom=line_vec[_atom_site["label_atom_id"]];
                    if (atom[0]=='"') atom=atom.substr(1);
                    if (atom.size() && atom[atom.size()-1]=='"')
                        atom=atom.substr(0,atom.size()-1);
                    if      (atom.size()==0) atom="    ";
                    else if (atom.size()==1) atom=" "+atom+"  ";
                    else if (atom.size()==2) atom=" "+atom+" ";
                    else if (atom.size()==3) atom=" "+atom;
                    else if (atom.size()>=5) atom=atom.substr(0,4);
            
                    AA=line_vec[_atom_site["label_comp_id"]]; // residue name
                    if      (AA.size()==1) AA="  "+AA;
                    else if (AA.size()==2) AA=" " +AA;
                    else if (AA.size()>=4) AA=AA.substr(0,3);
                
                    if (_atom_site.count("auth_seq_id"))
                        resi=line_vec[_atom_site["auth_seq_id"]];
                    else resi=line_vec[_atom_site["label_seq_id"]];
                    while (resi.size()<4) resi=' '+resi;
                    if (resi.size()>4) resi=resi.substr(0,4);
                
                    inscode=' ';
                    if (_atom_site.count("pdbx_PDB_ins_code") && 
                        line_vec[_atom_site["pdbx_PDB_ins_code"]]!="?")
                        inscode=line_vec[_atom_site["pdbx_PDB_ins_code"]][0];

                    if (_atom_site.count("auth_asym_id"))
                    {
                        if (chain1_sele.size()) after_ter
                            =line_vec[_atom_site["auth_asym_id"]]!=chain1_sele;
                        else if (ter_opt>=2 && ca_idx1 && asym_id.size() && 
                            asym_id!=line_vec[_atom_site["auth_asym_id"]])
                            after_ter=true;
                        asym_id=line_vec[_atom_site["auth_asym_id"]];
                    }
                    else if (_atom_site.count("label_asym_id"))
                    {
                        if (chain1_sele.size()) after_ter
                            =line_vec[_atom_site["label_asym_id"]]!=chain1_sele;
                        if (ter_opt>=2 && ca_idx1 && asym_id.size() && 
                            asym_id!=line_vec[_atom_site["label_asym_id"]])
                            after_ter=true;
                        asym_id=line_vec[_atom_site["label_asym_id"]];
                    }
                    //buf_pdb<<left<<setw(6)
                        //<<line_vec[_atom_site["group_PDB"]]<<right
                        //<<setw(5)<<lig_idx1%100000<<' '<<atom<<' '
                        //<<AA<<" "<<asym_id[asym_id.size()-1]
                        //<<resi<<inscode<<"   "
                        //<<setiosflags(ios::fixed)<<setprecision(3)
                        //<<setw(8)<<x1[0]
                        //<<setw(8)<<x1[1]
                        //<<setw(8)<<x1[2]<<'\n';

                    if (after_ter==false ||
                        line_vec[_atom_site["group_pdb"]]=="HETATM")
                    {
                        lig_idx1++;
                        buf_all_atm_lig<<left<<setw(6)
                            <<line_vec[_atom_site["group_PDB"]]<<right
                            <<setw(5)<<lig_idx1%100000<<' '<<atom<<' '
                            <<AA<<" A"<<resi<<inscode<<"   "
                            <<setiosflags(ios::fixed)<<setprecision(3)
                            <<setw(8)<<x1[0]
                            <<setw(8)<<x1[1]
                            <<setw(8)<<x1[2]<<'\n';
                        if (after_ter==false &&
                            line_vec[_atom_site["group_PDB"]]=="ATOM")
                        {
                            buf_all_atm<<"ATOM  "<<setw(6)
                                <<setw(5)<<lig_idx1%100000<<' '<<atom<<' '
                                <<AA<<" A"<<resi<<inscode<<"   "
                                <<setiosflags(ios::fixed)<<setprecision(3)
                                <<setw(8)<<x1[0]
                                <<setw(8)<<x1[1]
                                <<setw(8)<<x1[2]<<'\n';
                            if (!mm_opt && find(resi_aln1.begin(),
                                resi_aln1.end(),resi)!=resi_aln1.end())
                            {
                                buf_atm<<"ATOM  "<<setw(6)
                                    <<setw(5)<<lig_idx1%100000<<' '
                                    <<atom<<' '<<AA<<" A"<<resi<<inscode<<"   "
                                    <<setiosflags(ios::fixed)<<setprecision(3)
                                    <<setw(8)<<x1[0]
                                    <<setw(8)<<x1[1]
                                    <<setw(8)<<x1[2]<<'\n';
                            }
                            if (atom==" CA " || atom==" C3'")
                            {
                                ca_idx1++;
            //mm_opt, split_opt, mirror_opt, chainID1,chainID2);
                                buf_all<<"ATOM  "<<setw(6)
                                    <<setw(5)<<ca_idx1%100000<<' '<<atom<<' '
                                    <<AA<<" A"<<resi<<inscode<<"   "
                                    <<setiosflags(ios::fixed)<<setprecision(3)
                                    <<setw(8)<<x1[0]
                                    <<setw(8)<<x1[1]
                                    <<setw(8)<<x1[2]<<'\n';
                                if (!mm_opt && find(resi_aln1.begin(),
                                    resi_aln1.end(),resi)!=resi_aln1.end())
                                {
                                    buf<<"ATOM  "<<setw(6)
                                    <<setw(5)<<ca_idx1%100000<<' '<<atom<<' '
                                    <<AA<<" A"<<resi<<inscode<<"   "
                                    <<setiosflags(ios::fixed)<<setprecision(3)
                                    <<setw(8)<<x1[0]
                                    <<setw(8)<<x1[1]
                                    <<setw(8)<<x1[2]<<'\n';
                                    idx_vec.push_back(ca_idx1);
                                }
                            }
                        }
                    }
                }

                while(1)
                {
                    if (fin.good()) getline(fin, line);
                    else break;
                    if (line.size()) break;
                }
            }
        }
        else if (line.size() && is_mmcif==false)
        {
            //buf_pdb<<line<<'\n';
            if (ter_opt>=1 && line.compare(0,3,"END")==0) break;
        }
    }
    fin.close();
    if (!mm_opt) buf<<"TER\n";
    buf_all<<"TER\n";
    if (!mm_opt) buf_atm<<"TER\n";
    buf_all_atm<<"TER\n";
    buf_all_atm_lig<<"TER\n";
    for (i=1;i<ca_idx1;i++) buf_all<<"CONECT"
        <<setw(5)<<i%100000<<setw(5)<<(i+1)%100000<<'\n';
    if (!mm_opt) for (i=1;i<idx_vec.size();i++) buf<<"CONECT"
        <<setw(5)<<idx_vec[i-1]%100000<<setw(5)<<idx_vec[i]%100000<<'\n';
    idx_vec.clear();

    /* read second file */
    after_ter=false;
    asym_id="";
    fin.open(yname.c_str());
    while (fin.good())
    {
        getline(fin, line);
        if (ter_opt>=3 && line.compare(0,3,"TER")==0) after_ter=true;
        if (line.size()>=54 && (line.compare(0, 6, "ATOM  ")==0 ||
            line.compare(0, 6, "HETATM")==0)) // PDB format
        {
            if (line[16]!='A' && line[16]!=' ') continue;
            if (after_ter && line.compare(0,6,"ATOM  ")==0) continue;
            lig_idx2++;
            buf_all_atm_lig<<line.substr(0,6)<<setw(5)<<lig_idx1+lig_idx2
                <<line.substr(11,9)<<" B"<<line.substr(22,32)<<'\n';
            if (chain1_sele.size() && line[21]!=chain1_sele[0]) continue;
            if (after_ter || line.compare(0,6,"ATOM  ")) continue;
            if (ter_opt>=2)
            {
                if (ca_idx2 && asym_id.size() && asym_id!=line.substr(21,1))
                {
                    after_ter=true;
                    continue;
                }
                asym_id=line[21];
            }
            buf_all_atm<<"ATOM  "<<setw(5)<<lig_idx1+lig_idx2
                <<line.substr(11,9)<<" B"<<line.substr(22,32)<<'\n';
            if (!mm_opt && find(resi_aln2.begin(),resi_aln2.end(),
                line.substr(22,4))!=resi_aln2.end())
            {
                buf_atm<<"ATOM  "<<setw(5)<<lig_idx1+lig_idx2
                    <<line.substr(11,9)<<" B"<<line.substr(22,32)<<'\n';
            }
            if (line.substr(12,4)!=" CA " && line.substr(12,4)!=" C3'") continue;
            ca_idx2++;
            buf_all<<"ATOM  "<<setw(5)<<ca_idx1+ca_idx2<<' '<<line.substr(12,4)
                <<' '<<line.substr(17,3)<<" B"<<line.substr(22,32)<<'\n';
            if (find(resi_aln2.begin(),resi_aln2.end(),line.substr(22,4)
                )==resi_aln2.end()) continue;
            if (!mm_opt) buf<<"ATOM  "<<setw(5)<<ca_idx1+ca_idx2<<' '
                <<line.substr(12,4)<<' '<<line.substr(17,3)<<" B"
                <<line.substr(22,32)<<'\n';
            idx_vec.push_back(ca_idx1+ca_idx2);
        }
        else if (line.compare(0,5,"loop_")==0) // PDBx/mmCIF
        {
            while(1)
            {
                if (fin.good()) getline(fin, line);
                else PrintErrorAndQuit("ERROR! Unexpected end of "+yname);
                if (line.size()) break;
            }
            if (line.compare(0,11,"_atom_site.")) continue;
            _atom_site.clear();
            atom_site_pos=0;
            _atom_site[line.substr(11,line.size()-12)]=atom_site_pos;
            while(1)
            {
                if (fin.good()) getline(fin, line);
                else PrintErrorAndQuit("ERROR! Unexpected end of "+yname);
                if (line.size()==0) continue;
                if (line.compare(0,11,"_atom_site.")) break;
                _atom_site[line.substr(11,line.size()-12)]=++atom_site_pos;
            }

            while(1)
            {
                line_vec.clear();
                split(line,line_vec);
                if (line_vec[_atom_site["group_PDB"]]!="ATOM" &&
                    line_vec[_atom_site["group_PDB"]]!="HETATM") break;
                if (_atom_site.count("pdbx_PDB_model_num"))
                {
                    if (model_index.size() && model_index!=
                        line_vec[_atom_site["pdbx_PDB_model_num"]])
                        break;
                    model_index=line_vec[_atom_site["pdbx_PDB_model_num"]];
                }

                if (_atom_site.count("label_alt_id")==0 || 
                    line_vec[_atom_site["label_alt_id"]]=="." ||
                    line_vec[_atom_site["label_alt_id"]]=="A")
                {
                    atom=line_vec[_atom_site["label_atom_id"]];
                    if (atom[0]=='"') atom=atom.substr(1);
                    if (atom.size() && atom[atom.size()-1]=='"')
                        atom=atom.substr(0,atom.size()-1);
                    if      (atom.size()==0) atom="    ";
                    else if (atom.size()==1) atom=" "+atom+"  ";
                    else if (atom.size()==2) atom=" "+atom+" ";
                    else if (atom.size()==3) atom=" "+atom;
                    else if (atom.size()>=5) atom=atom.substr(0,4);
            
                    AA=line_vec[_atom_site["label_comp_id"]]; // residue name
                    if      (AA.size()==1) AA="  "+AA;
                    else if (AA.size()==2) AA=" " +AA;
                    else if (AA.size()>=4) AA=AA.substr(0,3);
                
                    if (_atom_site.count("auth_seq_id"))
                        resi=line_vec[_atom_site["auth_seq_id"]];
                    else resi=line_vec[_atom_site["label_seq_id"]];
                    while (resi.size()<4) resi=' '+resi;
                    if (resi.size()>4) resi=resi.substr(0,4);
                
                    inscode=' ';
                    if (_atom_site.count("pdbx_PDB_ins_code") && 
                        line_vec[_atom_site["pdbx_PDB_ins_code"]]!="?")
                        inscode=line_vec[_atom_site["pdbx_PDB_ins_code"]][0];
                    
                    if (_atom_site.count("auth_asym_id"))
                    {
                        if (chain2_sele.size()) after_ter
                            =line_vec[_atom_site["auth_asym_id"]]!=chain2_sele;
                        if (ter_opt>=2 && ca_idx2 && asym_id.size() && 
                            asym_id!=line_vec[_atom_site["auth_asym_id"]])
                            after_ter=true;
                        asym_id=line_vec[_atom_site["auth_asym_id"]];
                    }
                    else if (_atom_site.count("label_asym_id"))
                    {
                        if (chain2_sele.size()) after_ter
                            =line_vec[_atom_site["label_asym_id"]]!=chain2_sele;
                        if (ter_opt>=2 && ca_idx2 && asym_id.size() && 
                            asym_id!=line_vec[_atom_site["label_asym_id"]])
                            after_ter=true;
                        asym_id=line_vec[_atom_site["label_asym_id"]];
                    }
                    if (after_ter==false || 
                        line_vec[_atom_site["group_PDB"]]=="HETATM")
                    {
                        lig_idx2++;
                        buf_all_atm_lig<<left<<setw(6)
                            <<line_vec[_atom_site["group_PDB"]]<<right
                            <<setw(5)<<(lig_idx1+lig_idx2)%100000<<' '
                            <<atom<<' '<<AA<<" B"<<resi<<inscode<<"   "
                            <<setw(8)<<line_vec[_atom_site["Cartn_x"]]
                            <<setw(8)<<line_vec[_atom_site["Cartn_y"]]
                            <<setw(8)<<line_vec[_atom_site["Cartn_z"]]
                            <<'\n';
                        if (after_ter==false &&
                            line_vec[_atom_site["group_PDB"]]=="ATOM")
                        {
                            buf_all_atm<<"ATOM  "<<setw(6)
                                <<setw(5)<<(lig_idx1+lig_idx2)%100000<<' '
                                <<atom<<' '<<AA<<" B"<<resi<<inscode<<"   "
                                <<setw(8)<<line_vec[_atom_site["Cartn_x"]]
                                <<setw(8)<<line_vec[_atom_site["Cartn_y"]]
                                <<setw(8)<<line_vec[_atom_site["Cartn_z"]]
                                <<'\n';
                            if (!mm_opt && find(resi_aln2.begin(),
                                resi_aln2.end(),resi)!=resi_aln2.end())
                            {
                                buf_atm<<"ATOM  "<<setw(6)
                                    <<setw(5)<<(lig_idx1+lig_idx2)%100000<<' '
                                    <<atom<<' '<<AA<<" B"<<resi<<inscode<<"   "
                                    <<setw(8)<<line_vec[_atom_site["Cartn_x"]]
                                    <<setw(8)<<line_vec[_atom_site["Cartn_y"]]
                                    <<setw(8)<<line_vec[_atom_site["Cartn_z"]]
                                    <<'\n';
                            }
                            if (atom==" CA " || atom==" C3'")
                            {
                                ca_idx2++;
                                buf_all<<"ATOM  "<<setw(6)
                                    <<setw(5)<<(ca_idx1+ca_idx2)%100000
                                    <<' '<<atom<<' '<<AA<<" B"<<resi<<inscode<<"   "
                                    <<setw(8)<<line_vec[_atom_site["Cartn_x"]]
                                    <<setw(8)<<line_vec[_atom_site["Cartn_y"]]
                                    <<setw(8)<<line_vec[_atom_site["Cartn_z"]]
                                    <<'\n';
                                if (!mm_opt && find(resi_aln2.begin(),
                                    resi_aln2.end(),resi)!=resi_aln2.end())
                                {
                                    buf<<"ATOM  "<<setw(6)
                                    <<setw(5)<<(ca_idx1+ca_idx2)%100000
                                    <<' '<<atom<<' '<<AA<<" B"<<resi<<inscode<<"   "
                                    <<setw(8)<<line_vec[_atom_site["Cartn_x"]]
                                    <<setw(8)<<line_vec[_atom_site["Cartn_y"]]
                                    <<setw(8)<<line_vec[_atom_site["Cartn_z"]]
                                    <<'\n';
                                    idx_vec.push_back(ca_idx1+ca_idx2);
                                }
                            }
                        }
                    }
                }

                if (fin.good()) getline(fin, line);
                else break;
            }
        }
        else if (line.size())
        {
            if (ter_opt>=1 && line.compare(0,3,"END")==0) break;
        }
    }
    fin.close();
    if (!mm_opt) buf<<"TER\n";
    buf_all<<"TER\n";
    if (!mm_opt) buf_atm<<"TER\n";
    buf_all_atm<<"TER\n";
    buf_all_atm_lig<<"TER\n";
    for (i=ca_idx1+1;i<ca_idx1+ca_idx2;i++) buf_all<<"CONECT"
        <<setw(5)<<i%100000<<setw(5)<<(i+1)%100000<<'\n';
    for (i=1;i<idx_vec.size();i++) buf<<"CONECT"
        <<setw(5)<<idx_vec[i-1]%100000<<setw(5)<<idx_vec[i]%100000<<'\n';
    idx_vec.clear();

    /* write pymol script */
    ofstream fp;
    /*
    stringstream buf_pymol;
    vector<string> pml_list;
    pml_list.push_back(fname_super+"");
    pml_list.push_back(fname_super+"_atm");
    pml_list.push_back(fname_super+"_all");
    pml_list.push_back(fname_super+"_all_atm");
    pml_list.push_back(fname_super+"_all_atm_lig");
    for (i=0;i<pml_list.size();i++)
    {
        buf_pymol<<"#!/usr/bin/env pymol\n"
            <<"load "<<pml_list[i]<<"\n"
            <<"hide all\n"
            <<((i==0 || i==2)?("show stick\n"):("show cartoon\n"))
            <<"color blue, chain A\n"
            <<"color red, chain B\n"
            <<"set ray_shadow, 0\n"
            <<"set stick_radius, 0.3\n"
            <<"set sphere_scale, 0.25\n"
            <<"show stick, not polymer\n"
            <<"show sphere, not polymer\n"
            <<"bg_color white\n"
            <<"set transparency=0.2\n"
            <<"zoom polymer\n"
            <<endl;
        fp.open((pml_list[i]+".pml").c_str());
        fp<<buf_pymol.str();
        fp.close();
        buf_pymol.str(string());
        pml_list[i].clear();
    }
    pml_list.clear();
    */
    
    /* write rasmol script */
    if (!mm_opt)
    {
        fp.open((fname_super).c_str());
        fp<<buf.str();
        fp.close();
    }
    fp.open((fname_super+"_all").c_str());
    fp<<buf_all.str();
    fp.close();
    if (!mm_opt)
    {
        fp.open((fname_super+"_atm").c_str());
        fp<<buf_atm.str();
        fp.close();
    }
    fp.open((fname_super+"_all_atm").c_str());
    fp<<buf_all_atm.str();
    fp.close();
    fp.open((fname_super+"_all_atm_lig").c_str());
    fp<<buf_all_atm_lig.str();
    fp.close();
    //fp.open((fname_super+".pdb").c_str());
    //fp<<buf_pdb.str();
    //fp.close();

    /* clear stream */
    buf.str(string());
    buf_all.str(string());
    buf_atm.str(string());
    buf_all_atm.str(string());
    buf_all_atm_lig.str(string());
    //buf_pdb.str(string());
    buf_tm.str(string());
    resi_aln1.clear();
    resi_aln2.clear();
    asym_id.clear();
    line_vec.clear();
    atom.clear();
    AA.clear();
    resi.clear();
    inscode.clear();
    model_index.clear();
}

void output_flexalign_pymol(const string xname, const string yname,
    const string fname_super, const vector<vector<double> >&tu_vec,
    double t[3], double u[3][3], const int ter_opt, 
    const int mm_opt, const int split_opt, const int mirror_opt,
    const char *seqM, const char *seqxA, const char *seqyA,
    const vector<string>&resi_vec1, const vector<string>&resi_vec2,
    const string chainID1, const string chainID2)
{
    int compress_type=0; // uncompressed file
    ifstream fin;
#ifndef REDI_PSTREAM_H_SEEN
    ifstream fin_gz;
#else
    redi::ipstream fin_gz; // if file is compressed
    if (xname.size()>=3 && 
        xname.substr(xname.size()-3,3)==".gz")
    {
        fin_gz.open("gunzip -c "+xname);
        compress_type=1;
    }
    else if (xname.size()>=4 && 
        xname.substr(xname.size()-4,4)==".bz2")
    {
        fin_gz.open("bzcat "+xname);
        compress_type=2;
    }
    else
#endif
        fin.open(xname.c_str());

    map<string,int> resi2hinge_dict;
    int r,i,j;
    j=-1;
    char hinge_char=0;
    int xlen=resi_vec1.size();
    int ali_len=strlen(seqM);
    for (r=0;r<strlen(seqxA);r++)
    {
        if (seqxA[r]=='-') continue;
        j++;
        hinge_char=seqM[r];
        if (hinge_char==' ')
        {
            for (i=1;i<ali_len;i++)
            {
                if (r-i>=0 && seqM[r-i]!=' ')
                    hinge_char=seqM[r-i];
                else if (r+i<xlen && seqM[r+i]!=' ')
                    hinge_char=seqM[r+i];
                if (hinge_char!=' ') break;
            }
        }
        resi2hinge_dict[resi_vec1[j]]=hinge_char-'0';
    }
    string resi=resi_vec1[0];
    int read_resi=resi.size()-4;

    stringstream buf;
    stringstream buf_pymol;
    string line;
    double x[3];  // before transform
    double x1[3]; // after transform

    /* for PDBx/mmCIF only */
    map<string,int> _atom_site;
    size_t atom_site_pos;
    vector<string> line_vec;
    int infmt=-1; // 0 - PDB, 3 - PDBx/mmCIF
    int hinge=0;
    string asym_id="."; // this is similar to chainID, except that
                        // chainID is char while asym_id is a string
                        // with possibly multiple char
    while (compress_type?fin_gz.good():fin.good())
    {
        if (compress_type) getline(fin_gz, line);
        else               getline(fin, line);
        if (line.compare(0, 6, "ATOM  ")==0 || 
            line.compare(0, 6, "HETATM")==0) // PDB format
        {
            infmt=0;
            x[0]=atof(line.substr(30,8).c_str());
            x[1]=atof(line.substr(38,8).c_str());
            x[2]=atof(line.substr(46,8).c_str());
            if (mirror_opt) x[2]=-x[2];
            if (read_resi==1) resi=line.substr(22,5);
            else resi=line.substr(22,5)+line[21];
            hinge=0;
            if (resi2hinge_dict.count(resi)) hinge=resi2hinge_dict[resi];
            tu2t_u(tu_vec[hinge],t,u);
            transform(t, u, x, x1);
            buf<<line.substr(0,30)<<setiosflags(ios::fixed)
                <<setprecision(3)
                <<setw(8)<<x1[0] <<setw(8)<<x1[1] <<setw(8)<<x1[2]
                <<line.substr(54)<<'\n';
        }
        else if (line.compare(0,5,"loop_")==0) // PDBx/mmCIF
        {
            infmt=3;
            buf<<line<<'\n';
            while(1)
            {
                if (compress_type) 
                {
                    if (fin_gz.good()) getline(fin_gz, line);
                    else PrintErrorAndQuit("ERROR! Unexpected end of "+xname);
                }
                else
                {
                    if (fin.good()) getline(fin, line);
                    else PrintErrorAndQuit("ERROR! Unexpected end of "+xname);
                }
                if (line.size()) break;
            }
            buf<<line<<'\n';
            if (line.compare(0,11,"_atom_site.")) continue;
            _atom_site.clear();
            atom_site_pos=0;
            _atom_site[Trim(line.substr(11))]=atom_site_pos;
            while(1)
            {
                while(1)
                {
                    if (compress_type) 
                    {
                        if (fin_gz.good()) getline(fin_gz, line);
                        else PrintErrorAndQuit("ERROR! Unexpected end of "+xname);
                    }
                    else
                    {
                        if (fin.good()) getline(fin, line);
                        else PrintErrorAndQuit("ERROR! Unexpected end of "+xname);
                    }
                    if (line.size()) break;
                }
                if (line.compare(0,11,"_atom_site.")) break;
                _atom_site[Trim(line.substr(11))]=++atom_site_pos;
                buf<<line<<'\n';
            }

            if (_atom_site.count("group_PDB")*
                _atom_site.count("Cartn_x")*
                _atom_site.count("Cartn_y")*
                _atom_site.count("Cartn_z")==0)
            {
                buf<<line<<'\n';
                cerr<<"Warning! Missing one of the following _atom_site data items: group_PDB, Cartn_x, Cartn_y, Cartn_z"<<endl;
                continue;
            }

            while(1)
            {
                line_vec.clear();
                split(line,line_vec);
                if (line_vec[_atom_site["group_PDB"]]!="ATOM" &&
                    line_vec[_atom_site["group_PDB"]]!="HETATM") break;

                x[0]=atof(line_vec[_atom_site["Cartn_x"]].c_str());
                x[1]=atof(line_vec[_atom_site["Cartn_y"]].c_str());
                x[2]=atof(line_vec[_atom_site["Cartn_z"]].c_str());
                if (mirror_opt) x[2]=-x[2];



                if (_atom_site.count("auth_seq_id"))
                    resi=line_vec[_atom_site["auth_seq_id"]];
                else resi=line_vec[_atom_site["label_seq_id"]];
                if (_atom_site.count("pdbx_PDB_ins_code") && 
                    line_vec[_atom_site["pdbx_PDB_ins_code"]]!="?")
                    resi+=line_vec[_atom_site["pdbx_PDB_ins_code"]][0];
                else resi+=" ";
                if (read_resi>=2)
                {
                    if (_atom_site.count("auth_asym_id"))
                        asym_id=line_vec[_atom_site["auth_asym_id"]];
                    else asym_id=line_vec[_atom_site["label_asym_id"]];
                    if (asym_id==".") asym_id=" ";
                    resi+=asym_id[0];
                }
                hinge=0;
                if (resi2hinge_dict.count(resi)) hinge=resi2hinge_dict[resi];
                tu2t_u(tu_vec[hinge],t,u);
                transform(t, u, x, x1);

                for (atom_site_pos=0; atom_site_pos<_atom_site.size(); atom_site_pos++)
                {
                    if (atom_site_pos==_atom_site["Cartn_x"])
                        buf<<setiosflags(ios::fixed)<<setprecision(3)
                           <<setw(8)<<x1[0]<<' ';
                    else if (atom_site_pos==_atom_site["Cartn_y"])
                        buf<<setiosflags(ios::fixed)<<setprecision(3)
                           <<setw(8)<<x1[1]<<' ';
                    else if (atom_site_pos==_atom_site["Cartn_z"])
                        buf<<setiosflags(ios::fixed)<<setprecision(3)
                           <<setw(8)<<x1[2]<<' ';
                    else buf<<line_vec[atom_site_pos]<<' ';
                }
                buf<<'\n';

                if (compress_type && fin_gz.good()) getline(fin_gz, line);
                else if (!compress_type && fin.good()) getline(fin, line);
                else break;
            }
            if (compress_type?fin_gz.good():fin.good()) buf<<line<<'\n';
        }
        else if (line.size())
        {
            buf<<line<<'\n';
            if (ter_opt>=1 && line.compare(0,3,"END")==0) break;
        }
    }
    if (compress_type) fin_gz.close();
    else               fin.close();

    string fname_super_full=fname_super;
    if (infmt==0)      fname_super_full+=".pdb";
    else if (infmt==3) fname_super_full+=".cif";
    ofstream fp;
    fp.open(fname_super_full.c_str());
    fp<<buf.str();
    fp.close();
    buf.str(string()); // clear stream

    string chain1_sele;
    string chain2_sele;
    if (!mm_opt)
    {
        if (split_opt==2 && ter_opt>=1) // align one chain from model 1
        {
            chain1_sele=" and c. "+chainID1.substr(1);
            chain2_sele=" and c. "+chainID2.substr(1);
        }
        else if (split_opt==2 && ter_opt==0) // align one chain from each model
        {
            for (i=1;i<chainID1.size();i++) if (chainID1[i]==',') break;
            chain1_sele=" and c. "+chainID1.substr(i+1);
            for (i=1;i<chainID2.size();i++) if (chainID2[i]==',') break;
            chain2_sele=" and c. "+chainID2.substr(i+1);
        }
    }

    /* extract aligned region */
    int i1=-1;
    int i2=-1;
    string resi1_sele;
    string resi2_sele;
    string resi1_bond;
    string resi2_bond;
    string prev_resi1;
    string prev_resi2;
    string curr_resi1;
    string curr_resi2;
    if (mm_opt)
    {
        ;
    }
    else
    {
        for (i=0;i<strlen(seqM);i++)
        {
            i1+=(seqxA[i]!='-' && seqxA[i]!='*');
            i2+=(seqyA[i]!='-');
            if (seqM[i]==' ' || seqxA[i]=='*') continue;
            curr_resi1=resi_vec1[i1].substr(0,4);
            curr_resi2=resi_vec2[i2].substr(0,4);
            if (resi1_sele.size()==0)
                resi1_sele =    "i. "+curr_resi1;
            else
            {
                resi1_sele+=" or i. "+curr_resi1;
                resi1_bond+="bond structure1 and i. "+prev_resi1+
                                              ", i. "+curr_resi1+"\n";
            }
            if (resi2_sele.size()==0)
                resi2_sele =    "i. "+curr_resi2;
            else
            {
                resi2_sele+=" or i. "+curr_resi2;
                resi2_bond+="bond structure2 and i. "+prev_resi2+
                                              ", i. "+curr_resi2+"\n";
            }
            prev_resi1=curr_resi1;
            prev_resi2=curr_resi2;
            //if (seqM[i]!=':') continue;
        }
        if (resi1_sele.size()) resi1_sele=" and ( "+resi1_sele+")";
        if (resi2_sele.size()) resi2_sele=" and ( "+resi2_sele+")";
    }

    /* write pymol script */
    vector<string> pml_list;
    pml_list.push_back(fname_super+"");
    pml_list.push_back(fname_super+"_atm");
    pml_list.push_back(fname_super+"_all");
    pml_list.push_back(fname_super+"_all_atm");
    pml_list.push_back(fname_super+"_all_atm_lig");

    for (int p=0;p<pml_list.size();p++)
    {
        if (mm_opt && p<=1) continue;
        buf_pymol
            <<"#!/usr/bin/env pymol\n"
            <<"cmd.load(\""<<fname_super_full<<"\", \"structure1\")\n"
            <<"cmd.load(\""<<yname<<"\", \"structure2\")\n"
            <<"hide all\n"
            <<"set all_states, "<<((ter_opt==0)?"on":"off")<<'\n';
        if (p==0) // .pml
        {
            if (chain1_sele.size()) buf_pymol
                <<"remove structure1 and not "<<chain1_sele.substr(4)<<"\n";
            if (chain2_sele.size()) buf_pymol
                <<"remove structure2 and not "<<chain2_sele.substr(4)<<"\n";
            buf_pymol
                <<"remove not n. CA and not n. C3'\n"
                <<resi1_bond
                <<resi2_bond
                <<"show stick, structure1"<<chain1_sele<<resi1_sele<<"\n"
                <<"show stick, structure2"<<chain2_sele<<resi2_sele<<"\n";
        }
        else if (p==1) // _atm.pml
        {
            buf_pymol
                <<"show cartoon, structure1"<<chain1_sele<<resi1_sele<<"\n"
                <<"show cartoon, structure2"<<chain2_sele<<resi2_sele<<"\n";
        }
        else if (p==2) // _all.pml
        {
            buf_pymol
                <<"show ribbon, structure1"<<chain1_sele<<"\n"
                <<"show ribbon, structure2"<<chain2_sele<<"\n";
        }
        else if (p==3) // _all_atm.pml
        {
            buf_pymol
                <<"show cartoon, structure1"<<chain1_sele<<"\n"
                <<"show cartoon, structure2"<<chain2_sele<<"\n";
        }
        else if (p==4) // _all_atm_lig.pml
        {
            buf_pymol
                <<"show cartoon, structure1\n"
                <<"show cartoon, structure2\n"
                <<"show stick, not polymer\n"
                <<"show sphere, not polymer\n";
        }
        buf_pymol
            <<"color blue, structure1\n"
            <<"color red, structure2\n"
            <<"set ribbon_width, 6\n"
            <<"set stick_radius, 0.3\n"
            <<"set sphere_scale, 0.25\n"
            <<"set ray_shadow, 0\n"
            <<"bg_color white\n"
            <<"set transparency=0.2\n"
            <<"zoom polymer and ((structure1"<<chain1_sele
            <<") or (structure2"<<chain2_sele<<"))\n"
            <<endl;

        fp.open((pml_list[p]+".pml").c_str());
        fp<<buf_pymol.str();
        fp.close();
        buf_pymol.str(string());
    }

    /* clean up */
    pml_list.clear();
    
    resi1_sele.clear();
    resi2_sele.clear();
    
    resi1_bond.clear();
    resi2_bond.clear();
    
    prev_resi1.clear();
    prev_resi2.clear();

    curr_resi1.clear();
    curr_resi2.clear();

    chain1_sele.clear();
    chain2_sele.clear();
    resi2hinge_dict.clear();
}

//output the final results
void output_flexalign_results(const string xname, const string yname,
    const string chainID1, const string chainID2,
    const int xlen, const int ylen, double t[3], double u[3][3],
    const vector<vector<double> >&tu_vec, const double TM1, const double TM2,
    const double TM3, const double TM4, const double TM5,
    const double rmsd, const double d0_out, const char *seqM,
    const char *seqxA, const char *seqyA, const double Liden,
    const int n_ali8, const int L_ali, const double TM_ali,
    const double rmsd_ali, const double TM_0, const double d0_0,
    const double d0A, const double d0B, const double Lnorm_ass,
    const double d0_scale, const double d0a, const double d0u,
    const char* fname_matrix, const int outfmt_opt, const int ter_opt,
    const int mm_opt, const int split_opt, const int o_opt,
    const string fname_super, const int i_opt, const int a_opt,
    const bool u_opt, const bool d_opt, const int mirror_opt,
    const vector<string>&resi_vec1, const vector<string>&resi_vec2)
{
    if (outfmt_opt<=0)
    {
        printf("\nName of Structure_1: %s%s (to be superimposed onto Structure_2)\n",
            xname.c_str(), chainID1.c_str());
        printf("Name of Structure_2: %s%s\n", yname.c_str(), chainID2.c_str());
        printf("Length of Structure_1: %d residues\n", xlen);
        printf("Length of Structure_2: %d residues\n\n", ylen);

        if (i_opt)
            printf("User-specified initial alignment: TM/Lali/rmsd = %7.5lf, %4d, %6.3lf\n", TM_ali, L_ali, rmsd_ali);

        printf("Aligned length= %d, RMSD= %6.2f, Seq_ID=n_identical/n_aligned= %4.3f\n", n_ali8, rmsd, (n_ali8>0)?Liden/n_ali8:0);
        printf("TM-score= %6.5f (normalized by length of Structure_1: L=%d, d0=%.2f)\n", TM2, xlen, d0B);
        printf("TM-score= %6.5f (normalized by length of Structure_2: L=%d, d0=%.2f)\n", TM1, ylen, d0A);

        if (a_opt==1)
            printf("TM-score= %6.5f (if normalized by average length of two structures: L=%.1f, d0=%.2f)\n", TM3, (xlen+ylen)*0.5, d0a);
        if (u_opt)
            printf("TM-score= %6.5f (normalized by user-specified L=%.2f and d0=%.2f)\n", TM4, Lnorm_ass, d0u);
        if (d_opt)
            printf("TM-score= %6.5f (scaled by user-specified d0=%.2f, and L=%d)\n", TM5, d0_scale, ylen);
        printf("(You should use TM-score normalized by length of the reference structure)\n");
    
        //output alignment
        printf("\n([0-9] denote different aligned fragment pairs separated by different hinges)\n");
        printf("%s\n", seqxA);
        printf("%s\n", seqM);
        printf("%s\n", seqyA);
    }
    else if (outfmt_opt==1)
    {
        printf(">%s%s\tL=%d\td0=%.2f\tseqID=%.3f\tTM-score=%.5f\n",
            xname.c_str(), chainID1.c_str(), xlen, d0B, Liden/xlen, TM2);
        printf("%s\n", seqxA);
        printf(">%s%s\tL=%d\td0=%.2f\tseqID=%.3f\tTM-score=%.5f\n",
            yname.c_str(), chainID2.c_str(), ylen, d0A, Liden/ylen, TM1);
        printf("%s\n", seqyA);

        printf("# Lali=%d\tRMSD=%.2f\tseqID_ali=%.3f\n",
            n_ali8, rmsd, (n_ali8>0)?Liden/n_ali8:0);

        if (i_opt)
            printf("# User-specified initial alignment: TM=%.5lf\tLali=%4d\trmsd=%.3lf\n", TM_ali, L_ali, rmsd_ali);

        if(a_opt)
            printf("# TM-score=%.5f (normalized by average length of two structures: L=%.1f\td0=%.2f)\n", TM3, (xlen+ylen)*0.5, d0a);

        if(u_opt)
            printf("# TM-score=%.5f (normalized by user-specified L=%.2f\td0=%.2f)\n", TM4, Lnorm_ass, d0u);

        if(d_opt)
            printf("# TM-score=%.5f (scaled by user-specified d0=%.2f\tL=%d)\n", TM5, d0_scale, ylen);

        printf("$$$$\n");
    }
    else if (outfmt_opt==2)
    {
        printf("%s%s\t%s%s\t%.4f\t%.4f\t%.2f\t%4.3f\t%4.3f\t%4.3f\t%d\t%d\t%d",
            xname.c_str(), chainID1.c_str(), yname.c_str(), chainID2.c_str(),
            TM2, TM1, rmsd, Liden/xlen, Liden/ylen, (n_ali8>0)?Liden/n_ali8:0,
            xlen, ylen, n_ali8);
    }
    cout << endl;

    if (strlen(fname_matrix)) output_flexalign_rotation_matrix(
            fname_matrix, tu_vec, t, u);

    if (o_opt==1) output_flexalign_pymol(xname, yname, fname_super, tu_vec,
            t, u, ter_opt, mm_opt, split_opt, mirror_opt, seqM, seqxA, seqyA,
            resi_vec1, resi_vec2, chainID1, chainID2);
    else if (o_opt==2)
        output_flexalign_rasmol(xname, yname, fname_super, tu_vec,
            t, u, ter_opt, mm_opt, split_opt, mirror_opt, seqM, seqxA, seqyA,
            resi_vec1, resi_vec2, chainID1, chainID2,
            xlen, ylen, d0A, n_ali8, rmsd, TM1, Liden);
}

#endif
