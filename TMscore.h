#include "TMalign.h"

void GDT_from_TMscore(double **xa, double **ya, int xlen, int ylen,
    string &seqM, const string &seqxA, const string &seqyA,
    string &seq_scale, int &L_lt_d, double &rmsd_d0_out,
    double GDT_TS_list[4], double GDT_HA_list[4],
    double t[3], double u[3][3], const double d0_out)
{
    int i,j;

    /* prepare seqM */
    string seqM_tmp=seqM;
    seq_scale.clear();
    seq_scale.assign(seqM.size(),' ');
    L_lt_d=0;
    for (i=0;i<seqM.size();i++)
    {
        seq_scale[i]='0'+(1+i)%10;
        if (seqM[i]!=':') seqM_tmp[i]=' ';
        else L_lt_d++;
    }
    seqM=seqM_tmp;
    seqM_tmp.clear();
    return;

    /* extract alignment */
    int L = min(seqxA.size(), seqyA.size());
    int *invmap = new int[ylen+1];
    for (j=0; j<ylen; j++) invmap[j] = -1;

    int i1 = -1;// in C version, index starts from zero, not from one
    int i2 = -1;
    for (int kk1 = 0; kk1 < L; kk1++)
    {
        if (seqxA[kk1] != '-') i1++;
        if (seqyA[kk1] != '-')
        {
            i2++;
            if (i2 >= ylen || i1 >= xlen) kk1 = L;
            else if (seqxA[kk1] != '-') invmap[i2] = i1;
        }
    }

    /* calculate GDT */
    double x1[3];
    double d;
    rmsd_d0_out=0;
    L_lt_d=0;
    for (j=0;j<ylen;j++)
    {
        i=invmap[j];
        if (i==-1) continue;
        transform(t, u, xa[i], x1);
        d=sqrt(dist(x1,ya[j]));

        if (d<d0_out)
        {
            L_lt_d++;
            rmsd_d0_out+=d*d;
        }
        if (d<8)
        {
            GDT_TS_list[3]+=1;
            if (d<4)
            {
                GDT_TS_list[2]+=1;
                GDT_HA_list[3]+=1;
                if (d<2)
                {
                    GDT_TS_list[1]+=1;
                    GDT_HA_list[2]+=1;
                    if (d<1)
                    {
                        GDT_TS_list[0]+=1;
                        GDT_HA_list[1]+=1;
                        if (d<0.5) GDT_HA_list[0]+=1;
                    }
                }
            }
        }
    }
    if (L_lt_d) rmsd_d0_out=sqrt(rmsd_d0_out/L_lt_d);
    for (i=0;i<4;i++)
    {
        GDT_TS_list[i]/=ylen;
        GDT_HA_list[i]/=ylen;
    }
}

void output_TMscore_results(
    const string xname, const string yname,
    const char *chainID1, const char *chainID2,
    const int xlen, const int ylen, double t[3], double u[3][3],
    const double TM1, const double TM2,
    const double TM3, const double TM4, const double TM5,
    const double rmsd, const double d0_out,
    const char *seqM, const char *seqxA, const char *seqyA, const double Liden,
    const int n_ali8, const int L_ali,
    const double TM_ali, const double rmsd_ali, const double TM_0,
    const double d0_0, const double d0A, const double d0B,
    const double Lnorm_ass, const double d0_scale, 
    const double d0a, const double d0u, const char* fname_matrix,
    const int outfmt_opt, const int ter_opt, const char *fname_super,
    const int a_opt, const bool u_opt, const bool d_opt, const int mirror_opt,
    const char *seq_scale, const int L_lt_d, const double rmsd_d0_out,
    double GDT_TS_list[4], double GDT_HA_list[4])
{
    if (outfmt_opt<=0)
    {
        printf("\nStructure1: %s%s    Length=%5d\n",
            xname.c_str(), chainID1, xlen);
        printf("Structure2: %s%s    Length=%5d (by which all scores are normalized)\n",
            yname.c_str(), chainID2, ylen);

        printf("Number of residues in common=%5d\n", n_ali8);
        printf("RMSD of  the common residues=%9.3f\n\n", rmsd);
        printf("TM-score    = %6.4f  (d0= %.2f)\n", TM1, d0A);

        //double gdt_ts_score=0;
        //double gdt_ha_score=0;
        //for (int i=0;i<4;i++)
        //{
            //gdt_ts_score+=GDT_TS_list[i];
            //gdt_ha_score+=GDT_HA_list[i];
        //}
        //gdt_ts_score/=4.;
        //gdt_ha_score/=4.;
        //printf("GDT-TS-score= %6.4f %(d<1)=%6.4f %(d<2)=%6.4f %(d<4)=%6.4f %(d<8)=%6.4f\n",
            //gdt_ts_score, GDT_TS_list[0], GDT_TS_list[1],
                          //GDT_TS_list[2], GDT_TS_list[3]);
        //printf("GDT-HA-score= %6.4f %(d<0.5)=%6.4f %(d<1)=%6.4f %(d<2)=%6.4f %(d<4)=%6.4f\n",
            //gdt_ha_score, GDT_HA_list[0], GDT_HA_list[1],
                          //GDT_HA_list[2], GDT_HA_list[3]);

        if (a_opt==1)
            printf("TM-score    = %5.4f  (if normalized by average length of two structures, i.e., LN= %.1f, d0= %.2f)\n", TM3, (xlen+ylen)*0.5, d0a);
        if (u_opt)
            printf("TM-score    = %5.4f  (if normalized by user-specified LN=%.2f and d0=%.2f)\n", TM4, Lnorm_ass, d0u);
        if (d_opt)
            printf("TM-score    = %5.5f  (if scaled by user-specified d0= %.2f, and LN= %d)\n", TM5, d0_scale, ylen);
    

        printf("\n -------- rotation matrix to rotate Chain-1 to Chain-2 ------\n");
        printf(" i          t(i)         u(i,1)         u(i,2)         u(i,3)\n");
        printf(" 1 %17.10f %14.10f %14.10f %14.10f\n",t[0],u[0][0],u[0][1],u[0][2]);
        printf(" 2 %17.10f %14.10f %14.10f %14.10f\n",t[1],u[1][0],u[1][1],u[1][2]);
        printf(" 3 %17.10f %14.10f %14.10f %14.10f\n",t[2],u[2][0],u[2][1],u[2][2]);

        //output alignment
        printf("\nSuperposition in the TM-score: Length(d<%3.1f)= %d\n", d0_out, L_lt_d);
        //printf("\nSuperposition in the TM-score: Length(d<%3.1f)= %d  RMSD=%6.2f\n", d0_out, L_lt_d, rmsd_d0_out);
        printf("(\":\" denotes the residue pairs of distance <%4.1f Angstrom)\n", d0_out);
        printf("%s\n", seqxA);
        printf("%s\n", seqM);
        printf("%s\n", seqyA);
        printf("%s\n", seq_scale);
    }
    else if (outfmt_opt==1)
    {
        printf(">%s%s\tL=%d\td0=%.2f\tseqID=%.3f\tTM-score=%.5f\n",
            xname.c_str(), chainID1, xlen, d0B, Liden/xlen, TM2);
        printf("%s\n", seqxA);
        printf(">%s%s\tL=%d\td0=%.2f\tseqID=%.3f\tTM-score=%.5f\n",
            yname.c_str(), chainID2, ylen, d0A, Liden/ylen, TM1);
        printf("%s\n", seqyA);

        printf("# Lali=%d\tRMSD=%.2f\tseqID_ali=%.3f\n",
            n_ali8, rmsd, (n_ali8>0)?Liden/n_ali8:0);

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
            xname.c_str(), chainID1, yname.c_str(), chainID2, TM2, TM1, rmsd,
            Liden/xlen, Liden/ylen, (n_ali8>0)?Liden/n_ali8:0,
            xlen, ylen, n_ali8);
    }
    cout << endl;

    if (strlen(fname_matrix)) 
        output_rotation_matrix(fname_matrix, t, u);
    if (strlen(fname_super))
        output_superpose(xname, fname_super, t, u, ter_opt, mirror_opt);
}
