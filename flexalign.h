/* Functions for the core TMalign algorithm, including the entry function
 * flexalign_main */
#ifndef flexalign_h
#define flexalign_h 1

#include "TMalign.h"

void t_u2tu(double t0[3], double u0[3][3], vector<double> &tu_tmp)
{
    int i, j, k;
    for (i = 0; i < 3; i++)
        tu_tmp[i] = t0[i];
    k = 3;
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
        {
            tu_tmp[k] = u0[i][j];
            k++;
        }
}

void tu2t_u(vector<double> tu_tmp, double t0[3], double u0[3][3])
{
    int i, j, k;
    for (i = 0; i < 3; i++)
        t0[i] = tu_tmp[i];
    k = 3;
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
        {
            u0[i][j] = tu_tmp[k];
            k++;
        }
}

void aln2invmap(const string &seqxA, const string &seqyA, int *invmap)
{
    int i, j, r;
    int ylen = 0;
    for (r = 0; r < seqyA.size(); r++)
        ylen += seqyA[r] != '-';
    for (j = 0; j < ylen; j++)
        invmap[j] = -1;

    i = j = -1;
    for (r = 0; r < seqxA.size(); r++)
    {
        i += seqxA[r] != '-';
        j += seqyA[r] != '-';
        if (seqxA[r] != '-' && seqyA[r] != '-')
            invmap[j] = i;
    }
}

int flexalign_main(double **xa, double **ya,
                   const char *seqx, const char *seqy, const char *secx, const char *secy,
                   double t0[3], double u0[3][3], vector<vector<double>> &tu_vec,
                   double &TM1, double &TM2, double &TM3, double &TM4, double &TM5,
                   double &d0_0, double &TM_0,
                   double &d0A, double &d0B, double &d0u, double &d0a, double &d0_out,
                   string &seqM, string &seqxA, string &seqyA, vector<double> &do_vec,
                   double &rmsd0, int &L_ali, double &Liden,
                   double &TM_ali, double &rmsd_ali, int &n_ali, int &n_ali8,
                   const int xlen, const int ylen,
                   const vector<string> sequence, const double Lnorm_ass,
                   const double d0_scale, const int i_opt, const int a_opt,
                   const bool u_opt, const bool d_opt, const bool fast_opt,
                   const int mol_type, const int hinge_opt, const int ss_opt)
{
    vector<double> tu_tmp(12, 0);
    int round2 = tu_vec.size();
    if (round2 == 0)
    {
        TMalign_main(xa, ya, seqx, seqy, secx, secy, t0, u0,
                     TM1, TM2, TM3, TM4, TM5, d0_0, TM_0,
                     d0A, d0B, d0u, d0a, d0_out, seqM, seqxA, seqyA, do_vec,
                     rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                     xlen, ylen, sequence, Lnorm_ass,
                     d0_scale, i_opt, a_opt, u_opt, d_opt, fast_opt, mol_type, -1, ss_opt);

        t_u2tu(t0, u0, tu_tmp);
        tu_vec.push_back(tu_tmp);
    }

    int i, j, r;
    int *invmap = new int[ylen + 1];
    for (j = 0; j < ylen + 1; j++)
        invmap[j] = -1;
    double **xt;
    NewArray(&xt, xlen, 3);
    do_rotation(xa, xt, xlen, t0, u0);

    TM1 = TM2 = TM3 = TM4 = TM5 = rmsd0 = 0;
    seqM = "";
    seqxA = "";
    seqyA = "";
    n_ali = n_ali8 = 0;
    se_main(xt, ya, seqx, seqy, TM1, TM2, TM3, TM4, TM5, d0_0, TM_0,
            d0A, d0B, d0u, d0a, d0_out, seqM, seqxA, seqyA, do_vec,
            rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
            xlen, ylen, sequence, Lnorm_ass, d0_scale, i_opt,
            a_opt, u_opt, d_opt, mol_type, 0, invmap, 1);
    if (round2)
    {
        /* aligned structure A vs unaligned structure B */
        int xlen_h = n_ali8;
        int ylen_h = ylen - n_ali8;
        char *seqx_h = new char[xlen + 1];
        char *seqy_h = new char[ylen + 1];
        char *secx_h = new char[xlen + 1];
        char *secy_h = new char[ylen + 1];
        seqx_h[xlen] = seqy_h[ylen] = 0;
        secx_h[xlen] = secy_h[ylen] = 0;
        double **xa_h, **ya_h;
        NewArray(&xa_h, xlen, 3);
        NewArray(&ya_h, ylen, 3);

        int r1, r2;
        i = j = -1;
        r1 = r2 = 0;
        for (r = 0; r < seqxA.size(); r++)
        {
            i += (seqxA[r] != '-');
            j += (seqyA[r] != '-');
            if (seqxA[r] != '-' && seqyA[r] != '-')
            {
                seqx_h[r1] = seqx[i];
                secx_h[r1] = secx[i];
                xa_h[r1][0] = xa[i][0];
                xa_h[r1][1] = xa[i][1];
                xa_h[r1][2] = xa[i][2];
                r1++;
            }
            if (seqxA[r] == '-')
            {
                seqy_h[r2] = seqx[j];
                secy_h[r2] = secx[j];
                ya_h[r2][0] = ya[j][0];
                ya_h[r2][1] = ya[j][1];
                ya_h[r2][2] = ya[j][2];
                r2++;
            }
        }

        double TM1_h, TM2_h;
        double TM3_h, TM4_h, TM5_h; // for a_opt, u_opt, d_opt
        double d0_0_h, TM_0_h;
        double d0A_h, d0B_h, d0u_h, d0a_h;
        double d0_out_h = 5.0;
        string seqM_h, seqxA_h, seqyA_h; // for output alignment
        double rmsd0_h = 0.0;
        int L_ali_h = 0; // Aligned length in standard_TMscore
        double Liden_h = 0;
        double TM_ali_h, rmsd_ali_h; // TMscore and rmsd in standard_TMscore
        int n_ali_h = 0;
        int n_ali8_h = 0;

        TMalign_main(xa_h, ya_h, seqx_h, seqy_h, secx_h, secy_h, t0, u0,
                     TM1_h, TM2_h, TM3_h, TM4_h, TM5_h, d0_0_h, TM_0_h, d0A_h, d0B_h,
                     d0u_h, d0a_h, d0_out_h, seqM_h, seqxA_h, seqyA_h, do_vec,
                     rmsd0_h, L_ali_h, Liden_h, TM_ali_h, rmsd_ali_h, n_ali_h, n_ali8_h,
                     xlen_h, ylen_h, sequence, Lnorm_ass,
                     d0_scale, i_opt, a_opt, u_opt, d_opt, fast_opt, mol_type, -1, ss_opt);

        do_rotation(xa, xt, xlen, t0, u0);
        t_u2tu(t0, u0, tu_vec[0]);

        int *invmap_h = new int[ylen + 1];
        for (j = 0; j < ylen + 1; j++)
            invmap_h[j] = -1;
        TM1_h = TM2_h = TM3_h = TM4_h = TM5_h = rmsd0_h = 0;
        seqM_h = "";
        seqxA_h = "";
        seqyA_h = "";
        n_ali_h = n_ali8_h = 0;
        se_main(xt, ya, seqx, seqy, TM1_h, TM2_h, TM3_h, TM4_h, TM5_h, d0_0,
                TM_0, d0A, d0B, d0u, d0a, d0_out, seqM_h, seqxA_h, seqyA_h, do_vec,
                rmsd0_h, L_ali, Liden, TM_ali, rmsd_ali, n_ali_h, n_ali8_h,
                xlen, ylen, sequence, Lnorm_ass, d0_scale, i_opt,
                a_opt, u_opt, d_opt, mol_type, 0, invmap_h, 1);

        /* unaligned structure A vs aligned structure B */
        xlen_h = xlen - n_ali8;
        ylen_h = n_ali8;

        i = j = -1;
        r1 = r2 = 0;
        for (r = 0; r < seqxA.size(); r++)
        {
            i += (seqxA[r] != '-');
            j += (seqyA[r] != '-');
            if (seqyA[r] == '-')
            {
                seqx_h[r1] = seqx[i];
                secx_h[r1] = secx[i];
                xa_h[r1][0] = xa[i][0];
                xa_h[r1][1] = xa[i][1];
                xa_h[r1][2] = xa[i][2];
                r1++;
            }
            if (seqxA[r] != '-' && seqyA[r] != '-')
            {
                seqy_h[r2] = seqx[j];
                secy_h[r2] = secx[j];
                ya_h[r2][0] = ya[j][0];
                ya_h[r2][1] = ya[j][1];
                ya_h[r2][2] = ya[j][2];
                r2++;
            }
        }

        d0_out_h = 5.0;
        L_ali_h = Liden_h = 0;
        TM1 = TM2 = TM3 = TM4 = TM5 = rmsd0 = 0;
        seqM = "";
        seqxA = "";
        seqyA = "";
        n_ali = n_ali8 = 0;

        TMalign_main(xa_h, ya_h, seqx_h, seqy_h, secx_h, secy_h, t0, u0,
                     TM1, TM2, TM3, TM4, TM5, d0_0_h, TM_0_h, d0A_h, d0B_h,
                     d0u_h, d0a_h, d0_out_h, seqM, seqxA, seqyA, do_vec,
                     rmsd0, L_ali_h, Liden_h, TM_ali_h, rmsd_ali_h, n_ali, n_ali8,
                     xlen_h, ylen_h, sequence, Lnorm_ass,
                     d0_scale, i_opt, a_opt, u_opt, d_opt, fast_opt, mol_type, -1, ss_opt);

        do_rotation(xa, xt, xlen, t0, u0);

        for (j = 0; j < ylen + 1; j++)
            invmap[j] = -1;
        TM1 = TM2 = TM3 = TM4 = TM5 = rmsd0 = 0;
        seqM = "";
        seqxA = "";
        seqyA = "";
        n_ali = n_ali8 = 0;
        se_main(xt, ya, seqx, seqy, TM1, TM2, TM3, TM4, TM5, d0_0,
                TM_0, d0A, d0B, d0u, d0a, d0_out, seqM, seqxA, seqyA, do_vec,
                rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                xlen, ylen, sequence, Lnorm_ass, d0_scale, i_opt,
                a_opt, u_opt, d_opt, mol_type, 0, invmap, 1);

        double TM_h = (TM1_h > TM2_h) ? TM1_h : TM2_h;
        double TM = (TM1 > TM2) ? TM1 : TM2;
        if (TM_h > TM)
        {
            TM1 = TM1_h;
            TM2 = TM2_h;
            TM3 = TM3_h;
            TM4 = TM4_h;
            TM5 = TM5_h;
            seqM = seqM_h;
            seqxA = seqxA_h;
            seqyA = seqyA_h;
            rmsd0 = rmsd0_h;
            n_ali = n_ali_h;
            n_ali8 = n_ali8_h;
            for (j = 0; j < ylen + 1; j++)
                invmap[j] = invmap_h[j];
        }
        else
            t_u2tu(t0, u0, tu_vec[0]);

        /* clean up */
        delete[] invmap_h;
        DeleteArray(&xa_h, xlen);
        DeleteArray(&ya_h, ylen);
        seqM_h.clear();
        seqxA_h.clear();
        seqyA_h.clear();
        delete[] seqx_h;
        delete[] secx_h;
        delete[] seqy_h;
        delete[] secy_h;
    }
    for (r = 0; r < seqM.size(); r++)
        if (seqM[r] == '1')
            seqM[r] = '0';

    int minlen = min(xlen, ylen);
    int hinge;
    for (hinge = 0; hinge < hinge_opt; hinge++)
    {
        if (minlen - n_ali8 < 5)
            break;
        int xlen_h = xlen - n_ali8;
        int ylen_h = ylen - n_ali8;
        char *seqx_h = new char[xlen_h + 1];
        char *seqy_h = new char[ylen_h + 1];
        char *secx_h = new char[xlen_h + 1];
        char *secy_h = new char[ylen_h + 1];
        seqx_h[xlen_h] = seqy_h[ylen_h] = 0;
        secx_h[xlen_h] = secy_h[ylen_h] = 0;
        double **xa_h, **ya_h;
        NewArray(&xa_h, xlen_h, 3);
        NewArray(&ya_h, ylen_h, 3);
        vector<int> r1toi(xlen_h, 0);
        vector<int> r2toj(ylen_h, 0);

        int r1, r2;
        i = j = -1;
        r1 = r2 = 0;
        for (r = 0; r < seqxA.size(); r++)
        {
            i += (seqxA[r] != '-');
            j += (seqyA[r] != '-');
            if (seqyA[r] == '-')
            {
                seqx_h[r1] = seqx[i];
                secx_h[r1] = secx[i];
                xa_h[r1][0] = xa[i][0];
                xa_h[r1][1] = xa[i][1];
                xa_h[r1][2] = xa[i][2];
                r1toi[r1] = i;
                r1++;
            }
            if (seqxA[r] == '-')
            {
                seqy_h[r2] = seqx[j];
                secy_h[r2] = secx[j];
                ya_h[r2][0] = ya[j][0];
                ya_h[r2][1] = ya[j][1];
                ya_h[r2][2] = ya[j][2];
                r2toj[r2] = j;
                r2++;
            }
        }

        double TM1_h, TM2_h;
        double TM3_h, TM4_h, TM5_h; // for a_opt, u_opt, d_opt
        double d0_0_h, TM_0_h;
        double d0A_h, d0B_h, d0u_h, d0a_h;
        double d0_out_h = 5.0;
        string seqM_h, seqxA_h, seqyA_h; // for output alignment
        double rmsd0_h = 0.0;
        int L_ali_h = 0; // Aligned length in standard_TMscore
        double Liden_h = 0;
        double TM_ali_h, rmsd_ali_h; // TMscore and rmsd in standard_TMscore
        int n_ali_h = 0;
        int n_ali8_h = 0;

        TMalign_main(xa_h, ya_h, seqx_h, seqy_h, secx_h, secy_h, t0, u0,
                     TM1_h, TM2_h, TM3_h, TM4_h, TM5_h, d0_0_h, TM_0_h, d0A_h, d0B_h,
                     d0u_h, d0a_h, d0_out_h, seqM_h, seqxA_h, seqyA_h, do_vec,
                     rmsd0_h, L_ali_h, Liden_h, TM_ali_h, rmsd_ali_h, n_ali_h, n_ali8_h,
                     xlen_h, ylen_h, sequence, Lnorm_ass,
                     d0_scale, i_opt, a_opt, u_opt, d_opt, fast_opt, mol_type, -1, ss_opt);

        do_rotation(xa, xt, xlen, t0, u0);

        TM1_h = TM1;
        TM2_h = TM2;
        TM3_h = TM3;
        TM4_h = TM4;
        TM5_h = TM5;
        seqM_h = seqM;
        seqxA_h = seqxA;
        seqyA_h = seqyA;
        rmsd0_h = rmsd0;
        n_ali_h = n_ali;
        n_ali8_h = n_ali8;
        int *invmap_h = new int[ylen + 1];
        for (j = 0; j < ylen + 1; j++)
            invmap_h[j] = invmap[j];
        se_main(xt, ya, seqx, seqy, TM1_h, TM2_h, TM3_h, TM4_h, TM5_h, d0_0, TM_0,
                d0A, d0B, d0u, d0a, d0_out, seqM_h, seqxA_h, seqyA_h, do_vec,
                rmsd0_h, L_ali, Liden, TM_ali, rmsd_ali, n_ali_h, n_ali8_h,
                xlen, ylen, sequence, Lnorm_ass, d0_scale, i_opt,
                a_opt, u_opt, d_opt, mol_type, 0, invmap_h, hinge + 1);
        int new_ali = 0;
        for (r = 0; r < seqM_h.size(); r++)
            new_ali += (seqM_h[r] == hinge + '1');
        if (n_ali8_h - n_ali8 < 5)
            new_ali = 0;
        if (new_ali >= 5)
        {
            TM1 = TM1_h;
            TM2 = TM2_h;
            TM3 = TM3_h;
            TM4 = TM4_h;
            TM5 = TM5_h;
            seqM = seqM_h;
            seqxA = seqxA_h;
            seqyA = seqyA_h;
            rmsd0 = rmsd0_h;
            n_ali = n_ali_h;
            n_ali8 = n_ali8_h;
            t_u2tu(t0, u0, tu_tmp);
            tu_vec.push_back(tu_tmp);
            for (j = 0; j < ylen + 1; j++)
                invmap[j] = invmap_h[j];
            // cout<<">hinge="<<hinge<<'\n'
            //<<seqxA<<'\n'<<seqM<<'\n'<<seqyA<<endl;
            // for (j=0;j<ylen;j++) if ((i=invmap[j])>=0) cout<<"("<<i<<","<<j<<")";
            // cout<<endl;
        }

        /* clean up */
        delete[] invmap_h;
        DeleteArray(&xa_h, xlen_h);
        DeleteArray(&ya_h, ylen_h);
        r1toi.clear();
        r2toj.clear();
        seqM_h.clear();
        seqxA_h.clear();
        seqyA_h.clear();
        delete[] seqx_h;
        delete[] secx_h;
        delete[] seqy_h;
        delete[] secy_h;
        if (new_ali < 5)
            break;
    }

    if (tu_vec.size() <= 1)
    {
        DeleteArray(&xt, xlen);
        delete[] invmap;
        return tu_vec.size();
    }

    /* re-derive alignment based on tu_vec */
    vector<char> seqM_char(ylen, ' ');
    vector<double> di_vec(ylen, -1);
    double d;
    for (hinge = tu_vec.size() - 1; hinge >= 0; hinge--)
    {
        tu2t_u(tu_vec[hinge], t0, u0);
        do_rotation(xa, xt, xlen, t0, u0);
        for (j = 0; j < ylen; j++)
        {
            i = invmap[j];
            if (i < 0)
                continue;
            d = sqrt(dist(xt[i], ya[j]));
            if (di_vec[j] < 0 || d <= di_vec[j])
            {
                di_vec[j] = d;
                seqM_char[j] = hinge + '0';
            }
        }
    }
    j = -1;
    for (r = 0; r < seqM.size(); r++)
    {
        if (seqyA[r] == '-')
            continue;
        j++;
        seqM[r] = seqM_char[j];
    }

    /* smooth out AFP assignment: remove singleton insert */
    for (hinge = tu_vec.size() - 1; hinge >= 0; hinge--)
    {
        j = -1;
        for (r = 0; r < seqM.size(); r++)
        {
            if (seqyA[r] == '-')
                continue;
            j++;
            if (seqM_char[j] != hinge + '0')
                continue;
            if (r < seqM.size() - 1 && (seqM[r + 1] == hinge + '0' || seqM[r + 1] == ' '))
                continue;
            if (r > 0 && (seqM[r - 1] == hinge + '0' || seqM[r - 1] == ' '))
                continue;
            if (r < seqM.size() - 1 && r > 0 && seqM[r - 1] != seqM[r + 1])
                continue;
            if (r > 0)
                seqM[r] = seqM_char[j] = seqM[r - 1];
            else
                seqM[r] = seqM_char[j] = seqM[r + 1];
        }
    }
    /* smooth out AFP assignment: remove singleton at the end of fragment */
    char left_hinge = ' ';
    char right_hinge = ' ';
    for (hinge = tu_vec.size() - 1; hinge >= 0; hinge--)
    {
        j = -1;
        for (r = 0; r < seqM.size(); r++)
        {
            if (seqyA[r] == '-')
                continue;
            j++;
            if (seqM[r] != hinge + '0')
                continue;
            if (r > 0 && seqM[r - 1] == ' ' && r < seqM.size() - 1 && seqM[r + 1] == ' ')
                continue;

            left_hinge = ' ';
            for (i = r - 1; i >= 0; i--)
            {
                if (seqM[i] == ' ')
                    continue;
                left_hinge = seqM[i];
                break;
            }
            if (left_hinge == hinge + '0')
                continue;

            right_hinge = ' ';
            for (i = r + 1; i < seqM.size(); i++)
            {
                if (seqM[i] == ' ')
                    continue;
                right_hinge = seqM[i];
                break;
            }
            if (right_hinge == hinge + '0')
                continue;
            if (left_hinge != right_hinge && left_hinge != ' ' && right_hinge != ' ')
                continue;

            if (right_hinge != ' ')
                seqM[r] = seqM_char[j] = right_hinge;
            else if (left_hinge != ' ')
                seqM[r] = seqM_char[j] = left_hinge;
        }
    }
    /* smooth out AFP assignment: remove dimer insert */
    for (hinge = tu_vec.size() - 1; hinge >= 0; hinge--)
    {
        j = -1;
        for (r = 0; r < seqM.size() - 1; r++)
        {
            if (seqyA[r] == '-')
                continue;
            j++;
            if (seqM[r] != hinge + '0' || seqM[r + 1] != hinge + '0')
                continue;

            if (r < seqM.size() - 2 && (seqM[r + 2] == ' ' || seqM[r + 2] == hinge + '0'))
                continue;
            if (r > 0 && (seqM[r - 1] == ' ' || seqM[r - 1] == hinge + '0'))
                continue;
            if (r < seqM.size() - 2 && r > 0 && seqM[r - 1] != seqM[r + 2])
                continue;

            if (r > 0)
                seqM[r] = seqM_char[j] = seqM[r + 1] = seqM_char[j + 1] = seqM[r - 1];
            else
                seqM[r] = seqM_char[j] = seqM[r + 1] = seqM_char[j + 1] = seqM[r + 2];
        }
    }
    /* smooth out AFP assignment: remove disconnected singleton */
    int i1, i2;
    for (hinge = tu_vec.size() - 1; hinge >= 0; hinge--)
    {
        j = -1;
        for (r = 0; r < seqM.size(); r++)
        {
            if (seqyA[r] == '-')
                continue;
            j++;
            if (seqM[r] != hinge + '0')
                continue;

            left_hinge = ' ';
            for (i = r - 1; i >= 0; i--)
            {
                if (seqM[i] == ' ')
                    continue;
                left_hinge = seqM[i];
                i1 = (r - i);
                break;
            }
            if (left_hinge == hinge + '0')
                continue;

            right_hinge = ' ';
            for (i = r + 1; i < seqM.size(); i++)
            {
                if (seqM[i] == ' ')
                    continue;
                right_hinge = seqM[i];
                i2 = (i - r);
                break;
            }
            if (right_hinge == hinge + '0')
                continue;

            if (right_hinge == ' ')
                seqM[r] = seqM_char[j] = left_hinge;
            else if (left_hinge == ' ')
                seqM[r] = seqM_char[j] = right_hinge;
            else
            {
                if (i1 < i2)
                    seqM[r] = seqM_char[j] = left_hinge;
                else
                    seqM[r] = seqM_char[j] = right_hinge;
            }
        }
    }

    /* recalculate all scores */
    for (hinge = tu_vec.size() - 1; hinge >= 0; hinge--)
    {
        tu2t_u(tu_vec[hinge], t0, u0);
        do_rotation(xa, xt, xlen, t0, u0);
        for (j = 0; j < ylen; j++)
        {
            i = invmap[j];
            if (i < 0)
                continue;
            if (seqM_char[j] != hinge + '0')
                continue;
            d = sqrt(dist(xt[i], ya[j]));
            if (di_vec[j] < 0 || d <= di_vec[j])
            {
                di_vec[j] = d;
                seqM_char[j] = hinge + '0';
            }
        }
    }
    rmsd0 = TM1 = TM2 = TM3 = TM4 = TM5 = 0;
    Liden = 0;
    for (r = 0; r < seqM.size(); r++)
        if (seqM[r] != ' ')
            Liden += seqxA[r] == seqyA[r];
    for (j = 0; j < ylen; j++)
    {
        i = invmap[j];
        if (i < 0)
            continue;
        {
            d = di_vec[j];
            TM2 += 1 / (1 + (d / d0B) * (d / d0B)); // chain_1
            TM1 += 1 / (1 + (d / d0A) * (d / d0A)); // chain_2
            if (a_opt)
                TM3 += 1 / (1 + (d / d0a) * (d / d0a)); // -a
            if (u_opt)
                TM4 += 1 / (1 + (d / d0u) * (d / d0u)); // -u
            if (d_opt)
                TM5 += 1 / (1 + (d / d0_scale) * (d / d0_scale)); // -d
            rmsd0 += d * d;
        }
    }
    TM2 /= xlen;
    TM1 /= ylen;
    TM3 /= (xlen + ylen) * 0.5;
    TM4 /= Lnorm_ass;
    TM5 /= ylen;
    if (n_ali8)
        rmsd0 = sqrt(rmsd0 / n_ali8);
    for (hinge = tu_vec.size() - 1; hinge > 0; hinge--)
    {
        int afp_len = 0;
        for (r = 0; r < seqM.size(); r++)
            afp_len += seqM[r] == hinge + '0';
        if (afp_len)
            break;
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
void output_flexalign_rotation_matrix(const char *fname_matrix,
                                      const vector<vector<double>> &tu_vec, double t[3], double u[3][3])
{
    stringstream ss;
    char dest[1000];
    for (int hinge = 0; hinge < tu_vec.size(); hinge++)
    {
        tu2t_u(tu_vec[hinge], t, u);
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
    if (strcmp(fname_matrix, (char *)("-")) == 0)
        cout << ss.str();
    else
    {
        fstream fout;
        fout.open(fname_matrix, ios::out | ios::trunc);
        if (fout)
        {
            fout << ss.str();
            fout.close();
        }
        else
            cout << "Open file to output rotation matrix fail.\n";
    }
    ss.str(string());
}

void output_flexalign_rasmol(const string xname, const string yname,
                             const string fname_super, const vector<vector<double>> &tu_vec,
                             double t[3], double u[3][3], const int ter_opt,
                             const int mm_opt, const int split_opt, const int mirror_opt,
                             const char *seqM, const char *seqxA, const char *seqyA,
                             const vector<string> &resi_vec1, const vector<string> &resi_vec2,
                             const string chainID1, const string chainID2,
                             const int xlen, const int ylen, const double d0A, const int n_ali8,
                             const double rmsd, const double TM1, const double Liden)
{
    stringstream buf;
    stringstream buf_all;
    stringstream buf_atm;
    stringstream buf_all_atm;
    stringstream buf_all_atm_lig;
    // stringstream buf_pdb;
    stringstream buf_tm;
    string line;
    double x[3];    // before transform
    double x1[3];   // after transform
    bool after_ter; // true if passed the "TER" line in PDB
    string asym_id; // chain ID

    map<string, int> resi2hinge_dict;
    int r, i, j;
    j = -1;
    char hinge_char = 0;
    int ali_len = strlen(seqM);
    for (r = 0; r < strlen(seqxA); r++)
    {
        if (seqxA[r] == '-')
            continue;
        j++;
        hinge_char = seqM[r];
        if (hinge_char == ' ')
        {
            for (i = 1; i < ali_len; i++)
            {
                if (r - i >= 0 && seqM[r - i] != ' ')
                    hinge_char = seqM[r - i];
                else if (r + i < xlen && seqM[r + i] != ' ')
                    hinge_char = seqM[r + i];
                if (hinge_char != ' ')
                    break;
            }
        }
        int hinge_idx = 0;
        if (hinge_char >= '0' && hinge_char <= '9')
        {
            hinge_idx = hinge_char - '0';
        }
        else if (hinge_char >= 'a' && hinge_char <= 'z')
        {
            hinge_idx = hinge_char - 'a' + 10;
        }
        else if (hinge_char >= 'A' && hinge_char <= 'Z')
        {
            hinge_idx = hinge_char - 'A' + 36;
        }
        resi2hinge_dict[resi_vec1[j]] = hinge_idx;
    }
    string resi = resi_vec1[0];
    int read_resi = resi.size() - 4;

    buf_tm << "REMARK US-align"
           << "\nREMARK Structure 1:" << setw(11) << left << xname + chainID1 << " Size= " << xlen
           << "\nREMARK Structure 2:" << setw(11) << yname + chainID2 << right << " Size= " << ylen
           << " (TM-score is normalized by " << setw(4) << ylen << ", d0="
           << setiosflags(ios::fixed) << setprecision(2) << setw(6) << d0A << ")"
           << "\nREMARK Aligned length=" << setw(4) << n_ali8 << ", RMSD="
           << setw(6) << setiosflags(ios::fixed) << setprecision(2) << rmsd
           << ", TM-score=" << setw(7) << setiosflags(ios::fixed) << setprecision(5) << TM1
           << ", ID=" << setw(5) << setiosflags(ios::fixed) << setprecision(3)
           << ((n_ali8 > 0) ? Liden / n_ali8 : 0) << endl;
    string rasmol_CA_header = "load inline\nselect *A\nwireframe .45\nselect *B\nwireframe .20\nselect all\ncolor white\n";
    string rasmol_cartoon_header = "load inline\nselect all\ncartoon\nselect *A\ncolor blue\nselect *B\ncolor red\nselect ligand\nwireframe 0.25\nselect solvent\nspacefill 0.25\nselect all\nexit\n" + buf_tm.str();
    if (!mm_opt)
        buf << rasmol_CA_header;
    buf_all << rasmol_CA_header;
    if (!mm_opt)
        buf_atm << rasmol_cartoon_header;
    buf_all_atm << rasmol_cartoon_header;
    buf_all_atm_lig << rasmol_cartoon_header;

    /* selecting chains for -mol */
    string chain1_sele;
    string chain2_sele;
    if (!mm_opt)
    {
        if (split_opt == 2 && ter_opt >= 1) // align one chain from model 1
        {
            chain1_sele = chainID1.substr(1);
            chain2_sele = chainID2.substr(1);
        }
        else if (split_opt == 2 && ter_opt == 0) // align one chain from each model
        {
            for (i = 1; i < chainID1.size(); i++)
                if (chainID1[i] == ',')
                    break;
            chain1_sele = chainID1.substr(i + 1);
            for (i = 1; i < chainID2.size(); i++)
                if (chainID2[i] == ',')
                    break;
            chain2_sele = chainID2.substr(i + 1);
        }
    }

    /* for PDBx/mmCIF only */
    map<string, int> _atom_site;
    int atom_site_pos;
    vector<string> line_vec;
    string atom;        // 4-character atom name
    string AA;          // 3-character residue name
    string inscode;     // 1-character insertion code
    string model_index; // model index
    bool is_mmcif = false;

    /* used for CONECT record of chain1 */
    int ca_idx1 = 0;  // all CA atoms
    int lig_idx1 = 0; // all atoms
    vector<int> idx_vec;

    /* used for CONECT record of chain2 */
    int ca_idx2 = 0;  // all CA atoms
    int lig_idx2 = 0; // all atoms

    /* extract aligned region */
    vector<string> resi_aln1;
    vector<string> resi_aln2;
    int i1 = -1;
    int i2 = -1;
    if (!mm_opt)
    {
        for (i = 0; i < strlen(seqM); i++)
        {
            i1 += (seqxA[i] != '-');
            i2 += (seqyA[i] != '-');
            if (seqM[i] == ' ')
                continue;
            resi_aln1.push_back(resi_vec1[i1].substr(0, 4));
            resi_aln2.push_back(resi_vec2[i2].substr(0, 4));
            if (seqM[i] != ':')
                continue;
            buf << "select " << resi_aln1.back() << ":A,"
                << resi_aln2.back() << ":B\ncolor red\n";
            buf_all << "select " << resi_aln1.back() << ":A,"
                    << resi_aln2.back() << ":B\ncolor red\n";
        }
        buf << "select all\nexit\n"
            << buf_tm.str();
    }
    buf_all << "select all\nexit\n"
            << buf_tm.str();

    ifstream fin;
    /* read first file */
    after_ter = false;
    asym_id = "";
    fin.open(xname.c_str());
    int hinge = 0;
    while (fin.good())
    {
        getline(fin, line);
        if (ter_opt >= 3 && line.compare(0, 3, "TER") == 0)
            after_ter = true;
        if (is_mmcif == false && line.size() >= 54 &&
            (line.compare(0, 6, "ATOM  ") == 0 ||
             line.compare(0, 6, "HETATM") == 0)) // PDB format
        {
            if (line[16] != 'A' && line[16] != ' ')
                continue;
            x[0] = atof(line.substr(30, 8).c_str());
            x[1] = atof(line.substr(38, 8).c_str());
            x[2] = atof(line.substr(46, 8).c_str());
            if (mirror_opt)
                x[2] = -x[2];
            if (read_resi == 1)
                resi = line.substr(22, 5);
            else
                resi = line.substr(22, 5) + line[21];
            hinge = 0;
            if (resi2hinge_dict.count(resi))
                hinge = resi2hinge_dict[resi];
            tu2t_u(tu_vec[hinge], t, u);
            transform(t, u, x, x1);
            // buf_pdb<<line.substr(0,30)<<setiosflags(ios::fixed)
            //<<setprecision(3)
            //<<setw(8)<<x1[0] <<setw(8)<<x1[1] <<setw(8)<<x1[2]
            //<<line.substr(54)<<'\n';

            if (after_ter && line.compare(0, 6, "ATOM  ") == 0)
                continue;
            lig_idx1++;
            buf_all_atm_lig << line.substr(0, 6) << setw(5) << lig_idx1
                            << line.substr(11, 9) << " A" << line.substr(22, 8)
                            << setiosflags(ios::fixed) << setprecision(3)
                            << setw(8) << x1[0] << setw(8) << x1[1] << setw(8) << x1[2] << '\n';
            if (chain1_sele.size() && line[21] != chain1_sele[0])
                continue;
            if (after_ter || line.compare(0, 6, "ATOM  "))
                continue;
            if (ter_opt >= 2)
            {
                if (ca_idx1 && asym_id.size() && asym_id != line.substr(21, 1))
                {
                    after_ter = true;
                    continue;
                }
                asym_id = line[21];
            }
            buf_all_atm << "ATOM  " << setw(5) << lig_idx1
                        << line.substr(11, 9) << " A" << line.substr(22, 8)
                        << setiosflags(ios::fixed) << setprecision(3)
                        << setw(8) << x1[0] << setw(8) << x1[1] << setw(8) << x1[2] << '\n';
            if (!mm_opt && find(resi_aln1.begin(), resi_aln1.end(),
                                line.substr(22, 4)) != resi_aln1.end())
            {
                buf_atm << "ATOM  " << setw(5) << lig_idx1
                        << line.substr(11, 9) << " A" << line.substr(22, 8)
                        << setiosflags(ios::fixed) << setprecision(3)
                        << setw(8) << x1[0] << setw(8) << x1[1] << setw(8) << x1[2] << '\n';
            }
            if (line.substr(12, 4) != " CA " && line.substr(12, 4) != " C3'")
                continue;
            ca_idx1++;
            buf_all << "ATOM  " << setw(5) << ca_idx1 << ' '
                    << line.substr(12, 4) << ' ' << line.substr(17, 3) << " A" << line.substr(22, 8)
                    << setiosflags(ios::fixed) << setprecision(3)
                    << setw(8) << x1[0] << setw(8) << x1[1] << setw(8) << x1[2] << '\n';
            if (find(resi_aln1.begin(), resi_aln1.end(),
                     line.substr(22, 4)) == resi_aln1.end())
                continue;
            if (!mm_opt)
                buf << "ATOM  " << setw(5) << ca_idx1 << ' '
                    << line.substr(12, 4) << ' ' << line.substr(17, 3) << " A" << line.substr(22, 8)
                    << setiosflags(ios::fixed) << setprecision(3)
                    << setw(8) << x1[0] << setw(8) << x1[1] << setw(8) << x1[2] << '\n';
            idx_vec.push_back(ca_idx1);
        }
        else if (line.compare(0, 5, "loop_") == 0) // PDBx/mmCIF
        {
            while (1)
            {
                if (fin.good())
                    getline(fin, line);
                else
                    PrintErrorAndQuit("ERROR! Unexpected end of " + xname);
                if (line.size())
                    break;
            }
            if (line.compare(0, 11, "_atom_site."))
                continue;
            _atom_site.clear();
            atom_site_pos = 0;
            _atom_site[line.substr(11, line.size() - 12)] = atom_site_pos;
            while (1)
            {
                if (fin.good())
                    getline(fin, line);
                else
                    PrintErrorAndQuit("ERROR! Unexpected end of " + xname);
                if (line.size() == 0)
                    continue;
                if (line.compare(0, 11, "_atom_site."))
                    break;
                _atom_site[line.substr(11, line.size() - 12)] = ++atom_site_pos;
            }

            if (is_mmcif == false)
            {
                // buf_pdb.str(string());
                is_mmcif = true;
            }

            while (1)
            {
                line_vec.clear();
                split(line, line_vec);
                if (line_vec[_atom_site["group_PDB"]] != "ATOM" &&
                    line_vec[_atom_site["group_PDB"]] != "HETATM")
                    break;
                if (_atom_site.count("pdbx_PDB_model_num"))
                {
                    if (model_index.size() && model_index !=
                                                  line_vec[_atom_site["pdbx_PDB_model_num"]])
                        break;
                    model_index = line_vec[_atom_site["pdbx_PDB_model_num"]];
                }

                x[0] = atof(line_vec[_atom_site["Cartn_x"]].c_str());
                x[1] = atof(line_vec[_atom_site["Cartn_y"]].c_str());
                x[2] = atof(line_vec[_atom_site["Cartn_z"]].c_str());
                if (mirror_opt)
                    x[2] = -x[2];

                if (_atom_site.count("auth_seq_id"))
                    resi = line_vec[_atom_site["auth_seq_id"]];
                else
                    resi = line_vec[_atom_site["label_seq_id"]];
                if (_atom_site.count("pdbx_PDB_ins_code") &&
                    line_vec[_atom_site["pdbx_PDB_ins_code"]] != "?")
                    resi += line_vec[_atom_site["pdbx_PDB_ins_code"]][0];
                else
                    resi += " ";
                if (read_resi >= 2)
                {
                    if (_atom_site.count("auth_asym_id"))
                        asym_id = line_vec[_atom_site["auth_asym_id"]];
                    else
                        asym_id = line_vec[_atom_site["label_asym_id"]];
                    if (asym_id == ".")
                        asym_id = " ";
                    resi += asym_id[0];
                }
                hinge = 0;
                if (resi2hinge_dict.count(resi))
                    hinge = resi2hinge_dict[resi];
                tu2t_u(tu_vec[hinge], t, u);
                transform(t, u, x, x1);

                if (_atom_site.count("label_alt_id") == 0 ||
                    line_vec[_atom_site["label_alt_id"]] == "." ||
                    line_vec[_atom_site["label_alt_id"]] == "A")
                {
                    atom = line_vec[_atom_site["label_atom_id"]];
                    if (atom[0] == '"')
                        atom = atom.substr(1);
                    if (atom.size() && atom[atom.size() - 1] == '"')
                        atom = atom.substr(0, atom.size() - 1);
                    if (atom.size() == 0)
                        atom = "    ";
                    else if (atom.size() == 1)
                        atom = " " + atom + "  ";
                    else if (atom.size() == 2)
                        atom = " " + atom + " ";
                    else if (atom.size() == 3)
                        atom = " " + atom;
                    else if (atom.size() >= 5)
                        atom = atom.substr(0, 4);

                    AA = line_vec[_atom_site["label_comp_id"]]; // residue name
                    if (AA.size() == 1)
                        AA = "  " + AA;
                    else if (AA.size() == 2)
                        AA = " " + AA;
                    else if (AA.size() >= 4)
                        AA = AA.substr(0, 3);

                    if (_atom_site.count("auth_seq_id"))
                        resi = line_vec[_atom_site["auth_seq_id"]];
                    else
                        resi = line_vec[_atom_site["label_seq_id"]];
                    while (resi.size() < 4)
                        resi = ' ' + resi;
                    if (resi.size() > 4)
                        resi = resi.substr(0, 4);

                    inscode = ' ';
                    if (_atom_site.count("pdbx_PDB_ins_code") &&
                        line_vec[_atom_site["pdbx_PDB_ins_code"]] != "?")
                        inscode = line_vec[_atom_site["pdbx_PDB_ins_code"]][0];

                    if (_atom_site.count("auth_asym_id"))
                    {
                        if (chain1_sele.size())
                            after_ter = line_vec[_atom_site["auth_asym_id"]] != chain1_sele;
                        else if (ter_opt >= 2 && ca_idx1 && asym_id.size() &&
                                 asym_id != line_vec[_atom_site["auth_asym_id"]])
                            after_ter = true;
                        asym_id = line_vec[_atom_site["auth_asym_id"]];
                    }
                    else if (_atom_site.count("label_asym_id"))
                    {
                        if (chain1_sele.size())
                            after_ter = line_vec[_atom_site["label_asym_id"]] != chain1_sele;
                        if (ter_opt >= 2 && ca_idx1 && asym_id.size() &&
                            asym_id != line_vec[_atom_site["label_asym_id"]])
                            after_ter = true;
                        asym_id = line_vec[_atom_site["label_asym_id"]];
                    }
                    // buf_pdb<<left<<setw(6)
                    //<<line_vec[_atom_site["group_PDB"]]<<right
                    //<<setw(5)<<lig_idx1%100000<<' '<<atom<<' '
                    //<<AA<<" "<<asym_id[asym_id.size()-1]
                    //<<resi<<inscode<<"   "
                    //<<setiosflags(ios::fixed)<<setprecision(3)
                    //<<setw(8)<<x1[0]
                    //<<setw(8)<<x1[1]
                    //<<setw(8)<<x1[2]<<'\n';

                    if (after_ter == false ||
                        line_vec[_atom_site["group_pdb"]] == "HETATM")
                    {
                        lig_idx1++;
                        buf_all_atm_lig << left << setw(6)
                                        << line_vec[_atom_site["group_PDB"]] << right
                                        << setw(5) << lig_idx1 % 100000 << ' ' << atom << ' '
                                        << AA << " A" << resi << inscode << "   "
                                        << setiosflags(ios::fixed) << setprecision(3)
                                        << setw(8) << x1[0]
                                        << setw(8) << x1[1]
                                        << setw(8) << x1[2] << '\n';
                        if (after_ter == false &&
                            line_vec[_atom_site["group_PDB"]] == "ATOM")
                        {
                            buf_all_atm << "ATOM  " << setw(6)
                                        << setw(5) << lig_idx1 % 100000 << ' ' << atom << ' '
                                        << AA << " A" << resi << inscode << "   "
                                        << setiosflags(ios::fixed) << setprecision(3)
                                        << setw(8) << x1[0]
                                        << setw(8) << x1[1]
                                        << setw(8) << x1[2] << '\n';
                            if (!mm_opt && find(resi_aln1.begin(),
                                                resi_aln1.end(), resi) != resi_aln1.end())
                            {
                                buf_atm << "ATOM  " << setw(6)
                                        << setw(5) << lig_idx1 % 100000 << ' '
                                        << atom << ' ' << AA << " A" << resi << inscode << "   "
                                        << setiosflags(ios::fixed) << setprecision(3)
                                        << setw(8) << x1[0]
                                        << setw(8) << x1[1]
                                        << setw(8) << x1[2] << '\n';
                            }
                            if (atom == " CA " || atom == " C3'")
                            {
                                ca_idx1++;
                                // mm_opt, split_opt, mirror_opt, chainID1,chainID2);
                                buf_all << "ATOM  " << setw(6)
                                        << setw(5) << ca_idx1 % 100000 << ' ' << atom << ' '
                                        << AA << " A" << resi << inscode << "   "
                                        << setiosflags(ios::fixed) << setprecision(3)
                                        << setw(8) << x1[0]
                                        << setw(8) << x1[1]
                                        << setw(8) << x1[2] << '\n';
                                if (!mm_opt && find(resi_aln1.begin(),
                                                    resi_aln1.end(), resi) != resi_aln1.end())
                                {
                                    buf << "ATOM  " << setw(6)
                                        << setw(5) << ca_idx1 % 100000 << ' ' << atom << ' '
                                        << AA << " A" << resi << inscode << "   "
                                        << setiosflags(ios::fixed) << setprecision(3)
                                        << setw(8) << x1[0]
                                        << setw(8) << x1[1]
                                        << setw(8) << x1[2] << '\n';
                                    idx_vec.push_back(ca_idx1);
                                }
                            }
                        }
                    }
                }

                while (1)
                {
                    if (fin.good())
                        getline(fin, line);
                    else
                        break;
                    if (line.size())
                        break;
                }
            }
        }
        else if (line.size() && is_mmcif == false)
        {
            // buf_pdb<<line<<'\n';
            if (ter_opt >= 1 && line.compare(0, 3, "END") == 0)
                break;
        }
    }
    fin.close();
    if (!mm_opt)
        buf << "TER\n";
    buf_all << "TER\n";
    if (!mm_opt)
        buf_atm << "TER\n";
    buf_all_atm << "TER\n";
    buf_all_atm_lig << "TER\n";
    for (i = 1; i < ca_idx1; i++)
        buf_all << "CONECT"
                << setw(5) << i % 100000 << setw(5) << (i + 1) % 100000 << '\n';
    if (!mm_opt)
        for (i = 1; i < idx_vec.size(); i++)
            buf << "CONECT"
                << setw(5) << idx_vec[i - 1] % 100000 << setw(5) << idx_vec[i] % 100000 << '\n';
    idx_vec.clear();

    /* read second file */
    after_ter = false;
    asym_id = "";
    fin.open(yname.c_str());
    while (fin.good())
    {
        getline(fin, line);
        if (ter_opt >= 3 && line.compare(0, 3, "TER") == 0)
            after_ter = true;
        if (line.size() >= 54 && (line.compare(0, 6, "ATOM  ") == 0 ||
                                  line.compare(0, 6, "HETATM") == 0)) // PDB format
        {
            if (line[16] != 'A' && line[16] != ' ')
                continue;
            if (after_ter && line.compare(0, 6, "ATOM  ") == 0)
                continue;
            lig_idx2++;
            buf_all_atm_lig << line.substr(0, 6) << setw(5) << lig_idx1 + lig_idx2
                            << line.substr(11, 9) << " B" << line.substr(22, 32) << '\n';
            if (chain1_sele.size() && line[21] != chain1_sele[0])
                continue;
            if (after_ter || line.compare(0, 6, "ATOM  "))
                continue;
            if (ter_opt >= 2)
            {
                if (ca_idx2 && asym_id.size() && asym_id != line.substr(21, 1))
                {
                    after_ter = true;
                    continue;
                }
                asym_id = line[21];
            }
            buf_all_atm << "ATOM  " << setw(5) << lig_idx1 + lig_idx2
                        << line.substr(11, 9) << " B" << line.substr(22, 32) << '\n';
            if (!mm_opt && find(resi_aln2.begin(), resi_aln2.end(),
                                line.substr(22, 4)) != resi_aln2.end())
            {
                buf_atm << "ATOM  " << setw(5) << lig_idx1 + lig_idx2
                        << line.substr(11, 9) << " B" << line.substr(22, 32) << '\n';
            }
            if (line.substr(12, 4) != " CA " && line.substr(12, 4) != " C3'")
                continue;
            ca_idx2++;
            buf_all << "ATOM  " << setw(5) << ca_idx1 + ca_idx2 << ' ' << line.substr(12, 4)
                    << ' ' << line.substr(17, 3) << " B" << line.substr(22, 32) << '\n';
            if (find(resi_aln2.begin(), resi_aln2.end(), line.substr(22, 4)) == resi_aln2.end())
                continue;
            if (!mm_opt)
                buf << "ATOM  " << setw(5) << ca_idx1 + ca_idx2 << ' '
                    << line.substr(12, 4) << ' ' << line.substr(17, 3) << " B"
                    << line.substr(22, 32) << '\n';
            idx_vec.push_back(ca_idx1 + ca_idx2);
        }
        else if (line.compare(0, 5, "loop_") == 0) // PDBx/mmCIF
        {
            while (1)
            {
                if (fin.good())
                    getline(fin, line);
                else
                    PrintErrorAndQuit("ERROR! Unexpected end of " + yname);
                if (line.size())
                    break;
            }
            if (line.compare(0, 11, "_atom_site."))
                continue;
            _atom_site.clear();
            atom_site_pos = 0;
            _atom_site[line.substr(11, line.size() - 12)] = atom_site_pos;
            while (1)
            {
                if (fin.good())
                    getline(fin, line);
                else
                    PrintErrorAndQuit("ERROR! Unexpected end of " + yname);
                if (line.size() == 0)
                    continue;
                if (line.compare(0, 11, "_atom_site."))
                    break;
                _atom_site[line.substr(11, line.size() - 12)] = ++atom_site_pos;
            }

            while (1)
            {
                line_vec.clear();
                split(line, line_vec);
                if (line_vec[_atom_site["group_PDB"]] != "ATOM" &&
                    line_vec[_atom_site["group_PDB"]] != "HETATM")
                    break;
                if (_atom_site.count("pdbx_PDB_model_num"))
                {
                    if (model_index.size() && model_index !=
                                                  line_vec[_atom_site["pdbx_PDB_model_num"]])
                        break;
                    model_index = line_vec[_atom_site["pdbx_PDB_model_num"]];
                }

                if (_atom_site.count("label_alt_id") == 0 ||
                    line_vec[_atom_site["label_alt_id"]] == "." ||
                    line_vec[_atom_site["label_alt_id"]] == "A")
                {
                    atom = line_vec[_atom_site["label_atom_id"]];
                    if (atom[0] == '"')
                        atom = atom.substr(1);
                    if (atom.size() && atom[atom.size() - 1] == '"')
                        atom = atom.substr(0, atom.size() - 1);
                    if (atom.size() == 0)
                        atom = "    ";
                    else if (atom.size() == 1)
                        atom = " " + atom + "  ";
                    else if (atom.size() == 2)
                        atom = " " + atom + " ";
                    else if (atom.size() == 3)
                        atom = " " + atom;
                    else if (atom.size() >= 5)
                        atom = atom.substr(0, 4);

                    AA = line_vec[_atom_site["label_comp_id"]]; // residue name
                    if (AA.size() == 1)
                        AA = "  " + AA;
                    else if (AA.size() == 2)
                        AA = " " + AA;
                    else if (AA.size() >= 4)
                        AA = AA.substr(0, 3);

                    if (_atom_site.count("auth_seq_id"))
                        resi = line_vec[_atom_site["auth_seq_id"]];
                    else
                        resi = line_vec[_atom_site["label_seq_id"]];
                    while (resi.size() < 4)
                        resi = ' ' + resi;
                    if (resi.size() > 4)
                        resi = resi.substr(0, 4);

                    inscode = ' ';
                    if (_atom_site.count("pdbx_PDB_ins_code") &&
                        line_vec[_atom_site["pdbx_PDB_ins_code"]] != "?")
                        inscode = line_vec[_atom_site["pdbx_PDB_ins_code"]][0];

                    if (_atom_site.count("auth_asym_id"))
                    {
                        if (chain2_sele.size())
                            after_ter = line_vec[_atom_site["auth_asym_id"]] != chain2_sele;
                        if (ter_opt >= 2 && ca_idx2 && asym_id.size() &&
                            asym_id != line_vec[_atom_site["auth_asym_id"]])
                            after_ter = true;
                        asym_id = line_vec[_atom_site["auth_asym_id"]];
                    }
                    else if (_atom_site.count("label_asym_id"))
                    {
                        if (chain2_sele.size())
                            after_ter = line_vec[_atom_site["label_asym_id"]] != chain2_sele;
                        if (ter_opt >= 2 && ca_idx2 && asym_id.size() &&
                            asym_id != line_vec[_atom_site["label_asym_id"]])
                            after_ter = true;
                        asym_id = line_vec[_atom_site["label_asym_id"]];
                    }
                    if (after_ter == false ||
                        line_vec[_atom_site["group_PDB"]] == "HETATM")
                    {
                        lig_idx2++;
                        buf_all_atm_lig << left << setw(6)
                                        << line_vec[_atom_site["group_PDB"]] << right
                                        << setw(5) << (lig_idx1 + lig_idx2) % 100000 << ' '
                                        << atom << ' ' << AA << " B" << resi << inscode << "   "
                                        << setw(8) << line_vec[_atom_site["Cartn_x"]]
                                        << setw(8) << line_vec[_atom_site["Cartn_y"]]
                                        << setw(8) << line_vec[_atom_site["Cartn_z"]]
                                        << '\n';
                        if (after_ter == false &&
                            line_vec[_atom_site["group_PDB"]] == "ATOM")
                        {
                            buf_all_atm << "ATOM  " << setw(6)
                                        << setw(5) << (lig_idx1 + lig_idx2) % 100000 << ' '
                                        << atom << ' ' << AA << " B" << resi << inscode << "   "
                                        << setw(8) << line_vec[_atom_site["Cartn_x"]]
                                        << setw(8) << line_vec[_atom_site["Cartn_y"]]
                                        << setw(8) << line_vec[_atom_site["Cartn_z"]]
                                        << '\n';
                            if (!mm_opt && find(resi_aln2.begin(),
                                                resi_aln2.end(), resi) != resi_aln2.end())
                            {
                                buf_atm << "ATOM  " << setw(6)
                                        << setw(5) << (lig_idx1 + lig_idx2) % 100000 << ' '
                                        << atom << ' ' << AA << " B" << resi << inscode << "   "
                                        << setw(8) << line_vec[_atom_site["Cartn_x"]]
                                        << setw(8) << line_vec[_atom_site["Cartn_y"]]
                                        << setw(8) << line_vec[_atom_site["Cartn_z"]]
                                        << '\n';
                            }
                            if (atom == " CA " || atom == " C3'")
                            {
                                ca_idx2++;
                                buf_all << "ATOM  " << setw(6)
                                        << setw(5) << (ca_idx1 + ca_idx2) % 100000
                                        << ' ' << atom << ' ' << AA << " B" << resi << inscode << "   "
                                        << setw(8) << line_vec[_atom_site["Cartn_x"]]
                                        << setw(8) << line_vec[_atom_site["Cartn_y"]]
                                        << setw(8) << line_vec[_atom_site["Cartn_z"]]
                                        << '\n';
                                if (!mm_opt && find(resi_aln2.begin(),
                                                    resi_aln2.end(), resi) != resi_aln2.end())
                                {
                                    buf << "ATOM  " << setw(6)
                                        << setw(5) << (ca_idx1 + ca_idx2) % 100000
                                        << ' ' << atom << ' ' << AA << " B" << resi << inscode << "   "
                                        << setw(8) << line_vec[_atom_site["Cartn_x"]]
                                        << setw(8) << line_vec[_atom_site["Cartn_y"]]
                                        << setw(8) << line_vec[_atom_site["Cartn_z"]]
                                        << '\n';
                                    idx_vec.push_back(ca_idx1 + ca_idx2);
                                }
                            }
                        }
                    }
                }

                if (fin.good())
                    getline(fin, line);
                else
                    break;
            }
        }
        else if (line.size())
        {
            if (ter_opt >= 1 && line.compare(0, 3, "END") == 0)
                break;
        }
    }
    fin.close();
    if (!mm_opt)
        buf << "TER\n";
    buf_all << "TER\n";
    if (!mm_opt)
        buf_atm << "TER\n";
    buf_all_atm << "TER\n";
    buf_all_atm_lig << "TER\n";
    for (i = ca_idx1 + 1; i < ca_idx1 + ca_idx2; i++)
        buf_all << "CONECT"
                << setw(5) << i % 100000 << setw(5) << (i + 1) % 100000 << '\n';
    for (i = 1; i < idx_vec.size(); i++)
        buf << "CONECT"
            << setw(5) << idx_vec[i - 1] % 100000 << setw(5) << idx_vec[i] % 100000 << '\n';
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
        fp << buf.str();
        fp.close();
    }
    fp.open((fname_super + "_all").c_str());
    fp << buf_all.str();
    fp.close();
    if (!mm_opt)
    {
        fp.open((fname_super + "_atm").c_str());
        fp << buf_atm.str();
        fp.close();
    }
    fp.open((fname_super + "_all_atm").c_str());
    fp << buf_all_atm.str();
    fp.close();
    fp.open((fname_super + "_all_atm_lig").c_str());
    fp << buf_all_atm_lig.str();
    fp.close();
    // fp.open((fname_super+".pdb").c_str());
    // fp<<buf_pdb.str();
    // fp.close();

    /* clear stream */
    buf.str(string());
    buf_all.str(string());
    buf_atm.str(string());
    buf_all_atm.str(string());
    buf_all_atm_lig.str(string());
    // buf_pdb.str(string());
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
                            const string fname_super, const vector<vector<double>> &tu_vec,
                            double t[3], double u[3][3], const int ter_opt,
                            const int mm_opt, const int split_opt, const int mirror_opt,
                            const char *seqM, const char *seqxA, const char *seqyA,
                            const vector<string> &resi_vec1, const vector<string> &resi_vec2,
                            const string chainID1, const string chainID2)
{
    int compress_type = 0; // uncompressed file
    ifstream fin;
#ifndef REDI_PSTREAM_H_SEEN
    ifstream fin_gz;
#else
    redi::ipstream fin_gz; // if file is compressed
    if (xname.size() >= 3 &&
        xname.substr(xname.size() - 3, 3) == ".gz")
    {
        fin_gz.open("gunzip -c " + xname);
        compress_type = 1;
    }
    else if (xname.size() >= 4 &&
             xname.substr(xname.size() - 4, 4) == ".bz2")
    {
        fin_gz.open("bzcat " + xname);
        compress_type = 2;
    }
    else
#endif
    fin.open(xname.c_str());

    map<string, int> resi2hinge_dict;
    int r, i, j;
    j = -1;
    char hinge_char = 0;
    int xlen = resi_vec1.size();
    int ali_len = strlen(seqM);
    for (r = 0; r < strlen(seqxA); r++)
    {
        if (seqxA[r] == '-')
            continue;
        j++;
        hinge_char = seqM[r];
        if (hinge_char == ' ')
        {
            for (i = 1; i < ali_len; i++)
            {
                if (r - i >= 0 && seqM[r - i] != ' ')
                    hinge_char = seqM[r - i];
                else if (r + i < xlen && seqM[r + i] != ' ')
                    hinge_char = seqM[r + i];
                if (hinge_char != ' ')
                    break;
            }
        }
        int hinge_idx = 0;
        if (hinge_char >= '0' && hinge_char <= '9')
        {
            hinge_idx = hinge_char - '0';
        }
        else if (hinge_char >= 'a' && hinge_char <= 'z')
        {
            hinge_idx = hinge_char - 'a' + 10;
        }
        else if (hinge_char >= 'A' && hinge_char <= 'Z')
        {
            hinge_idx = hinge_char - 'A' + 36;
        }
        resi2hinge_dict[resi_vec1[j]] = hinge_idx;
    }
    string resi = resi_vec1[0];
    int read_resi = resi.size() - 4;

    stringstream buf;
    stringstream buf_pymol;
    string line;
    double x[3];  // before transform
    double x1[3]; // after transform

    /* for PDBx/mmCIF only */
    map<string, int> _atom_site;
    size_t atom_site_pos;
    vector<string> line_vec;
    int infmt = -1; // 0 - PDB, 3 - PDBx/mmCIF
    int hinge = 0;
    string asym_id = "."; // this is similar to chainID, except that
                          // chainID is char while asym_id is a string
                          // with possibly multiple char
    while (compress_type ? fin_gz.good() : fin.good())
    {
        if (compress_type)
            getline(fin_gz, line);
        else
            getline(fin, line);
        if (line.compare(0, 6, "ATOM  ") == 0 ||
            line.compare(0, 6, "HETATM") == 0) // PDB format
        {
            infmt = 0;
            x[0] = atof(line.substr(30, 8).c_str());
            x[1] = atof(line.substr(38, 8).c_str());
            x[2] = atof(line.substr(46, 8).c_str());
            if (mirror_opt)
                x[2] = -x[2];
            if (read_resi == 1)
                resi = line.substr(22, 5);
            else
                resi = line.substr(22, 5) + line[21];
            hinge = 0;
            if (resi2hinge_dict.count(resi))
                hinge = resi2hinge_dict[resi];
            tu2t_u(tu_vec[hinge], t, u);
            transform(t, u, x, x1);
            buf << line.substr(0, 30) << setiosflags(ios::fixed)
                << setprecision(3)
                << setw(8) << x1[0] << setw(8) << x1[1] << setw(8) << x1[2]
                << line.substr(54) << '\n';
        }
        else if (line.compare(0, 5, "loop_") == 0) // PDBx/mmCIF
        {
            infmt = 3;
            buf << line << '\n';
            while (1)
            {
                if (compress_type)
                {
                    if (fin_gz.good())
                        getline(fin_gz, line);
                    else
                        PrintErrorAndQuit("ERROR! Unexpected end of " + xname);
                }
                else
                {
                    if (fin.good())
                        getline(fin, line);
                    else
                        PrintErrorAndQuit("ERROR! Unexpected end of " + xname);
                }
                if (line.size())
                    break;
            }
            buf << line << '\n';
            if (line.compare(0, 11, "_atom_site."))
                continue;
            _atom_site.clear();
            atom_site_pos = 0;
            _atom_site[Trim(line.substr(11))] = atom_site_pos;
            while (1)
            {
                while (1)
                {
                    if (compress_type)
                    {
                        if (fin_gz.good())
                            getline(fin_gz, line);
                        else
                            PrintErrorAndQuit("ERROR! Unexpected end of " + xname);
                    }
                    else
                    {
                        if (fin.good())
                            getline(fin, line);
                        else
                            PrintErrorAndQuit("ERROR! Unexpected end of " + xname);
                    }
                    if (line.size())
                        break;
                }
                if (line.compare(0, 11, "_atom_site."))
                    break;
                _atom_site[Trim(line.substr(11))] = ++atom_site_pos;
                buf << line << '\n';
            }

            if (_atom_site.count("group_PDB") *
                    _atom_site.count("Cartn_x") *
                    _atom_site.count("Cartn_y") *
                    _atom_site.count("Cartn_z") ==
                0)
            {
                buf << line << '\n';
                cerr << "Warning! Missing one of the following _atom_site data items: group_PDB, Cartn_x, Cartn_y, Cartn_z" << endl;
                continue;
            }

            while (1)
            {
                line_vec.clear();
                split(line, line_vec);
                if (line_vec[_atom_site["group_PDB"]] != "ATOM" &&
                    line_vec[_atom_site["group_PDB"]] != "HETATM")
                    break;

                x[0] = atof(line_vec[_atom_site["Cartn_x"]].c_str());
                x[1] = atof(line_vec[_atom_site["Cartn_y"]].c_str());
                x[2] = atof(line_vec[_atom_site["Cartn_z"]].c_str());
                if (mirror_opt)
                    x[2] = -x[2];

                if (_atom_site.count("auth_seq_id"))
                    resi = line_vec[_atom_site["auth_seq_id"]];
                else
                    resi = line_vec[_atom_site["label_seq_id"]];
                if (_atom_site.count("pdbx_PDB_ins_code") &&
                    line_vec[_atom_site["pdbx_PDB_ins_code"]] != "?")
                    resi += line_vec[_atom_site["pdbx_PDB_ins_code"]][0];
                else
                    resi += " ";
                if (read_resi >= 2)
                {
                    if (_atom_site.count("auth_asym_id"))
                        asym_id = line_vec[_atom_site["auth_asym_id"]];
                    else
                        asym_id = line_vec[_atom_site["label_asym_id"]];
                    if (asym_id == ".")
                        asym_id = " ";
                    resi += asym_id[0];
                }
                hinge = 0;
                if (resi2hinge_dict.count(resi))
                    hinge = resi2hinge_dict[resi];
                tu2t_u(tu_vec[hinge], t, u);
                transform(t, u, x, x1);

                for (atom_site_pos = 0; atom_site_pos < _atom_site.size(); atom_site_pos++)
                {
                    if (atom_site_pos == _atom_site["Cartn_x"])
                        buf << setiosflags(ios::fixed) << setprecision(3)
                            << setw(8) << x1[0] << ' ';
                    else if (atom_site_pos == _atom_site["Cartn_y"])
                        buf << setiosflags(ios::fixed) << setprecision(3)
                            << setw(8) << x1[1] << ' ';
                    else if (atom_site_pos == _atom_site["Cartn_z"])
                        buf << setiosflags(ios::fixed) << setprecision(3)
                            << setw(8) << x1[2] << ' ';
                    else
                        buf << line_vec[atom_site_pos] << ' ';
                }
                buf << '\n';

                if (compress_type && fin_gz.good())
                    getline(fin_gz, line);
                else if (!compress_type && fin.good())
                    getline(fin, line);
                else
                    break;
            }
            if (compress_type ? fin_gz.good() : fin.good())
                buf << line << '\n';
        }
        else if (line.size())
        {
            buf << line << '\n';
            if (ter_opt >= 1 && line.compare(0, 3, "END") == 0)
                break;
        }
    }
    if (compress_type)
        fin_gz.close();
    else
        fin.close();

    string fname_super_full = fname_super;
    if (infmt == 0)
        fname_super_full += ".pdb";
    else if (infmt == 3)
        fname_super_full += ".cif";
    ofstream fp;
    fp.open(fname_super_full.c_str());
    fp << buf.str();
    fp.close();
    buf.str(string()); // clear stream

    string chain1_sele;
    string chain2_sele;
    if (!mm_opt)
    {
        if (split_opt == 2 && ter_opt >= 1) // align one chain from model 1
        {
            chain1_sele = " and c. " + chainID1.substr(1);
            chain2_sele = " and c. " + chainID2.substr(1);
        }
        else if (split_opt == 2 && ter_opt == 0) // align one chain from each model
        {
            for (i = 1; i < chainID1.size(); i++)
                if (chainID1[i] == ',')
                    break;
            chain1_sele = " and c. " + chainID1.substr(i + 1);
            for (i = 1; i < chainID2.size(); i++)
                if (chainID2[i] == ',')
                    break;
            chain2_sele = " and c. " + chainID2.substr(i + 1);
        }
    }

    /* extract aligned region */
    int i1 = -1;
    int i2 = -1;
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
        for (i = 0; i < strlen(seqM); i++)
        {
            i1 += (seqxA[i] != '-' && seqxA[i] != '*');
            i2 += (seqyA[i] != '-');
            if (seqM[i] == ' ' || seqxA[i] == '*')
                continue;
            curr_resi1 = resi_vec1[i1].substr(0, 4);
            curr_resi2 = resi_vec2[i2].substr(0, 4);
            if (resi1_sele.size() == 0)
                resi1_sele = "i. " + curr_resi1;
            else
            {
                resi1_sele += " or i. " + curr_resi1;
                resi1_bond += "bond structure1 and i. " + prev_resi1 +
                              ", i. " + curr_resi1 + "\n";
            }
            if (resi2_sele.size() == 0)
                resi2_sele = "i. " + curr_resi2;
            else
            {
                resi2_sele += " or i. " + curr_resi2;
                resi2_bond += "bond structure2 and i. " + prev_resi2 +
                              ", i. " + curr_resi2 + "\n";
            }
            prev_resi1 = curr_resi1;
            prev_resi2 = curr_resi2;
            // if (seqM[i]!=':') continue;
        }
        if (resi1_sele.size())
            resi1_sele = " and ( " + resi1_sele + ")";
        if (resi2_sele.size())
            resi2_sele = " and ( " + resi2_sele + ")";
    }

    /* write pymol script */
    vector<string> pml_list;
    pml_list.push_back(fname_super + "");
    pml_list.push_back(fname_super + "_atm");
    pml_list.push_back(fname_super + "_all");
    pml_list.push_back(fname_super + "_all_atm");
    pml_list.push_back(fname_super + "_all_atm_lig");

    for (int p = 0; p < pml_list.size(); p++)
    {
        if (mm_opt && p <= 1)
            continue;
        buf_pymol
            << "#!/usr/bin/env pymol\n"
            << "cmd.load(\"" << fname_super_full << "\", \"structure1\")\n"
            << "cmd.load(\"" << yname << "\", \"structure2\")\n"
            << "hide all\n"
            << "set all_states, " << ((ter_opt == 0) ? "on" : "off") << '\n';
        if (p == 0) // .pml
        {
            if (chain1_sele.size())
                buf_pymol
                    << "remove structure1 and not " << chain1_sele.substr(4) << "\n";
            if (chain2_sele.size())
                buf_pymol
                    << "remove structure2 and not " << chain2_sele.substr(4) << "\n";
            buf_pymol
                << "remove not n. CA and not n. C3'\n"
                << resi1_bond
                << resi2_bond
                << "show stick, structure1" << chain1_sele << resi1_sele << "\n"
                << "show stick, structure2" << chain2_sele << resi2_sele << "\n";
        }
        else if (p == 1) // _atm.pml
        {
            buf_pymol
                << "show cartoon, structure1" << chain1_sele << resi1_sele << "\n"
                << "show cartoon, structure2" << chain2_sele << resi2_sele << "\n";
        }
        else if (p == 2) // _all.pml
        {
            buf_pymol
                << "show ribbon, structure1" << chain1_sele << "\n"
                << "show ribbon, structure2" << chain2_sele << "\n";
        }
        else if (p == 3) // _all_atm.pml
        {
            buf_pymol
                << "show cartoon, structure1" << chain1_sele << "\n"
                << "show cartoon, structure2" << chain2_sele << "\n";
        }
        else if (p == 4) // _all_atm_lig.pml
        {
            buf_pymol
                << "show cartoon, structure1\n"
                << "show cartoon, structure2\n"
                << "show stick, not polymer\n"
                << "show sphere, not polymer\n";
        }
        buf_pymol
            << "color blue, structure1\n"
            << "color red, structure2\n"
            << "set ribbon_width, 6\n"
            << "set stick_radius, 0.3\n"
            << "set sphere_scale, 0.25\n"
            << "set ray_shadow, 0\n"
            << "bg_color white\n"
            << "set transparency=0.2\n"
            << "zoom polymer and ((structure1" << chain1_sele
            << ") or (structure2" << chain2_sele << "))\n"
            << endl;

        fp.open((pml_list[p] + ".pml").c_str());
        fp << buf_pymol.str();
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

// output the final results
void output_flexalign_results(const string xname, const string yname,
                              const string chainID1, const string chainID2,
                              const int xlen, const int ylen, double t[3], double u[3][3],
                              const vector<vector<double>> &tu_vec, const double TM1, const double TM2,
                              const double TM3, const double TM4, const double TM5,
                              const double rmsd, const double d0_out, const char *seqM,
                              const char *seqxA, const char *seqyA, const double Liden,
                              const int n_ali8, const int L_ali, const double TM_ali,
                              const double rmsd_ali, const double TM_0, const double d0_0,
                              const double d0A, const double d0B, const double Lnorm_ass,
                              const double d0_scale, const double d0a, const double d0u,
                              const char *fname_matrix, const int outfmt_opt, const int ter_opt,
                              const int mm_opt, const int split_opt, const int o_opt,
                              const string fname_super, const int i_opt, const int a_opt,
                              const bool u_opt, const bool d_opt, const int mirror_opt,
                              const vector<string> &resi_vec1, const vector<string> &resi_vec2)
{
    if (outfmt_opt <= 0)
    {
        printf("\nName of Structure_1: %s%s (to be superimposed onto Structure_2)\n",
               xname.c_str(), chainID1.c_str());
        printf("Name of Structure_2: %s%s\n", yname.c_str(), chainID2.c_str());
        printf("Length of Structure_1: %d residues\n", xlen);
        printf("Length of Structure_2: %d residues\n\n", ylen);

        if (i_opt)
            printf("User-specified initial alignment: TM/Lali/rmsd = %7.5lf, %4d, %6.3lf\n", TM_ali, L_ali, rmsd_ali);

        printf("Aligned length= %d, RMSD= %6.2f, Seq_ID=n_identical/n_aligned= %4.3f\n", n_ali8, rmsd, (n_ali8 > 0) ? Liden / n_ali8 : 0);
        printf("TM-score= %6.5f (normalized by length of Structure_1: L=%d, d0=%.2f)\n", TM2, xlen, d0B);
        printf("TM-score= %6.5f (normalized by length of Structure_2: L=%d, d0=%.2f)\n", TM1, ylen, d0A);

        if (a_opt == 1)
            printf("TM-score= %6.5f (if normalized by average length of two structures: L=%.1f, d0=%.2f)\n", TM3, (xlen + ylen) * 0.5, d0a);
        if (u_opt)
            printf("TM-score= %6.5f (normalized by user-specified L=%.2f and d0=%.2f)\n", TM4, Lnorm_ass, d0u);
        if (d_opt)
            printf("TM-score= %6.5f (scaled by user-specified d0=%.2f, and L=%d)\n", TM5, d0_scale, ylen);
        printf("(You should use TM-score normalized by length of the reference structure)\n");

        // output alignment
        printf("\n([0-9,a-z,A-Z] denote different aligned fragment pairs separated by different hinges)\n");
        printf("%s\n", seqxA);
        printf("%s\n", seqM);
        printf("%s\n", seqyA);
    }
    else if (outfmt_opt == 1)
    {
        printf(">%s%s\tL=%d\td0=%.2f\tseqID=%.3f\tTM-score=%.5f\n",
               xname.c_str(), chainID1.c_str(), xlen, d0B, Liden / xlen, TM2);
        printf("%s\n", seqxA);
        printf(">%s%s\tL=%d\td0=%.2f\tseqID=%.3f\tTM-score=%.5f\n",
               yname.c_str(), chainID2.c_str(), ylen, d0A, Liden / ylen, TM1);
        printf("%s\n", seqyA);

        printf("# Lali=%d\tRMSD=%.2f\tseqID_ali=%.3f\n",
               n_ali8, rmsd, (n_ali8 > 0) ? Liden / n_ali8 : 0);

        if (i_opt)
            printf("# User-specified initial alignment: TM=%.5lf\tLali=%4d\trmsd=%.3lf\n", TM_ali, L_ali, rmsd_ali);

        if (a_opt)
            printf("# TM-score=%.5f (normalized by average length of two structures: L=%.1f\td0=%.2f)\n", TM3, (xlen + ylen) * 0.5, d0a);

        if (u_opt)
            printf("# TM-score=%.5f (normalized by user-specified L=%.2f\td0=%.2f)\n", TM4, Lnorm_ass, d0u);

        if (d_opt)
            printf("# TM-score=%.5f (scaled by user-specified d0=%.2f\tL=%d)\n", TM5, d0_scale, ylen);

        printf("$$$$\n");
    }
    else if (outfmt_opt == 2)
    {
        printf("%s%s\t%s%s\t%.4f\t%.4f\t%.2f\t%4.3f\t%4.3f\t%4.3f\t%d\t%d\t%d\t%d",
               xname.c_str(), chainID1.c_str(), yname.c_str(), chainID2.c_str(),
               TM2, TM1, rmsd, Liden / xlen, Liden / ylen, (n_ali8 > 0) ? Liden / n_ali8 : 0,
               xlen, ylen, n_ali8, (int)tu_vec.size());
    }
    cout << endl;

    if (strlen(fname_matrix))
        output_flexalign_rotation_matrix(
            fname_matrix, tu_vec, t, u);

    if (o_opt == 1)
        output_flexalign_pymol(xname, yname, fname_super, tu_vec,
                               t, u, ter_opt, mm_opt, split_opt, mirror_opt, seqM, seqxA, seqyA,
                               resi_vec1, resi_vec2, chainID1, chainID2);
    else if (o_opt == 2)
        output_flexalign_rasmol(xname, yname, fname_super, tu_vec,
                                t, u, ter_opt, mm_opt, split_opt, mirror_opt, seqM, seqxA, seqyA,
                                resi_vec1, resi_vec2, chainID1, chainID2,
                                xlen, ylen, d0A, n_ali8, rmsd, TM1, Liden);
}

// Data structure to hold outputs of flexalign_main to avoid parameter clutter
struct FlexAlignResult
{
    double t0[3];
    double u0[3][3];
    vector<vector<double>> tu_vec;
    double TM1, TM2, TM3, TM4, TM5;
    double d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out;
    string seqM, seqxA, seqyA;
    vector<double> do_vec;
    double rmsd0, Liden, TM_ali, rmsd_ali;
    int L_ali, n_ali, n_ali8, hingeNum;

    FlexAlignResult() : TM1(-1.0), TM2(-1.0), TM3(-1.0), TM4(-1.0), TM5(-1.0),
                        d0_0(0.0), TM_0(0.0), d0A(0.0), d0B(0.0), d0u(0.0), d0a(0.0), d0_out(5.0),
                        rmsd0(0.0), Liden(0.0), TM_ali(0.0), rmsd_ali(0.0),
                        L_ali(0), n_ali(0), n_ali8(0), hingeNum(0)
    {
        for (int i = 0; i < 3; i++)
        {
            t0[i] = 0.0;
            for (int j = 0; j < 3; j++)
                u0[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }
};

enum FlexAlignMode
{
    FLEX_STANDARD = 0,
    FLEX_BEST = 1,
    FLEX_USBCAT = 2
};

// Encapsulates the execution of flexalign_main and its fallback refinement logic
void execute_flexalign_with_fallback(
    double **xa, double **ya, char *seqx, char *seqy, char *secx, char *secy,
    int xlen, int ylen, vector<string> &sequence, const double Lnorm_ass, const double d0_scale,
    const int i_opt, const int a_opt, const bool u_opt, const bool d_opt, const bool force_fast_opt,
    const int mol_type, const int hinge_opt, const int ss_opt, FlexAlignResult &res)
{
    res.hingeNum = flexalign_main(
        xa, ya, seqx, seqy, secx, secy,
        res.t0, res.u0, res.tu_vec, res.TM1, res.TM2, res.TM3, res.TM4, res.TM5,
        res.d0_0, res.TM_0, res.d0A, res.d0B, res.d0u, res.d0a, res.d0_out,
        res.seqM, res.seqxA, res.seqyA, res.do_vec,
        res.rmsd0, res.L_ali, res.Liden, res.TM_ali, res.rmsd_ali, res.n_ali, res.n_ali8,
        xlen, ylen, sequence, Lnorm_ass, d0_scale,
        i_opt, a_opt, u_opt, d_opt, force_fast_opt,
        mol_type, hinge_opt, ss_opt);

    // Fallback compensation when too few hinges are found
    if (hinge_opt && res.hingeNum <= 1 && res.n_ali8 < 0.6 * getmin(xlen, ylen))
    {
        FlexAlignResult res_h;
        res_h.tu_vec.push_back(res.tu_vec[0]);
        tu2t_u(res.tu_vec[0], res_h.t0, res_h.u0);

        res_h.hingeNum = flexalign_main(
            xa, ya, seqx, seqy, secx, secy,
            res_h.t0, res_h.u0, res_h.tu_vec,
            res_h.TM1, res_h.TM2, res_h.TM3, res_h.TM4, res_h.TM5,
            res_h.d0_0, res_h.TM_0, res.d0A, res.d0B, res.d0u, res.d0a, res_h.d0_out,
            res_h.seqM, res_h.seqxA, res_h.seqyA, res_h.do_vec,
            res_h.rmsd0, res_h.L_ali, res_h.Liden, res_h.TM_ali, res_h.rmsd_ali,
            res_h.n_ali, res_h.n_ali8,
            xlen, ylen, sequence, Lnorm_ass, d0_scale, i_opt,
            a_opt, u_opt, d_opt, force_fast_opt,
            mol_type, hinge_opt, ss_opt);

        double TM = (res.TM1 > res.TM2) ? res.TM1 : res.TM2;
        double TM_h = (res_h.TM1 > res_h.TM2) ? res_h.TM1 : res_h.TM2;
        if (TM_h > TM)
        {
            res = res_h; // Safely overwrite with the better refined results
        }
    }
}

// ==========================================
// USBCAT Core Algorithm (flexalign_usbcat_main)
// ==========================================
struct USBCAT_AFP
{
    int i, j, len;
    double score;
};

int flexalign_usbcat_main(double **xa, double **ya,
                          const char *seqx, const char *seqy, const char *secx, const char *secy,
                          double t0[3], double u0[3][3], std::vector<std::vector<double>> &tu_vec,
                          double &TM1, double &TM2, double &TM3, double &TM4, double &TM5,
                          double &d0_0, double &TM_0,
                          double &d0A, double &d0B, double &d0u, double &d0a, double &d0_out,
                          std::string &seqM, std::string &seqxA, std::string &seqyA, std::vector<double> &do_vec,
                          double &rmsd0, int &L_ali, double &Liden,
                          double &TM_ali, double &rmsd_ali, int &n_ali, int &n_ali8,
                          const int xlen, const int ylen,
                          const std::vector<std::string> sequence, const double Lnorm_ass,
                          const double d0_scale, const int i_opt, const int a_opt,
                          const bool u_opt, const bool d_opt, const bool fast_opt,
                          const int mol_type, const int hinge_opt, const int ss_opt,
                          int sparse_val = 0, bool hinge_set = false)
{
    // ==========================================
    // TRUE flexalign_greedy BASELINE (Defender)
    // Run full sequence without generate_bounds slicing!
    // ==========================================
    double best_global_max_TM = -1.0;
    std::vector<std::vector<double>> best_tu_vec;
    double best_t0[3], best_u0[3][3];
    double best_TM1 = 0.0, best_TM2 = 0.0, best_TM3 = 0.0, best_TM4 = 0.0, best_TM5 = 0.0;
    double best_rmsd0 = 0.0, best_Liden = 0.0, best_TM_ali = 0.0, best_rmsd_ali = 0.0;
    int best_L_ali = 0, best_n_ali = 0, best_n_ali8 = 0;
    std::string best_seqM = "", best_seqxA = "", best_seqyA = "";
    std::vector<double> best_do_vec;
    double best_d0A = 0.0, best_d0B = 0.0, best_d0a = 0.0, best_d0u = 0.0;

    bool force_fast_opt_global = (std::min(xlen, ylen) > 1500) ? true : fast_opt;
    std::vector<std::string> local_sequence = sequence;

    for (int cur_ss_opt = 0; cur_ss_opt <= 1; cur_ss_opt++)
    {
        FlexAlignResult base_res;
        execute_flexalign_with_fallback(
            xa, ya, (char *)seqx, (char *)seqy, (char *)secx, (char *)secy,
            xlen, ylen, local_sequence, Lnorm_ass, d0_scale,
            i_opt, a_opt, u_opt, d_opt, force_fast_opt_global,
            mol_type, hinge_opt, cur_ss_opt, base_res);

        double cur_max_TM = (base_res.TM1 > base_res.TM2) ? base_res.TM1 : base_res.TM2;
        if (cur_max_TM > best_global_max_TM)
        {
            best_global_max_TM = cur_max_TM;
            for (int a = 0; a < 3; a++)
            {
                best_t0[a] = base_res.t0[a];
                for (int b = 0; b < 3; b++)
                    best_u0[a][b] = base_res.u0[a][b];
            }
            best_tu_vec = base_res.tu_vec;
            best_TM1 = base_res.TM1;
            best_TM2 = base_res.TM2;
            best_TM3 = base_res.TM3;
            best_TM4 = base_res.TM4;
            best_TM5 = base_res.TM5;
            best_rmsd0 = base_res.rmsd0;
            best_Liden = base_res.Liden;
            best_TM_ali = base_res.TM_ali;
            best_rmsd_ali = base_res.rmsd_ali;
            best_L_ali = base_res.L_ali;
            best_n_ali = base_res.n_ali;
            best_n_ali8 = base_res.n_ali8;
            best_seqM = base_res.seqM;
            best_seqxA = base_res.seqxA;
            best_seqyA = base_res.seqyA;
            best_do_vec = base_res.do_vec;
            best_d0A = base_res.d0A;
            best_d0B = base_res.d0B;
            best_d0a = base_res.d0a;
            best_d0u = base_res.d0u;
        }
    }

    if (best_global_max_TM >= 0.85)
    {
        TM1 = best_TM1;
        TM2 = best_TM2;
        TM3 = best_TM3;
        TM4 = best_TM4;
        TM5 = best_TM5;
        rmsd0 = best_rmsd0;
        Liden = best_Liden;
        TM_ali = best_TM_ali;
        rmsd_ali = best_rmsd_ali;
        L_ali = best_L_ali;
        n_ali = best_n_ali;
        n_ali8 = best_n_ali8;
        seqM = best_seqM;
        seqxA = best_seqxA;
        seqyA = best_seqyA;
        do_vec = best_do_vec;
        tu_vec = best_tu_vec;
        d0A = best_d0A;
        d0B = best_d0B;
        d0a = best_d0a;
        d0u = best_d0u;
        for (int a = 0; a < 3; a++)
        {
            t0[a] = best_t0[a];
            for (int b = 0; b < 3; b++)
                u0[a][b] = best_u0[a][b];
        }
        return tu_vec.size();
    }

    // ==========================================
    // Proceed to USBCAT sliced bounds logic...
    // ==========================================
    int fragLen = 8;
    double resScore = 3.0;
    double gap_ext = -0.5;
    double disCut = 5.0;
    double disSmooth = 4.0;
    double twist_pen = -25.0;
    int max_gap = 40;
    double max_penalty = -5.0;
    int misCut = 2 * fragLen;
    int maxGapFrag = fragLen + max_gap;
    double afp_dis_cut = fragLen * fragLen * (disCut * disCut);
    int max_twists = hinge_opt;

    // OPTIMIZATION 1: Precompute local intra-protein distance matrices
    int max_dist_window = max_gap + 2 * fragLen + 1;
    std::vector<std::vector<double>> disTable1(xlen, std::vector<double>(max_dist_window, 0.0));
    std::vector<std::vector<double>> disTable2(ylen, std::vector<double>(max_dist_window, 0.0));

    for (int i = 0; i < xlen; i++)
    {
        for (int j = i; j < std::min(xlen, i + max_dist_window); j++)
            disTable1[i][j - i] = std::sqrt(dist(xa[i], xa[j]));
    }
    for (int i = 0; i < ylen; i++)
    {
        for (int j = i; j < std::min(ylen, i + max_dist_window); j++)
            disTable2[i][j - i] = std::sqrt(dist(ya[i], ya[j]));
    }

    // Wrapper for generating bounds
    auto generate_bounds = [&](double cur_rmsdCut, double cur_badRmsd, double cur_local_badRmsd) -> std::pair<std::vector<int>, std::vector<int>>
    {
        // Step 1: Extract initial AFPs in batches
        std::vector<USBCAT_AFP> initial_afps;
        int step = sparse_val + 1;

        double r1_static[8][3], r2_static[8][3];
        double *r1[8], *r2[8];
        for (int k = 0; k < 8; k++)
        {
            r1[k] = r1_static[k];
            r2[k] = r2_static[k];
        }

        for (int i = 0; i <= xlen - fragLen; i += step)
        {
            for (int j = 0; j <= ylen - fragLen; j += step)
            {
                int d3_term = std::min(i, j) + std::min(xlen - (i + fragLen - 1), ylen - (j + fragLen)) + fragLen;
                double d3_cutoff = std::max(0.3, best_global_max_TM - 0.3);
                if (d3_term < d3_cutoff * std::min(xlen, ylen))
                    continue;

                double dist1 = disTable1[i][fragLen - 1];
                double dist2 = disTable2[j][fragLen - 1];

                if (std::fabs(dist1 - dist2) > 2.0 * cur_rmsdCut)
                    continue;

                for (int k = 0; k < fragLen; k++)
                {
                    r1[k][0] = xa[i + k][0];
                    r1[k][1] = xa[i + k][1];
                    r1[k][2] = xa[i + k][2];
                    r2[k][0] = ya[j + k][0];
                    r2[k][1] = ya[j + k][1];
                    r2[k][2] = ya[j + k][2];
                }

                double rms_sum_sq, t_tmp[3], u_tmp[3][3];
                Kabsch(r1, r2, fragLen, 0, &rms_sum_sq, t_tmp, u_tmp);
                double rmsd_tmp = std::sqrt(rms_sum_sq / fragLen);

                if (rmsd_tmp < cur_rmsdCut)
                {
                    USBCAT_AFP afp;
                    afp.i = i;
                    afp.j = j;
                    afp.len = fragLen;
                    afp.score = resScore * fragLen * (1.0 - (rmsd_tmp / cur_badRmsd) * (rmsd_tmp / cur_badRmsd));
                    initial_afps.push_back(afp);
                }
            }
        }

        // Step 2: Merge diagonal AFPs
        int max_diagonal_idx = xlen + ylen + 1;
        std::vector<std::vector<USBCAT_AFP>> diagonals(max_diagonal_idx);
        for (size_t k = 0; k < initial_afps.size(); k++)
        {
            diagonals[initial_afps[k].i - initial_afps[k].j + ylen].push_back(initial_afps[k]);
        }

        std::vector<USBCAT_AFP> merged_afps;
        int max_merge_len = std::min(xlen, ylen);
        double **r1_merge, **r2_merge;
        NewArray(&r1_merge, max_merge_len, 3);
        NewArray(&r2_merge, max_merge_len, 3);

        for (int d = 0; d < max_diagonal_idx; d++)
        {
            if (diagonals[d].empty())
                continue;
            std::vector<USBCAT_AFP> &group = diagonals[d];

            std::sort(group.begin(), group.end(), [](const USBCAT_AFP &a, const USBCAT_AFP &b)
                      { return a.i < b.i; });

            int n_group = group.size();
            std::vector<bool> invalid(n_group, false);
            for (int idx = 0; idx < n_group; idx++)
            {
                if (invalid[idx])
                    continue;
                USBCAT_AFP curr = group[idx];
                for (int nxt_idx = idx + 1; nxt_idx < n_group; nxt_idx++)
                {
                    USBCAT_AFP nxt = group[nxt_idx];
                    if (nxt.i > curr.i + curr.len)
                        break;

                    if (nxt.i + nxt.len > curr.i + curr.len)
                    {
                        int new_len = (nxt.i + nxt.len) - curr.i;
                        for (int k = 0; k < new_len; k++)
                        {
                            r1_merge[k][0] = xa[curr.i + k][0];
                            r1_merge[k][1] = xa[curr.i + k][1];
                            r1_merge[k][2] = xa[curr.i + k][2];
                            r2_merge[k][0] = ya[curr.j + k][0];
                            r2_merge[k][1] = ya[curr.j + k][1];
                            r2_merge[k][2] = ya[curr.j + k][2];
                        }

                        double rms_sum_sq, t_tmp[3], u_tmp[3][3];
                        Kabsch(r1_merge, r2_merge, new_len, 0, &rms_sum_sq, t_tmp, u_tmp);
                        double rmsd_tmp = std::sqrt(rms_sum_sq / new_len);

                        if (rmsd_tmp < cur_rmsdCut)
                        {
                            curr.len = new_len;
                            curr.score = resScore * new_len * (1.0 - (rmsd_tmp / cur_badRmsd) * (rmsd_tmp / cur_badRmsd));
                            invalid[nxt_idx] = true;
                        }
                    }
                }
                merged_afps.push_back(curr);
            }
        }
        DeleteArray(&r1_merge, max_merge_len);
        DeleteArray(&r2_merge, max_merge_len);

        std::sort(merged_afps.begin(), merged_afps.end(), [](const USBCAT_AFP &a, const USBCAT_AFP &b)
                  {
            if (a.i == b.i) return a.j < b.j;
            return a.i < b.i; });

        int n_afps = merged_afps.size();
        std::vector<int> ret_b1, ret_b2;
        if (n_afps == 0)
            return std::make_pair(ret_b1, ret_b2);

        // Step 3 & 4: Dual Dynamic Programming and Domain Splitting
        std::vector<int> afp_aft_index(xlen * ylen, -1);
        std::vector<int> afp_bef_index(xlen * ylen, -1);

        std::vector<std::vector<std::pair<int, int>>> i_to_j(xlen);
        for (int m = 0; m < n_afps; m++)
        {
            i_to_j[merged_afps[m].i].push_back(std::make_pair(merged_afps[m].j, m));
        }

        for (int i_val = 0; i_val < xlen; i_val++)
        {
            if (i_to_j[i_val].empty())
                continue;
            for (size_t p = 0; p < i_to_j[i_val].size(); p++)
            {
                int j_val = i_to_j[i_val][p].first;
                afp_aft_index[i_val * ylen + j_val] = i_to_j[i_val][p].second;
                afp_bef_index[i_val * ylen + j_val] = i_to_j[i_val][p].second;
            }
            int curr_bef = -1;
            for (int j_val = 0; j_val < ylen; j_val++)
            {
                if (afp_bef_index[i_val * ylen + j_val] != -1)
                    curr_bef = afp_bef_index[i_val * ylen + j_val];
                else
                    afp_bef_index[i_val * ylen + j_val] = curr_bef;
            }
            int curr_aft = -1;
            for (int j_val = ylen - 1; j_val >= 0; j_val--)
            {
                if (afp_aft_index[i_val * ylen + j_val] != -1)
                    curr_aft = afp_aft_index[i_val * ylen + j_val];
                else
                    afp_aft_index[i_val * ylen + j_val] = curr_aft;
            }
        }

        auto get_dvar = [&](const USBCAT_AFP &prv, const USBCAT_AFP &curr) -> double
        {
            double rms_sq = 0;
            for (int i_idx = 0; i_idx < fragLen; i_idx++)
            {
                for (int j_idx = 0; j_idx < fragLen; j_idx++)
                {
                    double dist1, dist2;
                    int idx1_a = curr.i + i_idx, idx1_b = prv.i + j_idx;
                    if (idx1_a >= idx1_b)
                        dist1 = disTable1[idx1_b][idx1_a - idx1_b];
                    else
                        dist1 = disTable1[idx1_a][idx1_b - idx1_a];

                    int idx2_a = curr.j + i_idx, idx2_b = prv.j + j_idx;
                    if (idx2_a >= idx2_b)
                        dist2 = disTable2[idx2_b][idx2_a - idx2_b];
                    else
                        dist2 = disTable2[idx2_a][idx2_b - idx2_a];

                    rms_sq += (dist1 - dist2) * (dist1 - dist2);
                }
            }
            if (rms_sq > afp_dis_cut)
                return 1e9;
            return std::sqrt(rms_sq / (fragLen * fragLen));
        };

        auto calc_block_rmsd = [&](const std::vector<USBCAT_AFP> &afp_list) -> double
        {
            std::vector<int> r1, r2;
            for (size_t a = 0; a < afp_list.size(); a++)
            {
                for (int l = 0; l < afp_list[a].len; l++)
                {
                    r1.push_back(afp_list[a].i + l);
                    r2.push_back(afp_list[a].j + l);
                }
            }
            int n = r1.size();
            if (n < 3)
                return 0.0;
            double **p1;
            NewArray(&p1, n, 3);
            double **p2;
            NewArray(&p2, n, 3);
            for (int i = 0; i < n; i++)
            {
                p1[i][0] = xa[r1[i]][0];
                p1[i][1] = xa[r1[i]][1];
                p1[i][2] = xa[r1[i]][2];
                p2[i][0] = ya[r2[i]][0];
                p2[i][1] = ya[r2[i]][1];
                p2[i][2] = ya[r2[i]][2];
            }
            double rms_sq_sum, t_tmp[3], u_tmp[3][3];
            Kabsch(p1, p2, n, 0, &rms_sq_sum, t_tmp, u_tmp);
            DeleteArray(&p1, n);
            DeleteArray(&p2, n);
            return std::sqrt(rms_sq_sum / n);
        };

        struct Region
        {
            int s1, e1, s2, e2;
        };

        std::vector<double> sco(n_afps);
        std::vector<int> twi(n_afps, 0);
        std::vector<int> pre(n_afps, -1);
        for (int m = 0; m < n_afps; m++)
            sco[m] = merged_afps[m].score;

        for (int m = 0; m < n_afps; m++)
        {
            int curr_i = merged_afps[m].i;
            int curr_j = merged_afps[m].j;
            int a3 = curr_i - fragLen;
            int a2 = std::max(0, a3 - misCut);
            int a1 = std::max(0, curr_i - maxGapFrag);
            int b3 = curr_j - fragLen;
            int b2 = std::max(0, b3 - misCut);
            int b1 = std::max(0, curr_j - maxGapFrag);

            std::vector<int> valid_prevs;
            for (int st = 0; st < 2; st++)
            {
                int a_s, a_e, b_s, b_e;
                if (st == 0)
                {
                    a_s = std::max(a1, 0);
                    a_e = std::min(a3, xlen - 1);
                    b_s = std::max(b2, 0);
                    b_e = std::min(b3, ylen - 1);
                }
                else
                {
                    a_s = std::max(a2, 0);
                    a_e = std::min(a3, xlen - 1);
                    b_s = std::max(b1, 0);
                    b_e = std::min(b2 - 1, ylen - 1);
                }

                if (b_s >= ylen || b_e < 0)
                    continue;
                for (int prev_i = a_s; prev_i <= a_e; prev_i++)
                {
                    int s1 = afp_aft_index[prev_i * ylen + b_s];
                    int s2 = afp_bef_index[prev_i * ylen + b_e];
                    if (s1 != -1 && s2 != -1 && s1 <= s2)
                        for (int s = s1; s <= s2; s++)
                            valid_prevs.push_back(s);
                }
            }

            double curr_sco = merged_afps[m].score;
            for (size_t v = 0; v < valid_prevs.size(); v++)
            {
                int prev = valid_prevs[v];
                int prev_twi = twi[prev];
                if (prev_twi > max_twists)
                    continue;

                int gap_i = curr_i - (merged_afps[prev].i + merged_afps[prev].len);
                int gap_j = curr_j - (merged_afps[prev].j + merged_afps[prev].len);
                int m_gap = std::max(gap_i, gap_j);

                double gp = 0.0;
                int m_mis = 0;
                if (gap_i < 0 || gap_j < 0)
                    m_mis = (gap_i < gap_j) ? -gap_i : -gap_j;
                gp = gap_ext * m_mis;
                if (m_gap > 0)
                    gp += gap_ext * m_gap;
                if (gp < max_penalty)
                    gp = max_penalty;

                double rms_sq = 0;
                for (int k = 0; k < fragLen; k++)
                {
                    for (int l = 0; l < fragLen; l++)
                    {
                        double dist1, dist2;
                        int idx1_a = curr_i + k, idx1_b = merged_afps[prev].i + l;
                        if (idx1_a >= idx1_b)
                            dist1 = disTable1[idx1_b][idx1_a - idx1_b];
                        else
                            dist1 = disTable1[idx1_a][idx1_b - idx1_a];

                        int idx2_a = curr_j + k, idx2_b = merged_afps[prev].j + l;
                        if (idx2_a >= idx2_b)
                            dist2 = disTable2[idx2_b][idx2_a - idx2_b];
                        else
                            dist2 = disTable2[idx2_a][idx2_b - idx2_a];

                        rms_sq += (dist1 - dist2) * (dist1 - dist2);
                    }
                }

                double tp = 0.0;
                int is_twist = 0;
                if (rms_sq >= afp_dis_cut)
                {
                    tp = twist_pen;
                    is_twist = 1;
                }
                else
                {
                    double dvar = std::sqrt(rms_sq / (fragLen * fragLen));
                    if (dvar > disCut - disSmooth)
                        tp = twist_pen * std::sqrt((dvar - disCut + disSmooth) / disSmooth);
                }

                if (prev_twi + is_twist > max_twists)
                    continue;

                double stmp = sco[prev] + curr_sco + tp + gp;
                if (stmp > sco[m])
                {
                    sco[m] = stmp;
                    pre[m] = prev;
                    twi[m] = prev_twi + is_twist;
                }
            }
        }

        int best_m = 0;
        for (int m = 1; m < n_afps; m++)
            if (sco[m] > sco[best_m])
                best_m = m;

        std::vector<int> path;
        int curr_m = best_m;
        while (curr_m != -1)
        {
            path.push_back(curr_m);
            curr_m = pre[curr_m];
        }
        std::reverse(path.begin(), path.end());

        if (path.empty())
            return std::make_pair(ret_b1, ret_b2);

        struct Block
        {
            std::vector<USBCAT_AFP> afps;
            std::vector<double> dvars;
        };
        std::vector<Block> candidate_blocks;
        Block curr_block;
        curr_block.afps.push_back(merged_afps[path[0]]);
        curr_block.dvars.push_back(0.0);

        for (size_t k = 1; k < path.size(); k++)
        {
            USBCAT_AFP curr = merged_afps[path[k]];
            USBCAT_AFP prv = merged_afps[path[k - 1]];
            double dvar = get_dvar(prv, curr);

            if (dvar >= disCut)
            {
                candidate_blocks.push_back(curr_block);
                curr_block.afps.clear();
                curr_block.dvars.clear();
                curr_block.afps.push_back(curr);
                curr_block.dvars.push_back(0.0);
            }
            else
            {
                curr_block.afps.push_back(curr);
                curr_block.dvars.push_back(dvar);
            }
        }
        if (!curr_block.afps.empty())
            candidate_blocks.push_back(curr_block);

        bool splitted = true;
        while (splitted && candidate_blocks.size() < (size_t)(max_twists + 1))
        {
            splitted = false;
            double max_rmsd = 0.0;
            int target_b = -1;

            for (size_t b = 0; b < candidate_blocks.size(); b++)
            {
                if (candidate_blocks[b].afps.size() > 2)
                {
                    double cur_rmsd = calc_block_rmsd(candidate_blocks[b].afps);
                    if (cur_rmsd > max_rmsd)
                    {
                        max_rmsd = cur_rmsd;
                        target_b = b;
                    }
                }
            }

            if (max_rmsd >= cur_local_badRmsd && target_b != -1)
            {
                double max_t = 0;
                int cut_idx = 0;
                for (size_t i = 1; i < candidate_blocks[target_b].afps.size(); i++)
                {
                    if (candidate_blocks[target_b].dvars[i] > max_t)
                    {
                        max_t = candidate_blocks[target_b].dvars[i];
                        cut_idx = i;
                    }
                }

                if (cut_idx > 0)
                {
                    Block right_blk;
                    right_blk.afps.assign(candidate_blocks[target_b].afps.begin() + cut_idx, candidate_blocks[target_b].afps.end());
                    right_blk.dvars.assign(candidate_blocks[target_b].dvars.begin() + cut_idx, candidate_blocks[target_b].dvars.end());
                    right_blk.dvars[0] = 0.0;
                    candidate_blocks[target_b].afps.erase(candidate_blocks[target_b].afps.begin() + cut_idx, candidate_blocks[target_b].afps.end());
                    candidate_blocks[target_b].dvars.erase(candidate_blocks[target_b].dvars.begin() + cut_idx, candidate_blocks[target_b].dvars.end());
                    candidate_blocks.insert(candidate_blocks.begin() + target_b + 1, right_blk);
                    splitted = true;
                }
            }
        }

        for (int b = 0; b < (int)candidate_blocks.size(); b++)
        {
            if (candidate_blocks[b].afps.size() <= 1)
            {
                int e1 = (b < (int)candidate_blocks.size() - 1) ? candidate_blocks[b + 1].afps.front().i : xlen;
                int e2 = (b < (int)candidate_blocks.size() - 1) ? candidate_blocks[b + 1].afps.front().j : ylen;
                int b1 = (b > 0) ? candidate_blocks[b - 1].afps.back().i + candidate_blocks[b - 1].afps.back().len : 0;
                int b2 = (b > 0) ? candidate_blocks[b - 1].afps.back().j + candidate_blocks[b - 1].afps.back().len : 0;
                int span = std::min(e1 - b1, e2 - b2);
                if (span < 2 * fragLen)
                {
                    candidate_blocks.erase(candidate_blocks.begin() + b);
                    b--;
                }
            }
        }

        bool merged = true;
        while (merged && candidate_blocks.size() > 1)
        {
            merged = false;
            double min_rmsd = 1e9;
            int min_b = -1;
            for (size_t b = 0; b < candidate_blocks.size() - 1; b++)
            {
                std::vector<USBCAT_AFP> temp_merged = candidate_blocks[b].afps;
                temp_merged.insert(temp_merged.end(), candidate_blocks[b + 1].afps.begin(), candidate_blocks[b + 1].afps.end());
                double cur_rmsd = calc_block_rmsd(temp_merged);
                if (cur_rmsd < min_rmsd)
                {
                    min_rmsd = cur_rmsd;
                    min_b = b;
                }
            }

            if (min_rmsd < cur_local_badRmsd && min_b != -1)
            {
                candidate_blocks[min_b].afps.insert(candidate_blocks[min_b].afps.end(), candidate_blocks[min_b + 1].afps.begin(), candidate_blocks[min_b + 1].afps.end());
                candidate_blocks.erase(candidate_blocks.begin() + min_b + 1);
                merged = true;
            }
        }

        std::vector<Region> usbcat_domains;
        int last_i = 0, last_j = 0;
        for (size_t b = 0; b < candidate_blocks.size(); b++)
        {
            int b_s1 = -1, b_e1 = -1, b_s2 = -1, b_e2 = -1;
            for (size_t a = 0; a < candidate_blocks[b].afps.size(); a++)
            {
                USBCAT_AFP afp = candidate_blocks[b].afps[a];
                int skip = std::max(std::max(last_i - afp.i, last_j - afp.j), 0);
                if (skip >= afp.len)
                    continue;

                int eff_i = afp.i + skip;
                int eff_j = afp.j + skip;
                int eff_L = afp.len - skip;
                if (b_s1 == -1)
                {
                    b_s1 = eff_i;
                    b_s2 = eff_j;
                }
                b_e1 = eff_i + eff_L;
                b_e2 = eff_j + eff_L;
                last_i = b_e1;
                last_j = b_e2;
            }
            if (b_s1 != -1)
            {
                if (b_e1 - b_s1 >= 4 && b_e2 - b_s2 >= 4)
                {
                    Region r = {b_s1, b_e1, b_s2, b_e2};
                    usbcat_domains.push_back(r);
                }
            }
        }

        if (usbcat_domains.empty())
            return std::make_pair(ret_b1, ret_b2);

        ret_b1.push_back(0);
        ret_b2.push_back(0);
        for (size_t k = 0; k < usbcat_domains.size() - 1; k++)
        {
            ret_b1.push_back((usbcat_domains[k].e1 + usbcat_domains[k + 1].s1) / 2);
            ret_b2.push_back((usbcat_domains[k].e2 + usbcat_domains[k + 1].s2) / 2);
        }
        ret_b1.push_back(xlen);
        ret_b2.push_back(ylen);

        return std::make_pair(ret_b1, ret_b2);
    };

    auto bounds_usbcat = generate_bounds(3.0, 4.0, 4.0);
    auto bounds_strict = generate_bounds(2.0, 3.0, 2.0);

    std::vector<std::pair<std::vector<int>, std::vector<int>>> all_bounds;
    all_bounds.push_back(bounds_usbcat);
    if (bounds_strict.first != bounds_usbcat.first || bounds_strict.second != bounds_usbcat.second)
    {
        all_bounds.push_back(bounds_strict);
    }

    for (size_t b_idx = 0; b_idx < all_bounds.size(); b_idx++)
    {
        std::vector<int> &bounds1 = all_bounds[b_idx].first;
        std::vector<int> &bounds2 = all_bounds[b_idx].second;

        // Skip if only one interval (block) is generated
        if (bounds1.size() <= 2)
            continue;

        // =========================================================================
        // Greedy Dynamic Hinge Budgeting based on dRMSD (Strain Energy)
        // =========================================================================
        int num_blocks = bounds1.size() - 1;

        // 1. Define execution node for out-of-order processing
        struct BlockMeta
        {
            int original_idx;
            double drmsd;
            int L1_sub, L2_sub;
        };
        std::vector<BlockMeta> block_queue;

        // 2. Calculate proxy dRMSD for each block to evaluate internal strain energy
        for (int k = 0; k < num_blocks; k++)
        {
            int x_s = bounds1[k], x_e = bounds1[k + 1];
            int y_s = bounds2[k], y_e = bounds2[k + 1];
            int L1_sub = x_e - x_s;
            int L2_sub = y_e - y_s;
            int min_L = std::min(L1_sub, L2_sub);

            double block_drmsd = 0.0;
            // Only calculate if the block is long enough
            if (min_L >= 2 * fragLen)
            {
                double rms_sq = 0.0;
                int count = 0;
                for (int i = 0; i < min_L; i++)
                {
                    // DOUBLE OPTIMIZATION:
                    // 1. Start from i + 2 to skip adjacent amino acids (peptide bond noise).
                    // 2. Cap j at i + max_dist_window to SAFELY reuse precomputed disTable!
                    //    This perfectly evaluates "local" strain energy and reduces time complexity to O(N).
                    int j_end = std::min((int)min_L, i + max_dist_window);

                    for (int j = i + 2; j < j_end; j++)
                    {
                        // Directly query the precomputed distance tables
                        double d1 = disTable1[x_s + i][j - i];
                        double d2 = disTable2[y_s + i][j - i];
                        rms_sq += (d1 - d2) * (d1 - d2);
                        count++;
                    }
                }
                if (count > 0)
                    block_drmsd = std::sqrt(rms_sq / count);
            }
            block_queue.push_back({k, block_drmsd, L1_sub, L2_sub});
        }

        // 3. Sort blocks by dRMSD descending (most twisted blocks get priority)
        if (hinge_set)
        {
            std::sort(block_queue.begin(), block_queue.end(), [](const BlockMeta &a, const BlockMeta &b)
                      { return a.drmsd > b.drmsd; });
        }

        // 4. Initialize global hinge pool
        // N blocks intrinsically use N-1 cut points, remaining hinges = hinge_opt + 1 - N
        int remaining_hinges = hinge_set ? std::max(0, hinge_opt + 1 - num_blocks) : 0;

        // Structure to store out-of-order execution results
        struct BlockResult
        {
            bool valid;
            double t0[3];
            double u0[3][3];
            std::string seqM, seqxA, seqyA;
            std::vector<std::vector<double>> tu_vec;
            BlockResult() : valid(false) {}
        };
        std::vector<BlockResult> block_results(num_blocks);

        // 5. Greedy execution: Allocate all available hinges to the current most twisted block
        for (size_t q = 0; q < block_queue.size(); q++)
        {
            int k = block_queue[q].original_idx;
            int L1_sub = block_queue[q].L1_sub;
            int L2_sub = block_queue[q].L2_sub;
            int x_s = bounds1[k], y_s = bounds2[k];

            // Skip invalid or too short blocks
            if (L1_sub < 3 || L2_sub < 3)
                continue;

            // CORE LOGIC: Pass all remaining budget to the current block
            int local_hinge_opt = 0;
            if (remaining_hinges > 0 && std::min(L1_sub, L2_sub) >= 2 * fragLen)
            {
                local_hinge_opt = remaining_hinges;
            }

            // Reuse variables from the original logic for sub-block allocation
            double **xa_sub, **ya_sub;
            NewArray(&xa_sub, L1_sub, 3);
            NewArray(&ya_sub, L2_sub, 3);
            char *seqx_sub = new char[L1_sub + 1];
            char *secx_sub = new char[L1_sub + 1];
            char *seqy_sub = new char[L2_sub + 1];
            char *secy_sub = new char[L2_sub + 1];

            for (int i = 0; i < L1_sub; i++)
            {
                xa_sub[i][0] = xa[x_s + i][0];
                xa_sub[i][1] = xa[x_s + i][1];
                xa_sub[i][2] = xa[x_s + i][2];
                seqx_sub[i] = seqx[x_s + i];
                secx_sub[i] = secx[x_s + i];
            }
            seqx_sub[L1_sub] = '\0';
            secx_sub[L1_sub] = '\0';

            for (int i = 0; i < L2_sub; i++)
            {
                ya_sub[i][0] = ya[y_s + i][0];
                ya_sub[i][1] = ya[y_s + i][1];
                ya_sub[i][2] = ya[y_s + i][2];
                seqy_sub[i] = seqy[y_s + i];
                secy_sub[i] = secy[y_s + i];
            }
            seqy_sub[L2_sub] = '\0';
            secy_sub[L2_sub] = '\0';

            bool force_fast_opt = (std::min(L1_sub, L2_sub) > 1500) ? true : fast_opt;
            double TM_best_max = -1.0;

            // Try both secondary structure configurations
            for (int cur_ss_opt = 0; cur_ss_opt <= 1; cur_ss_opt++)
            {
                FlexAlignResult cur_res;
                execute_flexalign_with_fallback(
                    xa_sub, ya_sub, seqx_sub, seqy_sub, secx_sub, secy_sub,
                    L1_sub, L2_sub, local_sequence, Lnorm_ass, d0_scale,
                    i_opt, a_opt, u_opt, d_opt, force_fast_opt,
                    mol_type, local_hinge_opt, cur_ss_opt, cur_res);

                double cur_max_TM = (cur_res.TM1 > cur_res.TM2) ? cur_res.TM1 : cur_res.TM2;
                if (cur_max_TM > TM_best_max)
                {
                    TM_best_max = cur_max_TM;
                    for (int a = 0; a < 3; a++)
                    {
                        block_results[k].t0[a] = cur_res.t0[a];
                        for (int b = 0; b < 3; b++)
                            block_results[k].u0[a][b] = cur_res.u0[a][b];
                    }
                    block_results[k].seqM = cur_res.seqM;
                    block_results[k].seqxA = cur_res.seqxA;
                    block_results[k].seqyA = cur_res.seqyA;
                    block_results[k].tu_vec = cur_res.tu_vec;
                    block_results[k].valid = true;
                }
            }

            // Deduct actually consumed hinges from the global budget
            if (block_results[k].valid && !block_results[k].tu_vec.empty())
            {
                int consumed_hinges = block_results[k].tu_vec.size() - 1;
                if (consumed_hinges > 0)
                {
                    remaining_hinges -= consumed_hinges;
                    if (remaining_hinges < 0)
                        remaining_hinges = 0; // Guard against negative budget
                }
            }

            // Clean up sub-block memory
            DeleteArray(&xa_sub, L1_sub);
            DeleteArray(&ya_sub, L2_sub);
            delete[] seqx_sub;
            delete[] seqy_sub;
            delete[] secx_sub;
            delete[] secy_sub;
        }

        // 6. Chronological Stitching: Reassemble results in spatial sequence order
        std::string cur_global_seqM = "", cur_global_seqxA = "", cur_global_seqyA = "";
        std::vector<std::vector<double>> cur_tu_vec;
        std::vector<int> cur_global_res_tu(xlen, -1);

        for (int k = 0; k < num_blocks; k++)
        {
            int L1_sub = bounds1[k + 1] - bounds1[k];
            int L2_sub = bounds2[k + 1] - bounds2[k];

            if (!block_results[k].valid)
            {
                // Fill gaps if block was invalid or bypassed
                for (int i = 0; i < L1_sub; i++)
                {
                    cur_global_seqxA += seqx[bounds1[k] + i];
                    cur_global_seqyA += '-';
                    cur_global_seqM += ' ';
                }
                for (int i = 0; i < L2_sub; i++)
                {
                    cur_global_seqxA += '-';
                    cur_global_seqyA += seqy[bounds2[k] + i];
                    cur_global_seqM += ' ';
                }
                continue;
            }

            BlockResult &res = block_results[k];
            if (res.tu_vec.empty())
            {
                std::vector<double> tu_tmp(12);
                t_u2tu(res.t0, res.u0, tu_tmp);
                res.tu_vec.push_back(tu_tmp);
            }

            int base_tu_idx = cur_tu_vec.size();
            for (size_t m = 0; m < res.tu_vec.size(); m++)
                cur_tu_vec.push_back(res.tu_vec[m]);

            int rx = bounds1[k];
            int current_global_idx = base_tu_idx;

            for (size_t i = 0; i < res.seqxA.length(); i++)
            {
                char c = res.seqM[i];
                if (c != ' ' && c != '.' && c != ':')
                {
                    int local_hinge_idx = -1;
                    if (c >= '0' && c <= '9')
                        local_hinge_idx = c - '0';
                    else if (c >= 'a' && c <= 'z')
                        local_hinge_idx = c - 'a' + 10;
                    else if (c >= 'A' && c <= 'Z')
                        local_hinge_idx = c - 'A' + 36;

                    if (local_hinge_idx >= 0 && local_hinge_idx < res.tu_vec.size())
                    {
                        current_global_idx = base_tu_idx + local_hinge_idx;
                    }
                }

                if (res.seqxA[i] != '-')
                {
                    cur_global_res_tu[rx] = current_global_idx;
                    rx++;
                }

                if (res.seqxA[i] != '-' && res.seqyA[i] != '-')
                {
                    if (c != ' ' && c != '.' && c != ':')
                    {
                        char global_c;
                        if (current_global_idx < 10)
                            global_c = '0' + current_global_idx;
                        else if (current_global_idx < 36)
                            global_c = 'a' + (current_global_idx - 10);
                        else if (current_global_idx < 62)
                            global_c = 'A' + (current_global_idx - 36);
                        else
                            global_c = '*';
                        res.seqM[i] = global_c;
                    }
                    else
                    {
                        res.seqM[i] = c;
                    }
                }
                else
                {
                    res.seqM[i] = ' ';
                }
            }

            cur_global_seqM += res.seqM;
            cur_global_seqxA += res.seqxA;
            cur_global_seqyA += res.seqyA;
        }

        // Step 7: Recalculate global metrics correctly for current DP boundary
        double dummy_D0_MIN, dummy_Lnorm, dummy_d0_search;
        double cur_d0A, cur_d0B, cur_d0a, cur_d0u = 0.0;

        parameter_set4final(ylen, dummy_D0_MIN, dummy_Lnorm, cur_d0A, dummy_d0_search, mol_type);
        parameter_set4final(xlen, dummy_D0_MIN, dummy_Lnorm, cur_d0B, dummy_d0_search, mol_type);
        parameter_set4final((xlen + ylen) * 0.5, dummy_D0_MIN, dummy_Lnorm, cur_d0a, dummy_d0_search, mol_type);

        if (u_opt)
        {
            parameter_set4final(Lnorm_ass, dummy_D0_MIN, dummy_Lnorm, cur_d0u, dummy_d0_search, mol_type);
        }

        double cur_TM1 = 0.0, cur_TM2 = 0.0, cur_TM3 = 0.0, cur_TM4 = 0.0, cur_TM5 = 0.0;
        double cur_rmsd0 = 0.0, cur_Liden = 0.0;
        int cur_n_ali8 = 0, cur_n_ali = 0;
        std::vector<double> cur_do_vec;

        int i_res = 0, j_res = 0;
        for (size_t r = 0; r < cur_global_seqxA.length(); r++)
        {
            bool x_valid = (cur_global_seqxA[r] != '-');
            bool y_valid = (cur_global_seqyA[r] != '-');

            if (x_valid && y_valid)
            {
                int matrix_idx = cur_global_res_tu[i_res];
                if (matrix_idx >= 0 && matrix_idx < cur_tu_vec.size())
                {
                    double t_k[3], u_k[3][3];
                    tu2t_u(cur_tu_vec[matrix_idx], t_k, u_k);

                    double x_rot[3];
                    transform(t_k, u_k, xa[i_res], x_rot);
                    double dist2 = dist(x_rot, ya[j_res]);
                    double d = std::sqrt(dist2);

                    cur_TM2 += 1.0 / (1.0 + dist2 / (cur_d0B * cur_d0B));
                    cur_TM1 += 1.0 / (1.0 + dist2 / (cur_d0A * cur_d0A));
                    if (a_opt)
                        cur_TM3 += 1.0 / (1.0 + dist2 / (cur_d0a * cur_d0a));
                    if (u_opt)
                        cur_TM4 += 1.0 / (1.0 + dist2 / (cur_d0u * cur_d0u));
                    if (d_opt)
                        cur_TM5 += 1.0 / (1.0 + dist2 / (d0_scale * d0_scale));

                    cur_n_ali++;
                    cur_do_vec.push_back(d);

                    if (d <= d0_out)
                    {
                        cur_rmsd0 += dist2;
                        cur_n_ali8++;
                        if (seqx[i_res] == seqy[j_res])
                            cur_Liden += 1.0;
                    }
                }
                else
                {
                    cur_do_vec.push_back(-1);
                }
            }
            else
            {
                cur_do_vec.push_back(-1);
            }

            if (x_valid)
                i_res++;
            if (y_valid)
                j_res++;
        }

        cur_TM2 /= xlen;
        cur_TM1 /= ylen;
        if (a_opt)
            cur_TM3 /= (xlen + ylen) * 0.5;
        if (u_opt)
            cur_TM4 /= Lnorm_ass;
        if (d_opt)
            cur_TM5 /= ylen;
        if (cur_n_ali8 > 0)
            cur_rmsd0 = std::sqrt(cur_rmsd0 / cur_n_ali8);
        else
            cur_rmsd0 = 0.0;

        double cur_global_max_TM = (cur_TM1 > cur_TM2) ? cur_TM1 : cur_TM2;

        if (cur_global_max_TM > best_global_max_TM)
        {
            best_global_max_TM = cur_global_max_TM;
            best_tu_vec = cur_tu_vec;
            best_TM1 = cur_TM1;
            best_TM2 = cur_TM2;
            best_TM3 = cur_TM3;
            best_TM4 = cur_TM4;
            best_TM5 = cur_TM5;
            best_rmsd0 = cur_rmsd0;
            best_Liden = cur_Liden;
            best_TM_ali = cur_TM1;
            best_rmsd_ali = cur_rmsd0;
            best_L_ali = cur_n_ali;
            best_n_ali = cur_n_ali;
            best_n_ali8 = cur_n_ali8;
            best_seqM = cur_global_seqM;
            best_seqxA = cur_global_seqxA;
            best_seqyA = cur_global_seqyA;
            best_do_vec = cur_do_vec;
            best_d0A = cur_d0A;
            best_d0B = cur_d0B;
            best_d0a = cur_d0a;
            best_d0u = cur_d0u;

            if (!best_tu_vec.empty())
            {
                tu2t_u(best_tu_vec[0], best_t0, best_u0);
            }
        }
    }

    // Safety check
    if (best_global_max_TM < 0)
        return 0;

    // Output best values back to the reference parameters
    TM1 = best_TM1;
    TM2 = best_TM2;
    TM3 = best_TM3;
    TM4 = best_TM4;
    TM5 = best_TM5;
    rmsd0 = best_rmsd0;
    Liden = best_Liden;
    TM_ali = best_TM_ali;
    rmsd_ali = best_rmsd_ali;
    L_ali = best_L_ali;
    n_ali = best_n_ali;
    n_ali8 = best_n_ali8;
    seqM = best_seqM;
    seqxA = best_seqxA;
    seqyA = best_seqyA;
    do_vec = best_do_vec;
    tu_vec = best_tu_vec;
    d0A = best_d0A;
    d0B = best_d0B;
    d0a = best_d0a;
    d0u = best_d0u;

    for (int a = 0; a < 3; a++)
    {
        t0[a] = best_t0[a];
        for (int b = 0; b < 3; b++)
            u0[a][b] = best_u0[a][b];
    }

    return tu_vec.size();
}
#endif
