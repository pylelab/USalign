#ifndef TMalign_HwRMSD_h
#define TMalign_HwRMSD_h 1
#include "NWalign.h"
#include "se.h"

double Kabsch_Superpose(double **r1, double **r2, double **xt,
    double **xa, double **ya, int xlen, int ylen, int invmap[],
    int& L_ali, double t[3], double u[3][3], const int mol_type)
{
    L_ali = 0;
    int i, j;
    for (j = 0; j<ylen; j++)
    {
        i = invmap[j];
        if (i >= 0)
        {
            r1[L_ali][0]  = xa[i][0];
            r1[L_ali][1]  = xa[i][1];
            r1[L_ali][2]  = xa[i][2];

            r2[L_ali][0]  = ya[j][0];
            r2[L_ali][1]  = ya[j][1];
            r2[L_ali][2]  = ya[j][2];

            L_ali++;
        }
        else if (i != -1) PrintErrorAndQuit("Wrong map!\n");
    }

    double RMSD = 0;
    Kabsch(r1, r2, L_ali, 1, &RMSD, t, u);
    RMSD = sqrt( RMSD/(1.0*L_ali) );

    //for (i=0;i<3;i++)
    //{
        //cout<<t[i];
        //for (j=0;j<3;j++)
            //cout<<'\t'<<u[i][j];
        //cout<<endl;
    //}

    for (i=0; i<xlen; i++)
    {
        xt[i][0] = xa[i][0];
        xt[i][1] = xa[i][1];
        xt[i][2] = xa[i][2];
    }
    do_rotation(xa, xt, xlen, t,u);
    return RMSD;
}

int HwRMSD_main(double **xa, double **ya, const char *seqx, const char *seqy,
    double t0[3], double u0[3][3], double &TM1, double &TM2,
    double &TM3, double &TM4, double &TM5, double &d0_0, double &TM_0,
    double &d0A, double &d0B, double &d0u, double &d0a, double &d0_out,
    string &seqM, string &seqxA, string &seqyA, double &rmsd0, int &L_ali,
    double &Liden, double &TM_ali, double &rmsd_ali, int &n_ali, 
    int &n_ali8, const int xlen, const int ylen, vector<string> sequence,
    const double Lnorm_ass, const double d0_scale, const bool i_opt,
    const bool I_opt, const int a_opt, const bool u_opt, const bool d_opt,
    const int mol_type, const int glocal=0)
{
    if (!I_opt && !i_opt)
    {
        sequence.push_back("");
        sequence.push_back("");
        NWalign(seqx, seqy, xlen, ylen, sequence[0], sequence[1],
                mol_type, glocal);
    }

    /***********************/
    /* allocate memory     */
    /***********************/
    double t[3], u[3][3]; //Kabsch translation vector and rotation matrix
    double **xt;          //for saving the superposed version of r_1 or xtm
    double **r1, **r2;    // for Kabsch rotation
    int minlen = min(xlen, ylen);
    NewArray(&xt, xlen, 3);
    NewArray(&r1, minlen, 3);
    NewArray(&r2, minlen, 3);
    int *invmap = new int[ylen+1];

    /*****************************/
    /*    get initial alignment  */
    /*****************************/
    for (int j = 0; j < ylen; j++)// Set aligned position to be "-1"
        invmap[j] = -1;

    int i1 = -1, i2 = -1;
    int L = min(sequence[0].size(), sequence[1].size());
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
    Kabsch_Superpose(r1, r2, xt, xa, ya, xlen, ylen, invmap,
        L_ali, t, u, mol_type);

    se_main(xt, ya, seqx, seqy, TM1, TM2, TM3, TM4, TM5,
        d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
        seqM, seqxA, seqyA, rmsd0, L_ali, Liden,
        TM_ali, rmsd_ali, n_ali, n_ali8, xlen, ylen, sequence,
        Lnorm_ass, d0_scale, I_opt, a_opt, u_opt, d_opt, mol_type);

    /* clean up */
    delete [] invmap;
    DeleteArray(&xt, xlen);
    DeleteArray(&r1, minlen);
    DeleteArray(&r2, minlen);
    return 0;
}
#endif
