/* Isolating and simplifying the SOIalign_main function so that I can create a 
 * python wrapper implementation for monomeric protein pair alignment using SOI
 * alignment methods (sNS and fNS) */

#include "SOIalign.h"	// 
#include "TMalign.h"	//

#include "param_set.h"	// 
#include "basic_fun.h"	//

#include "NW.h"		//

#include "Kabsch.h"	//

#include "NWalign.h"	//
#include "BLOSUM.h"	//

#include <math.h>	//

/////////////////////////////////////////////////////////////////////////////
/* TMalign.h code */
/////////////////////////////////////////////////////////////////////////////

double standard_TMscore(double **r1, double **r2, double **xtm, double **ytm,
    double **xt, double **x, double **y, int xlen, int ylen, int invmap[],
    int& L_ali, double& RMSD, double D0_MIN, double Lnorm, double d0,
    double d0_search, double score_d8, double t[3], double u[3][3],
    const int mol_type)
{
    D0_MIN = 0.5;
    Lnorm = ylen; // !!! WHY? this doesn't make sense to me
    /* !!! here's a check for RNA molecule type*/
    if (mol_type>0) // RNA
    {
        if     (Lnorm<=11) d0=0.3; 
        else if(Lnorm>11 && Lnorm<=15) d0=0.4;
        else if(Lnorm>15 && Lnorm<=19) d0=0.5;
        else if(Lnorm>19 && Lnorm<=23) d0=0.6;
        else if(Lnorm>23 && Lnorm<30)  d0=0.7;
        else d0=(0.6*pow((Lnorm*1.0-0.5), 1.0/2)-2.5);
    }
    else
    {
        if (Lnorm > 21) d0=(1.24*pow((Lnorm*1.0-15), 1.0/3) -1.8);
        else d0 = D0_MIN;
        if (d0 < D0_MIN) d0 = D0_MIN;
    }
    double d0_input = d0;// Scaled by seq_min

    double tmscore;// collected aligned residues from invmap
    int n_al = 0;
    int i;
    for (int j = 0; j<ylen; j++)
    {
        i = invmap[j];
        if (i >= 0)
        {
            xtm[n_al][0] = x[i][0];
            xtm[n_al][1] = x[i][1];
            xtm[n_al][2] = x[i][2];

            ytm[n_al][0] = y[j][0];
            ytm[n_al][1] = y[j][1];
            ytm[n_al][2] = y[j][2];

            r1[n_al][0] = x[i][0];
            r1[n_al][1] = x[i][1];
            r1[n_al][2] = x[i][2];

            r2[n_al][0] = y[j][0];
            r2[n_al][1] = y[j][1];
            r2[n_al][2] = y[j][2];

            n_al++;
        }

        /* !!! ERROR HANDLING */
        else if (i != -1) PrintErrorAndQuit("Wrong map!\n");

    }
    L_ali = n_al;

    // calculate the optimal rotation between the mobile and target
    Kabsch(r1, r2, n_al, 0, &RMSD, t, u);
    // finish calculating the RMSD
    RMSD = sqrt( RMSD/(1.0*n_al) );
    
    // initiate variables for the scoring method
    int temp_simplify_step = 1;
    int temp_score_sum_method = 0;
    d0_search = d0_input;
    double rms = 0.0;
    // run a TMscore calculation
    tmscore = TMscore8_search_standard(r1, r2, xtm, ytm, xt, n_al, t, u,
        temp_simplify_step, temp_score_sum_method, &rms, d0_input,
        score_d8, d0);
    // finish calculating the TMscore value
    tmscore = tmscore * n_al / (1.0*Lnorm);

    return tmscore;
}

/* Entry function for TM-align. Return TM-score calculation status:
 * 0   - full TM-score calculation 
 * 1   - terminated due to exception
 * 2-7 - pre-terminated due to low TM-score */
int TMalign_main(double **xa, double **ya,
    const char *seqx, const char *seqy, const char *secx, const char *secy,
    double t0[3], double u0[3][3],
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
    if (Lnorm <= 40) ddcc=0.1;   //Lnorm was set in parameter_set4search
    double local_d0_search = d0_search;

/* !!! so many boolean checks. Need to determine which 
    parameters are actually useful*/


    //************************************************//
    //    get initial alignment from user's input:    //
    //    Stick to the initial alignment              //
    //************************************************//
    /* THIS CHUNK OF CODE IS ONLY USED IF THE USER HAS READ IN A FASTA FILE 
       OF A SEQUENCE ALIGNMENT. REMOVING THIS OPTION FOR THE PYTHON BINDINGS
    
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
    } */

    /******************************************************/
    /*   !!! get initial alignment with gapless threading    */
    /******************************************************/
    if (i_opt<=1) // i_opt will always equal 0 if -i or -I are never input; basically, no initial alignment is provided
    {
        get_initial(r1, r2, xtm, ytm, xa, ya, xlen, ylen, invmap0, d0,
            d0_search, fast_opt, t, u);
        TM = detailed_search(r1, r2, xtm, ytm, xt, xa, ya, xlen, ylen, invmap0,
            t, u, simplify_step, score_sum_method, local_d0_search, Lnorm,
            score_d8, d0);
        if (TM>TMmax) TMmax = TM;
        if (TMcut>0) copy_t_u(t, u, t0, u0);
        //run dynamic programing iteratively to find the best alignment
        TM = DP_iter(r1, r2, xtm, ytm, xt, path, val, xa, ya, xlen, ylen,
             t, u, invmap, 0, 2, (fast_opt)?2:30, local_d0_search,
             D0_MIN, Lnorm, d0, score_d8);
        if (TM>TMmax)
        {
            TMmax = TM;
            for (int i = 0; i<ylen; i++) invmap0[i] = invmap[i];
            if (TMcut>0) copy_t_u(t, u, t0, u0);
        }

        // !!! calls a clean up function
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
        get_initial_ss(path, val, secx, secy, xlen, ylen, invmap);
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
            TM = DP_iter(r1, r2, xtm, ytm, xt, path, val, xa, ya,
                xlen, ylen, t, u, invmap, 0, 2, (fast_opt)?2:30,
                local_d0_search, D0_MIN, Lnorm, d0, score_d8);
            if (TM>TMmax)
            {
                TMmax = TM;
                for (int i = 0; i<ylen; i++) invmap0[i] = invmap[i];
                if (TMcut>0) copy_t_u(t, u, t0, u0);
            }
        }

        // !!! calls a clean up function
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
        if (get_initial5( r1, r2, xtm, ytm, path, val, xa, ya,
            xlen, ylen, invmap, d0, d0_search, fast_opt, D0_MIN))
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
                TM = DP_iter(r1, r2, xtm, ytm, xt, path, val, xa, ya,
                    xlen, ylen, t, u, invmap, 0, 2, 2, local_d0_search,
                    D0_MIN, Lnorm, d0, score_d8);
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

        // !!! calls a clean up function
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
        get_initial_ssplus(r1, r2, score, path, val, secx, secy, xa, ya,
            xlen, ylen, invmap0, invmap, D0_MIN, d0);
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
            TM = DP_iter(r1, r2, xtm, ytm, xt, path, val, xa, ya,
                xlen, ylen, t, u, invmap, 0, 2, (fast_opt)?2:30,
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
            TM = DP_iter(r1, r2, xtm, ytm, xt, path, val, xa, ya,
                xlen, ylen, t, u, invmap, 1, 2, 2, local_d0_search, D0_MIN,
                Lnorm, d0, score_d8);
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
    }

    //************************************************//
    //    get initial alignment from user's input:    //
    //************************************************//
    /* !!! */
    if (i_opt>=1 && i_opt<=2)// if input has set parameter for "-i"
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
        TM = DP_iter(r1, r2, xtm, ytm, xt, path, val, xa, ya,
            xlen, ylen, t, u, invmap, 0, 2, (fast_opt)?2:30,
            local_d0_search, D0_MIN, Lnorm, d0, score_d8);
        if (TM>TMmax)
        {
            TMmax = TM;
            for (i = 0; i<ylen; i++) invmap0[i] = invmap[i];
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

    double Lnorm_d0; // !!! ???

    /* !!! normalizing by different values than protein lengths */
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
    Liden=0;
    //double SO=0;
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
        //SO+=(d<3.5);
        kk++;  
        i_old=m1[k]+1;
        j_old=m2[k]+1;
    }
    //SO/=getmin(xlen,ylen);
    //cout<<n_ali8<<'\t'
        //<<rmsd0<<'\t'
        //<<100.*SO<<endl;


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

/* entry function for TM-align with circular permutation
 * i_opt, a_opt, u_opt, d_opt, TMcut are not implemented yet */
int CPalign_main(double **xa, double **ya,
    const char *seqx, const char *seqy, const char *secx, const char *secy,
    double t0[3], double u0[3][3],
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
    const int mol_type, const double TMcut=-1)
{
    char   *seqx_cp; // for the protein sequence 
    char   *secx_cp; // for the secondary structure 
    double **xa_cp;   // coordinates
    string seqxA_cp,seqyA_cp;  // alignment
    int    i,r;
    int    cp_point=0;    // position of circular permutation
    int    cp_aln_best=0; // amount of aligned residue in sliding window
    int    cp_aln_current;// amount of aligned residue in sliding window

    /* duplicate structure, double the seq and sec representation of x 
       since the TMalign_main function only matches with residues down
       stream. */
    NewArray(&xa_cp, xlen*2, 3);
    seqx_cp = new char[xlen*2 + 1];
    secx_cp = new char[xlen*2 + 1];
    for (r=0;r<xlen;r++)
    {
        xa_cp[r+xlen][0]=xa_cp[r][0]=xa[r][0];
        xa_cp[r+xlen][1]=xa_cp[r][1]=xa[r][1];
        xa_cp[r+xlen][2]=xa_cp[r][2]=xa[r][2];
        seqx_cp[r+xlen]=seqx_cp[r]=seqx[r];
        secx_cp[r+xlen]=secx_cp[r]=secx[r];
    }
    seqx_cp[2*xlen]=0;
    secx_cp[2*xlen]=0;
    
    /* fTM-align alignment */
    double TM1_cp,TM2_cp,TM4_cp;
    const double Lnorm_tmp=getmin(xlen,ylen);
    TMalign_main(xa_cp, ya, seqx_cp, seqy, secx_cp, secy,
        t0, u0, TM1_cp, TM2_cp, TM3, TM4_cp, TM5,
        d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out, seqM, seqxA_cp, seqyA_cp,
        rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
        xlen*2, ylen, sequence, Lnorm_tmp, d0_scale,
        0, false, true, false, true, mol_type, -1);

    /* delete gap in seqxA_cp */
    r=0;
    seqxA=seqxA_cp;
    seqyA=seqyA_cp;
    for (i=0;i<seqxA_cp.size();i++)
    {
        if (seqxA_cp[i]!='-')
        {
            seqxA[r]=seqxA_cp[i];
            seqyA[r]=seqyA_cp[i];
            r++;
        }
    }
    seqxA=seqxA.substr(0,r);
    seqyA=seqyA.substr(0,r);

    /* count the number of aligned residues in each window
     * r - residue index in the original unaligned sequence 
     * i - position in the alignment */
    for (r=0;r<xlen-1;r++)
    {
        cp_aln_current=0;
        for (i=r;i<r+xlen;i++) cp_aln_current+=(seqyA[i]!='-');

        if (cp_aln_current>cp_aln_best)
        {
            cp_aln_best=cp_aln_current;
            cp_point=r;
        }
    }
    seqM.clear();
    seqxA.clear();
    seqyA.clear();
    seqxA_cp.clear();
    seqyA_cp.clear();
    rmsd0=Liden=n_ali=n_ali8=0;

    /* fTM-align alignment */
    TMalign_main(xa, ya, seqx, seqy, secx, secy,
        t0, u0, TM1, TM2, TM3, TM4, TM5,
        d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out, seqM, seqxA, seqyA,
        rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
        xlen, ylen, sequence, Lnorm_tmp, d0_scale,
        0, false, true, false, true, mol_type, -1);

    /* do not use circular permutation of number of aligned residues is not
     * larger than sequence-order dependent alignment */
    //cout<<"cp: aln="<<cp_aln_best<<"\tTM="<<TM4_cp<<endl;
    //cout<<"TM: aln="<<n_ali8<<"\tTM="<<TM4<<endl;
    if (n_ali8>=cp_aln_best || TM4>=TM4_cp) cp_point=0;

    /* prepare structure for final alignment */
    seqM.clear();
    seqxA.clear();
    seqyA.clear();
    rmsd0=Liden=n_ali=n_ali8=0;
    if (cp_point!=0)
    {
        for (r=0;r<xlen;r++)
        {
            xa_cp[r][0]=xa_cp[r+cp_point][0];
            xa_cp[r][1]=xa_cp[r+cp_point][1];
            xa_cp[r][2]=xa_cp[r+cp_point][2];
            seqx_cp[r]=seqx_cp[r+cp_point];
            secx_cp[r]=secx_cp[r+cp_point];
        }
    }
    seqx_cp[xlen]=0;
    secx_cp[xlen]=0;

    /* test another round of alignment as concatenated alignment can
     * inflate the number of aligned residues and TM-score. e.g. 1yadA 2duaA */
    if (cp_point!=0)
    {
        TMalign_main(xa_cp, ya, seqx_cp, seqy, secx_cp, secy,
            t0, u0, TM1_cp, TM2_cp, TM3, TM4_cp, TM5,
            d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out, seqM, seqxA_cp, seqyA_cp,
            rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, cp_aln_best,
            xlen, ylen, sequence, Lnorm_tmp, d0_scale,
            0, false, true, false, true, mol_type, -1);
        //cout<<"cp: aln="<<cp_aln_best<<"\tTM="<<TM4_cp<<endl;
        if (n_ali8>=cp_aln_best || TM4>=TM4_cp)
        {
            cp_point=0;
            for (r=0;r<xlen;r++)
            {
                xa_cp[r][0]=xa[r][0];
                xa_cp[r][1]=xa[r][1];
                xa_cp[r][2]=xa[r][2];
                seqx_cp[r]=seqx[r];
                secx_cp[r]=secx[r];
            }
        }
    }

    /* full TM-align */
    TMalign_main(xa_cp, ya, seqx_cp, seqy, secx_cp, secy,
        t0, u0, TM1, TM2, TM3, TM4, TM5,
        d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out, seqM, seqxA_cp, seqyA_cp,
        rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
        xlen, ylen, sequence, Lnorm_ass, d0_scale,
        i_opt, a_opt, u_opt, d_opt, fast_opt, mol_type, TMcut);

    /* correct alignment
     * r - residue index in the original unaligned sequence 
     * i - position in the alignment */
    if (cp_point>0)
    {
        r=0;
        for (i=0;i<seqxA_cp.size();i++)
        {
            r+=(seqxA_cp[i]!='-');
            if (r>=(xlen-cp_point)) 
            {
                i++;
                break;
            }
        }
        seqxA=seqxA_cp.substr(0,i)+'*'+seqxA_cp.substr(i);
        seqM =seqM.substr(0,i)    +' '+seqM.substr(i);
        seqyA=seqyA_cp.substr(0,i)+'-'+seqyA_cp.substr(i);
    }
    else
    {
        seqxA=seqxA_cp;
        seqyA=seqyA_cp;
    }

    /* clean up */
    delete[]seqx_cp;
    delete[]secx_cp;
    DeleteArray(&xa_cp,xlen*2);
    seqxA_cp.clear();
    seqyA_cp.clear();
    return cp_point;
}


// NEED TO MAKE PYTHON BINDINGS AROUND THE FUNCTION TO CALC THE SECONDARY STRUCTURE FUNCTION


/* NEED TO MAKE PYTHON BINDINGS AROUND THIS FUNCTION */
int SOIalign_main(double **xa, double **ya, // prepare cartesian coordinates
    double **xk, double **yk, // prepare the nearest-neighbors array
    const int closeK_opt, // tell the code how many nearest-neigbors are considered
    const char *seqx, const char *seqy, // prepare the sequence string
    const char *secx, const char *secy, // prepare the 2ndary structure string
    // variables to be filled:
    double t0[3], double u0[3][3], // final translation and rotation matrices
    double &TM1, double &TM2, // two relevant TMscores, I'm assuming TM1 is mobile and TM2 is target
    //double &TM3, double &TM4, double &TM5, // various containers for TMscore values; likely unused since a_opt,d_opt,u_opt aren't used here
    double &d0_0, double &TM_0, // final d_{0}, TMscore value BUT TM_0 is not used/output 
    double &d0A, double &d0B, // did we switch from x,y or 1,2 to B,A??? YES YES WE DID
    //double &d0u, double &d0a, // only used if a_opt, u_opt
    //double d0_out, // deleted and recreated, probably messed this up;
    string &seqM, string &seqxA, string &seqyA, // sequence map and aligned sequences; very confusing naming... grumble grumble grumble; A stands for aligned here
    int *invmap, // mapping of target to mobile; defined as `int *invmap = new int[ylen+1]` in `SOIalign` in USalign.cpp; not used external of this function...
    double &rmsd0, // final rmsd value
    //int *L_ali, double *Liden; // deleted and recreated, probably messed this up;
    //double &TM_ali, double &rmsd_ali, int &n_ali, //unused outside of SOIalign_main
    int &n_ali8, // final number of aligned residues... not sure why the 8 is there now...
    const int xlen, const int ylen, // !!! prepare these values; maybe move these up so they sit with the other prepared variables
    //const vector<string> sequence, // !!! only used in CPalign_main call. Need to investigate; might be associated with i_opt
    //const double Lnorm_ass, const double d0_scale, // only used if d_opt
    //const int i_opt, 
    //const int a_opt,
    //const bool u_opt, 
    //const bool d_opt, 
    //const bool fast_opt,
    //const int mol_type, 
    double *dist_list, // vector of distance values... not sure why it's shape is set as ylen+1 in `SOIalign` 
    int **secx_bond, int **secy_bond,  // not const
    const int mm_opt)
{
    /* defining local variables that were originally arg inputs to 
     * `SOIalign_main` pulled out to hardcode their values and/or not pass them
     * upon return */
    double TM3; // only used if a_opt or u_opt, can't remember
    double TM4; // only used if a_opt or u_opt, can't remember
    double TM5; // only used if a_opt or u_opt, can't remember
    double d0u; // only used if u_opt
    double d0a; // only used if a_opt
    
    double d0_out {5.0};// was hard coded as 5.0 /AA in `SOIalign` in 
			// USalign.cpp; distance used to mark residues as 
			// aligned or not in the seqM string
    int L_ali; 		//declared in `SOIalign`
    double Liden {0}; 	// percent aligned? unused 
    double TM_ali;	// ...
    double rmsd_ali;	// ...
    int n_ali {0}; 	// hardcoded in `SOIalign` 
    
    vector<string> sequence; // only used if i_opt
    
    double Lnorm_ass;	// ...
    double d0_scale;	// ...
    // switch variables 
    int i_opt     {0};     // don't do any initial alignment stuff
    int a_opt     {0};     // don't use average of seq lengths to normalize 
			   // TMscores
    /* these variables are using different object types than other _opt 
    * variables */
    bool u_opt    {false}; // don't use user-defined Lnorm_ass value to 
			   // normalize TMscores
    bool d_opt    {false}; // don't use user-defined d_{0} when calculating 
			   // pair weights for TMscore calculations
    bool fast_opt {false}; // don't perform a fast but inaccurate alignment
    int mol_type  {-2};    // for aln btw 2 proteins, mol_type = -2; seems like
			   //  most tests for this are `if >0` then RNA scoring

    // original local variable declarations
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
    int minlen = (xlen<ylen)?xlen:ylen; // min(xlen, ylen);
    int maxlen = (xlen>ylen)?xlen:ylen; // max(xlen, ylen);
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
    parameter_set4search(xlen, ylen, D0_MIN, Lnorm, score_d8, d0, d0_search, dcu0);
    
    // !!! these are being defined well before they are actually used
    // !!! they are actually redefined before they are used; commenting out this instance
    //int simplify_step    = 40; //for simplified search engine
    //int score_sum_method = 8;  //for scoring method, whether only sum over pairs with dis<score_d8

    // prep the inverse and foward maps for residue-mapping
    int i,j;
    int *fwdmap0         = new int[xlen+1]; // mobile to target map
    int *invmap0         = new int[ylen+1]; // target to mobile map
    // set all initial values to -1
    for(i=0; i<xlen; i++) fwdmap0[i]=-1;
    for(j=0; j<ylen; j++) invmap0[j]=-1;
    
    double TMmax=-1, TM=-1;
    double local_d0_search = d0_search;
    
    // !!! maybe iteration_max should be a user-defined value instead of this
    int iteration_max=(fast_opt)?2:30; // fast_opt != 0, iteration_max = 2; fast_opt == 0, iteration_max = 30
    //if (mm_opt==6) iteration_max=1;

    /*************************************************************/
    /* 'initial' alignment with sequence order dependent alignment */
    /*************************************************************/
    CPalign_main(
        xa, ya, seqx, seqy, secx, secy,
        t0, u0, TM1, TM2, TM3, TM4, TM5,
        d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out, 
        seqM, seqxA, seqyA,
        rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
        xlen, ylen, sequence, Lnorm_ass, d0_scale,
        i_opt, a_opt, u_opt, d_opt, fast_opt,
        mol_type,-1);
    
    // sNS alignment-only
    if (mm_opt==6)
    {
        i=0; // mobile seq position counter
        j=0; // target seq position counter
        // loop over length of mobile's aligned sequence space 
        for (int r=0;r<seqxA.size();r++)
        {
            // check for circular permutation point at mobile's seq pos `r`
            if (seqxA[r]=='*') // circular permutation point
            {
                // loop over ...
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
        // update the fwdmap0 variable to mirror the invmap0 values
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

    /*if (a_opt>0)
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
    }*/

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


