#include "TMalign.h"

int score_fun8( double **xa, double **ya, int n_ali, double d, int i_ali[],
    double *score1, int score_sum_method, const double Lnorm, 
    const double score_d8, const double d0,
    double GDT_list_tmp[5], double &maxsub_tmp)
{
    double score_sum=0, di;
    double d_tmp=d*d;
    double d02=d0*d0;
    double score_d8_cut = score_d8*score_d8;
    
    int i, n_cut, inc=0;

    while(1)
    {
        for (i=0;i<5;i++) GDT_list_tmp[i]=0;
        maxsub_tmp=0;

        n_cut=0;
        score_sum=0;
        for(i=0; i<n_ali; i++)
        {
            di = dist(xa[i], ya[i]);
            if(di<d_tmp)
            {
                i_ali[n_cut]=i;
                n_cut++;
            }
            if(score_sum_method==8)
            {                
                if(di<=score_d8_cut) score_sum += 1/(1+di/d02);
            }
            else score_sum += 1/(1+di/d02);

            /* for maxsub score */
            //maxsub_tmp+=1/(1+di/12.25);
            if (di<64) // 8*8=64
            {
                GDT_list_tmp[4]+=1;
                if (di<16) // 4*4=16
                {
                    GDT_list_tmp[3]+=1;
                    if (di<12.25) // 3.5^2=12.25
                    {
                        maxsub_tmp+=1/(1+di/12.25);
                        if (di<4) // 2*2=4
                        {
                            GDT_list_tmp[2]+=1;
                            if (di<1) // 1*1=1
                            {
                                GDT_list_tmp[1]+=1;
                                if (di<0.25) // 0.5*0.5=0.25
                                    GDT_list_tmp[0]+=1;
                            }
                        }
                    }
                }
            }
        }
        //there are not enough feasible pairs, relieve the threshold         
        if(n_cut<3 && n_ali>3)
        {
            inc++;
            double dinc=(d+inc*0.5);
            d_tmp = dinc * dinc;
        }
        else break;
    }  

    *score1=score_sum/Lnorm;
    return n_cut;
}

int score_fun8_standard(double **xa, double **ya, int n_ali, double d,
    int i_ali[], double *score1, int score_sum_method,
    double score_d8, double d0, double GDT_list_tmp[5], double &maxsub_tmp)
{
    double score_sum = 0, di;
    double d_tmp = d*d;
    double d02 = d0*d0;
    double score_d8_cut = score_d8*score_d8;

    int i, n_cut, inc = 0;
    while (1)
    {
        for (i=0;i<5;i++) GDT_list_tmp[i]=0;
        maxsub_tmp=0;
        n_cut = 0;
        score_sum = 0;
        for (i = 0; i<n_ali; i++)
        {
            di = dist(xa[i], ya[i]);
            if (di<d_tmp)
            {
                i_ali[n_cut] = i;
                n_cut++;
            }
            if (score_sum_method == 8)
            {
                if (di <= score_d8_cut) score_sum += 1 / (1 + di / d02);
            }
            else
            {
                score_sum += 1 / (1 + di / d02);
            }

            /* for maxsub score */
            //maxsub_tmp+=1/(1+di/12.25);
            if (di<64) // 8*8=64
            {
                GDT_list_tmp[4]+=1;
                if (di<16) // 4*4=16
                {
                    GDT_list_tmp[3]+=1;
                    if (di<12.25) // 3.5^2=12.25
                    {
                        maxsub_tmp+=1/(1+di/12.25);
                        if (di<4) // 2*2=4
                        {
                            GDT_list_tmp[2]+=1;
                            if (di<1) // 1*1=1
                            {
                                GDT_list_tmp[1]+=1;
                                if (di<0.25) // 0.5*0.5=0.25
                                    GDT_list_tmp[0]+=1;
                            }
                        }
                    }
                }
            }
        }
        //there are not enough feasible pairs, relieve the threshold         
        if (n_cut<3 && n_ali>3)
        {
            inc++;
            double dinc = (d + inc*0.5);
            d_tmp = dinc * dinc;
        }
        else break;
    }

    *score1 = score_sum / n_ali;
    return n_cut;
}

double TMscore8_search(double **r1, double **r2, double **xtm, double **ytm,
    double **xt, int Lali, double t0[3], double u0[3][3], int simplify_step,
    int score_sum_method, double *Rcomm, double local_d0_search, double Lnorm,
    double score_d8, double d0, double GDT_list[5], double &maxsub)
{
    double GDT_list_tmp[5]={0,0,0,0,0};
    double maxsub_tmp=0;
    int i, m;
    double score_max, score, rmsd;    
    const int kmax=Lali;    
    int k_ali[kmax], ka, k;
    double t[3];
    double u[3][3];
    double d;
    

    //iterative parameters
    int n_it=20;            //maximum number of iterations
    int n_init_max=6; //maximum number of different fragment length 
    int L_ini[n_init_max];  //fragment lengths, Lali, Lali/2, Lali/4 ... 4   
    int L_ini_min=4;
    if(Lali<L_ini_min) L_ini_min=Lali;   

    int n_init=0, i_init;      
    for(i=0; i<n_init_max-1; i++)
    {
        n_init++;
        L_ini[i]=(int) (Lali/pow(2.0, (double) i));
        if(L_ini[i]<=L_ini_min)
        {
            L_ini[i]=L_ini_min;
            break;
        }
    }
    if(i==n_init_max-1)
    {
        n_init++;
        L_ini[i]=L_ini_min;
    }
    
    score_max=-1;
    //find the maximum score starting from local structures superposition
    int i_ali[kmax], n_cut;
    int L_frag; //fragment length
    int iL_max; //maximum starting postion for the fragment
    
    for(i_init=0; i_init<n_init; i_init++)
    {
        L_frag=L_ini[i_init];
        iL_max=Lali-L_frag;
      
        i=0;   
        while(1)
        {
            //extract the fragment starting from position i 
            ka=0;
            for(k=0; k<L_frag; k++)
            {
                int kk=k+i;
                r1[k][0]=xtm[kk][0];  
                r1[k][1]=xtm[kk][1]; 
                r1[k][2]=xtm[kk][2];   
                
                r2[k][0]=ytm[kk][0];  
                r2[k][1]=ytm[kk][1]; 
                r2[k][2]=ytm[kk][2];
                
                k_ali[ka]=kk;
                ka++;
            }
            
            //extract rotation matrix based on the fragment
            Kabsch(r1, r2, L_frag, 1, &rmsd, t, u);
            if (simplify_step != 1)
                *Rcomm = 0;
            do_rotation(xtm, xt, Lali, t, u);
            
            //get subsegment of this fragment
            d = local_d0_search - 1;
            n_cut=score_fun8(xt, ytm, Lali, d, i_ali, &score, 
                score_sum_method, Lnorm, score_d8, d0, 
                GDT_list_tmp, maxsub_tmp);
            if(score>score_max)
            {
                score_max=score;
                
                //save the rotation matrix
                for(k=0; k<3; k++)
                {
                    t0[k]=t[k];
                    u0[k][0]=u[k][0];
                    u0[k][1]=u[k][1];
                    u0[k][2]=u[k][2];
                }
            }
            if (maxsub_tmp>maxsub) maxsub=maxsub_tmp;
            for (k=0;k<5;k++)
                if (GDT_list_tmp[k]>GDT_list[k])
                    GDT_list[k]=GDT_list_tmp[k];
            
            //try to extend the alignment iteratively            
            d = local_d0_search + 1;
            for(int it=0; it<n_it; it++)            
            {
                ka=0;
                for(k=0; k<n_cut; k++)
                {
                    m=i_ali[k];
                    r1[k][0]=xtm[m][0];  
                    r1[k][1]=xtm[m][1]; 
                    r1[k][2]=xtm[m][2];
                    
                    r2[k][0]=ytm[m][0];  
                    r2[k][1]=ytm[m][1]; 
                    r2[k][2]=ytm[m][2];
                    
                    k_ali[ka]=m;
                    ka++;
                } 
                //extract rotation matrix based on the fragment                
                Kabsch(r1, r2, n_cut, 1, &rmsd, t, u);
                do_rotation(xtm, xt, Lali, t, u);
                n_cut=score_fun8(xt, ytm, Lali, d, i_ali, &score, 
                    score_sum_method, Lnorm, score_d8, d0);
                if(score>score_max)
                {
                    score_max=score;

                    //save the rotation matrix
                    for(k=0; k<3; k++)
                    {
                        t0[k]=t[k];
                        u0[k][0]=u[k][0];
                        u0[k][1]=u[k][1];
                        u0[k][2]=u[k][2];
                    }                     
                }
                if (maxsub_tmp>maxsub) maxsub=maxsub_tmp;
                for (k=0;k<5;k++)
                    if (GDT_list_tmp[k]>GDT_list[k])
                        GDT_list[k]=GDT_list_tmp[k];
                
                //check if it converges            
                if(n_cut==ka)
                {                
                    for(k=0; k<n_cut; k++)
                    {
                        if(i_ali[k]!=k_ali[k]) break;
                    }
                    if(k==n_cut) break;
                }                                                               
            } //for iteration            

            if(i<iL_max)
            {
                i=i+simplify_step; //shift the fragment        
                if(i>iL_max) i=iL_max;  //do this to use the last missed fragment
            }
            else if(i>=iL_max) break;
        }//while(1)
        //end of one fragment
    }//for(i_init
    return score_max;
}


double TMscore8_search_standard( double **r1, double **r2,
    double **xtm, double **ytm, double **xt, int Lali,
    double t0[3], double u0[3][3], int simplify_step, int score_sum_method,
    double *Rcomm, double local_d0_search, double score_d8, double d0,
    double GDT_list[5], double &maxsub)
{
    double GDT_list_tmp[5]={0,0,0,0,0};
    double maxsub_tmp=0;
    int i, m;
    double score_max, score, rmsd;
    const int kmax = Lali;
    int k_ali[kmax], ka, k;
    double t[3];
    double u[3][3];
    double d;

    //iterative parameters
    int n_it = 20;            //maximum number of iterations
    int n_init_max = 6; //maximum number of different fragment length 
    int L_ini[n_init_max];  //fragment lengths, Lali, Lali/2, Lali/4 ... 4   
    int L_ini_min = 4;
    if (Lali<L_ini_min) L_ini_min = Lali;

    int n_init = 0, i_init;
    for (i = 0; i<n_init_max - 1; i++)
    {
        n_init++;
        L_ini[i] = (int)(Lali / pow(2.0, (double)i));
        if (L_ini[i] <= L_ini_min)
        {
            L_ini[i] = L_ini_min;
            break;
        }
    }
    if (i == n_init_max - 1)
    {
        n_init++;
        L_ini[i] = L_ini_min;
    }

    score_max = -1;
    //find the maximum score starting from local structures superposition
    int i_ali[kmax], n_cut;
    int L_frag; //fragment length
    int iL_max; //maximum starting position for the fragment

    for (i_init = 0; i_init<n_init; i_init++)
    {
        L_frag = L_ini[i_init];
        iL_max = Lali - L_frag;

        i = 0;
        while (1)
        {
            //extract the fragment starting from position i 
            ka = 0;
            for (k = 0; k<L_frag; k++)
            {
                int kk = k + i;
                r1[k][0] = xtm[kk][0];
                r1[k][1] = xtm[kk][1];
                r1[k][2] = xtm[kk][2];

                r2[k][0] = ytm[kk][0];
                r2[k][1] = ytm[kk][1];
                r2[k][2] = ytm[kk][2];

                k_ali[ka] = kk;
                ka++;
            }
            //extract rotation matrix based on the fragment
            Kabsch(r1, r2, L_frag, 1, &rmsd, t, u);
            if (simplify_step != 1)
                *Rcomm = 0;
            do_rotation(xtm, xt, Lali, t, u);

            //get subsegment of this fragment
            d = local_d0_search - 1;
            n_cut = score_fun8_standard(xt, ytm, Lali, d, i_ali, &score,
                score_sum_method, score_d8, d0, GDT_list_tmp, maxsub_tmp);

            if (score>score_max)
            {
                score_max = score;

                //save the rotation matrix
                for (k = 0; k<3; k++)
                {
                    t0[k] = t[k];
                    u0[k][0] = u[k][0];
                    u0[k][1] = u[k][1];
                    u0[k][2] = u[k][2];
                }
            }
            if (maxsub_tmp>maxsub) maxsub=maxsub_tmp;
            for (k=0;k<5;k++)
                if (GDT_list_tmp[k]>GDT_list[k])
                    GDT_list[k]=GDT_list_tmp[k];

            //try to extend the alignment iteratively            
            d = local_d0_search + 1;
            for (int it = 0; it<n_it; it++)
            {
                ka = 0;
                for (k = 0; k<n_cut; k++)
                {
                    m = i_ali[k];
                    r1[k][0] = xtm[m][0];
                    r1[k][1] = xtm[m][1];
                    r1[k][2] = xtm[m][2];

                    r2[k][0] = ytm[m][0];
                    r2[k][1] = ytm[m][1];
                    r2[k][2] = ytm[m][2];

                    k_ali[ka] = m;
                    ka++;
                }
                //extract rotation matrix based on the fragment                
                Kabsch(r1, r2, n_cut, 1, &rmsd, t, u);
                do_rotation(xtm, xt, Lali, t, u);
                n_cut = score_fun8_standard(xt, ytm, Lali, d, i_ali, &score,
                    score_sum_method, score_d8, d0, GDT_list_tmp, maxsub_tmp);
                if (score>score_max)
                {
                    score_max = score;

                    //save the rotation matrix
                    for (k = 0; k<3; k++)
                    {
                        t0[k] = t[k];
                        u0[k][0] = u[k][0];
                        u0[k][1] = u[k][1];
                        u0[k][2] = u[k][2];
                    }
                }
                if (maxsub_tmp>maxsub) maxsub=maxsub_tmp;
                for (k=0;k<5;k++)
                    if (GDT_list_tmp[k]>GDT_list[k])
                        GDT_list[k]=GDT_list_tmp[k];

                //check if it converges            
                if (n_cut == ka)
                {
                    for (k = 0; k<n_cut; k++)
                    {
                        if (i_ali[k] != k_ali[k]) break;
                    }
                    if (k == n_cut) break;
                }
            } //for iteration            

            if (i<iL_max)
            {
                i = i + simplify_step; //shift the fragment        
                if (i>iL_max) i = iL_max;  //do this to use the last missed fragment
            }
            else if (i >= iL_max) break;
        }//while(1)
        //end of one fragment
    }//for(i_init
    return score_max;
}

double detailed_search_standard( double **r1, double **r2,
    double **xtm, double **ytm, double **xt, double **x, double **y,
    int xlen, int ylen, int invmap0[], double t[3], double u[3][3],
    int simplify_step, int score_sum_method, double local_d0_search,
    const bool& bNormalize, double Lnorm, double score_d8, double d0,
    double GDT_list[5], double &maxsub)
{
    //x is model, y is template, try to superpose onto y
    int i, j, k;     
    double tmscore;
    double rmsd;

    k=0;
    for(i=0; i<ylen; i++) 
    {
        j=invmap0[i];
        if(j>=0) //aligned
        {
            xtm[k][0]=x[j][0];
            xtm[k][1]=x[j][1];
            xtm[k][2]=x[j][2];
                
            ytm[k][0]=y[i][0];
            ytm[k][1]=y[i][1];
            ytm[k][2]=y[i][2];
            k++;
        }
    }

    //detailed search 40-->1
    tmscore = TMscore8_search_standard( r1, r2, xtm, ytm, xt, k, t, u,
        simplify_step, score_sum_method, &rmsd, local_d0_search, score_d8, d0,
        GDT_list, maxsub);
    if (bNormalize)// "-i", to use standard_TMscore, then bNormalize=true, else bNormalize=false; 
        tmscore = tmscore * k / Lnorm;

    return tmscore;
}

/* Entry function for TM-score. Return TM-score calculation status:
 * 0   - full TM-score calculation 
 * 1   - terminated due to exception
 * 2-7 - pre-terminated due to low TM-score */
int TMscore_main(double **xa, double **ya,
    const char *seqx, const char *seqy, double t0[3], double u0[3][3],
    double &TM1, double &TM2, double &TM3, double &TM4, double &TM5,
    double &d0_0, double &TM_0,
    double &d0A, double &d0B, double &d0u, double &d0a, double &d0_out,
    string &seqM, string &seqxA, string &seqyA,
    double &rmsd0, int &L_ali, double &Liden,
    double &TM_ali, double &rmsd_ali, int &n_ali, int &n_ali8,
    const int xlen, const int ylen,
    const vector<string> sequence, const double Lnorm_ass,
    const double d0_scale, const int a_opt,
    const bool u_opt, const bool d_opt, const bool fast_opt,
    const int mol_type, double GDT_list[5], double &maxsub,
    const double TMcut=-1)
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
    //    Stick to the initial alignment              //
    //************************************************//
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
        false, Lnorm, score_d8, d0,
        GDT_list, maxsub);

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
        score_sum_method, &rmsd, local_d0_search, Lnorm, score_d8, d0,
        GDT_list, maxsub);
    TM_0 = TM1;

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
        //if(d<d0_out) seqM[kk]=':';
        //else         seqM[kk]='.';
        if(d<5) seqM[kk]=':';
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

void output_TMscore_results(
    const string xname, const string yname,
    const string chainID1, const string chainID2,
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
    int L_lt_d, const double rmsd_d0_out,
    double GDT_list[5], double maxsub, const int split_opt,
    const vector<string>&resi_vec1, const vector<string>&resi_vec2)
{
    if (outfmt_opt<=0)
    {
        printf("\nStructure1: %s%s    Length=%5d\n",
            xname.c_str(), chainID1.c_str(), xlen);
        printf("Structure2: %s%s    Length=%5d (by which all scores are normalized)\n",
            yname.c_str(), chainID2.c_str(), ylen);

        printf("Number of residues in common=%5d\n", n_ali8);
        printf("RMSD of  the common residues=%9.3f\n\n", rmsd);
        printf("TM-score    = %6.4f  (d0= %.2f)\n", TM1, d0A);
        printf("MaxSub-score= %6.4f  (d0= 3.50)\n", maxsub/ylen);

        double gdt_ts_score=0;
        double gdt_ha_score=0;
        int i;
        for (i=0;i<4;i++)
        {
            gdt_ts_score+=GDT_list[i+1];
            gdt_ha_score+=GDT_list[i];
        }
        gdt_ts_score/=(4*ylen);
        gdt_ha_score/=(4*ylen);
        printf("GDT-TS-score= %6.4f %%(d<1)=%6.4f %%(d<2)=%6.4f %%(d<4)=%6.4f %%(d<8)=%6.4f\n",
            gdt_ts_score, GDT_list[1]/ylen, GDT_list[2]/ylen,
                          GDT_list[3]/ylen, GDT_list[4]/ylen);
        printf("GDT-HA-score= %6.4f %%(d<0.5)=%6.4f %%(d<1)=%6.4f %%(d<2)=%6.4f %%(d<4)=%6.4f\n",
            gdt_ha_score, GDT_list[0]/ylen, GDT_list[1]/ylen,
                          GDT_list[2]/ylen, GDT_list[3]/ylen);

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
        string seq_scale=seqM;
        for (i=0;i<strlen(seqM);i++)
        {
            L_lt_d+=seqM[i]==':';
            seq_scale[i]=(i+1)%10+'0';
        }
        printf("\nSuperposition in the TM-score: Length(d<%3.1f)= %d\n", d0_out, L_lt_d);
        //printf("\nSuperposition in the TM-score: Length(d<%3.1f)= %d  RMSD=%6.2f\n", d0_out, L_lt_d, rmsd_d0_out);
        printf("(\":\" denotes the residue pairs of distance <%4.1f Angstrom)\n", d0_out);
        printf("%s\n", seqxA);
        printf("%s\n", seqM);
        printf("%s\n", seqyA);
        printf("%s\n", seq_scale.c_str());
        seq_scale.clear();
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

    if (strlen(fname_matrix)) 
        output_rotation_matrix(fname_matrix, t, u);
    if (strlen(fname_super))
        output_pymol(xname, yname, fname_super, t, u, ter_opt, 
            0, split_opt, mirror_opt, seqM, seqxA, seqyA,
            resi_vec1, resi_vec2, chainID1, chainID2);
}
