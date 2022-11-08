/* Functions for the core TMalign algorithm, including the entry function
 * TMalign_main */
#ifndef TMalign_h
#define TMalign_h 1

#include "param_set.h"
#include "NW.h"
#include "Kabsch.h"
#include "NWalign.h"

//     1, collect those residues with dis<d;
//     2, calculate TMscore
int score_fun8( double **xa, double **ya, int n_ali, double d, int i_ali[],
    double *score1, int score_sum_method, const double Lnorm, 
    const double score_d8, const double d0)
{
    double score_sum=0, di;
    double d_tmp=d*d;
    double d02=d0*d0;
    double score_d8_cut = score_d8*score_d8;
    
    int i, n_cut, inc=0;

    while(1)
    {
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
    double score_d8, double d0)
{
    double score_sum = 0, di;
    double d_tmp = d*d;
    double d02 = d0*d0;
    double score_d8_cut = score_d8*score_d8;

    int i, n_cut, inc = 0;
    while (1)
    {
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
    double score_d8, double d0)
{
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
    int iL_max; //maximum starting position for the fragment
    
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
    double *Rcomm, double local_d0_search, double score_d8, double d0)
{
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
                score_sum_method, score_d8, d0);

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
                    score_sum_method, score_d8, d0);
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

//Comprehensive TMscore search engine
// input:   two vector sets: x, y
//          an alignment invmap0[] between x and y
//          simplify_step: 1 or 40 or other integers
//          score_sum_method: 0 for score over all pairs
//                            8 for socre over the pairs with dist<score_d8
// output:  the best rotaion matrix t, u that results in highest TMscore
double detailed_search(double **r1, double **r2, double **xtm, double **ytm,
    double **xt, double **x, double **y, int xlen, int ylen, 
    int invmap0[], double t[3], double u[3][3], int simplify_step,
    int score_sum_method, double local_d0_search, double Lnorm,
    double score_d8, double d0)
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
    tmscore = TMscore8_search(r1, r2, xtm, ytm, xt, k, t, u, simplify_step,
        score_sum_method, &rmsd, local_d0_search, Lnorm, score_d8, d0);
    return tmscore;
}

double detailed_search_standard( double **r1, double **r2,
    double **xtm, double **ytm, double **xt, double **x, double **y,
    int xlen, int ylen, int invmap0[], double t[3], double u[3][3],
    int simplify_step, int score_sum_method, double local_d0_search,
    const bool& bNormalize, double Lnorm, double score_d8, double d0)
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
        simplify_step, score_sum_method, &rmsd, local_d0_search, score_d8, d0);
    if (bNormalize)// "-i", to use standard_TMscore, then bNormalize=true, else bNormalize=false; 
        tmscore = tmscore * k / Lnorm;

    return tmscore;
}

//compute the score quickly in three iterations
double get_score_fast( double **r1, double **r2, double **xtm, double **ytm,
    double **x, double **y, int xlen, int ylen, int invmap[],
    double d0, double d0_search, double t[3], double u[3][3])
{
    double rms, tmscore, tmscore1, tmscore2;
    int i, j, k;

    k=0;
    for(j=0; j<ylen; j++)
    {
        i=invmap[j];
        if(i>=0)
        {
            r1[k][0]=x[i][0];
            r1[k][1]=x[i][1];
            r1[k][2]=x[i][2];

            r2[k][0]=y[j][0];
            r2[k][1]=y[j][1];
            r2[k][2]=y[j][2];
            
            xtm[k][0]=x[i][0];
            xtm[k][1]=x[i][1];
            xtm[k][2]=x[i][2];
            
            ytm[k][0]=y[j][0];
            ytm[k][1]=y[j][1];
            ytm[k][2]=y[j][2];                  
            
            k++;
        }
        else if(i!=-1) PrintErrorAndQuit("Wrong map!\n");
    }
    Kabsch(r1, r2, k, 1, &rms, t, u);
    
    //evaluate score   
    double di;
    const int len=k;
    double dis[len];    
    double d00=d0_search;
    double d002=d00*d00;
    double d02=d0*d0;
    
    int n_ali=k;
    double xrot[3];
    tmscore=0;
    for(k=0; k<n_ali; k++)
    {
        transform(t, u, &xtm[k][0], xrot);        
        di=dist(xrot, &ytm[k][0]);
        dis[k]=di;
        tmscore += 1/(1+di/d02);
    }
    
   
   
   //second iteration 
    double d002t=d002;
    vector<double> dis_vec(dis, dis+n_ali);
    sort(dis_vec.begin(), dis_vec.end());
    if (d002t<dis_vec[2]) d002t=dis_vec[2];
    dis_vec.clear();
    while(1)
    {
        j=0;
        for(k=0; k<n_ali; k++)
        {            
            if(dis[k]<=d002t)
            {
                r1[j][0]=xtm[k][0];
                r1[j][1]=xtm[k][1];
                r1[j][2]=xtm[k][2];
                
                r2[j][0]=ytm[k][0];
                r2[j][1]=ytm[k][1];
                r2[j][2]=ytm[k][2];
                
                j++;
            }
        }
        //there are not enough feasible pairs, relieve the threshold 
        if(j<3 && n_ali>3) d002t += 0.5;
        else break;
    }
    
    if(n_ali!=j)
    {
        Kabsch(r1, r2, j, 1, &rms, t, u);
        tmscore1=0;
        for(k=0; k<n_ali; k++)
        {
            transform(t, u, &xtm[k][0], xrot);        
            di=dist(xrot, &ytm[k][0]);
            dis[k]=di;
            tmscore1 += 1/(1+di/d02);
        }
        
        //third iteration
        d002t=d002+1;
        vector<double> dis_vec(dis, dis+n_ali);
        sort(dis_vec.begin(), dis_vec.end());
        if (d002t<dis_vec[2]) d002t=dis_vec[2];
        dis_vec.clear();
        while(1)
        {
            j=0;
            for(k=0; k<n_ali; k++)
            {            
                if(dis[k]<=d002t)
                {
                    r1[j][0]=xtm[k][0];
                    r1[j][1]=xtm[k][1];
                    r1[j][2]=xtm[k][2];
                    
                    r2[j][0]=ytm[k][0];
                    r2[j][1]=ytm[k][1];
                    r2[j][2]=ytm[k][2];
                                        
                    j++;
                }
            }
            //there are not enough feasible pairs, relieve the threshold 
            if(j<3 && n_ali>3) d002t += 0.5;
            else break;
        }

        //evaluate the score
        Kabsch(r1, r2, j, 1, &rms, t, u);
        tmscore2=0;
        for(k=0; k<n_ali; k++)
        {
            transform(t, u, &xtm[k][0], xrot);
            di=dist(xrot, &ytm[k][0]);
            tmscore2 += 1/(1+di/d02);
        }    
    }
    else
    {
        tmscore1=tmscore;
        tmscore2=tmscore;
    }
    
    if(tmscore1>=tmscore) tmscore=tmscore1;
    if(tmscore2>=tmscore) tmscore=tmscore2;
    return tmscore; // no need to normalize this score because it will not be used for latter scoring
}


//perform gapless threading to find the best initial alignment
//input: x, y, xlen, ylen
//output: y2x0 stores the best alignment: e.g., 
//y2x0[j]=i means:
//the jth element in y is aligned to the ith element in x if i>=0 
//the jth element in y is aligned to a gap in x if i==-1
double get_initial(double **r1, double **r2, double **xtm, double **ytm,
    double **x, double **y, int xlen, int ylen, int *y2x,
    double d0, double d0_search, const bool fast_opt,
    double t[3], double u[3][3])
{
    int min_len=getmin(xlen, ylen);
    if(min_len<3) PrintErrorAndQuit("Sequence is too short <3!\n");
    
    int min_ali= min_len/2;              //minimum size of considered fragment 
    if(min_ali<=5)  min_ali=5;    
    int n1, n2;
    n1 = -ylen+min_ali; 
    n2 = xlen-min_ali;

    int i, j, k, k_best;
    double tmscore, tmscore_max=-1;

    k_best=n1;
    for(k=n1; k<=n2; k+=(fast_opt)?5:1)
    {
        //get the map
        for(j=0; j<ylen; j++)
        {
            i=j+k;
            if(i>=0 && i<xlen) y2x[j]=i;
            else y2x[j]=-1;
        }
        
        //evaluate the map quickly in three iterations
        //this is not real tmscore, it is used to evaluate the goodness of the initial alignment
        tmscore=get_score_fast(r1, r2, xtm, ytm,
            x, y, xlen, ylen, y2x, d0,d0_search, t, u);
        if(tmscore>=tmscore_max)
        {
            tmscore_max=tmscore;
            k_best=k;
        }
    }
    
    //extract the best map
    k=k_best;
    for(j=0; j<ylen; j++)
    {
        i=j+k;
        if(i>=0 && i<xlen) y2x[j]=i;
        else y2x[j]=-1;
    }    

    return tmscore_max;
}

void smooth(int *sec, int len)
{
    int i, j;
    //smooth single  --x-- => -----
    for (i=2; i<len-2; i++)
    {
        if(sec[i]==2 || sec[i]==4)
        {
            j=sec[i];
            if (sec[i-2]!=j && sec[i-1]!=j && sec[i+1]!=j && sec[i+2]!=j)
                sec[i]=1;
        }
    }

    //   smooth double 
    //   --xx-- => ------
    for (i=0; i<len-5; i++)
    {
        //helix
        if (sec[i]!=2   && sec[i+1]!=2 && sec[i+2]==2 && sec[i+3]==2 &&
            sec[i+4]!=2 && sec[i+5]!= 2)
        {
            sec[i+2]=1;
            sec[i+3]=1;
        }

        //beta
        if (sec[i]!=4   && sec[i+1]!=4 && sec[i+2]==4 && sec[i+3]==4 &&
            sec[i+4]!=4 && sec[i+5]!= 4)
        {
            sec[i+2]=1;
            sec[i+3]=1;
        }
    }

    //smooth connect
    for (i=0; i<len-2; i++)
    {        
        if (sec[i]==2 && sec[i+1]!=2 && sec[i+2]==2) sec[i+1]=2;
        else if(sec[i]==4 && sec[i+1]!=4 && sec[i+2]==4) sec[i+1]=4;
    }

}

char sec_str(double dis13, double dis14, double dis15,
            double dis24, double dis25, double dis35)
{
    char s='C';
    
    double delta=2.1;
    if (fabs(dis15-6.37)<delta && fabs(dis14-5.18)<delta && 
        fabs(dis25-5.18)<delta && fabs(dis13-5.45)<delta &&
        fabs(dis24-5.45)<delta && fabs(dis35-5.45)<delta)
    {
        s='H'; //helix                        
        return s;
    }

    delta=1.42;
    if (fabs(dis15-13  )<delta && fabs(dis14-10.4)<delta &&
        fabs(dis25-10.4)<delta && fabs(dis13-6.1 )<delta &&
        fabs(dis24-6.1 )<delta && fabs(dis35-6.1 )<delta)
    {
        s='E'; //strand
        return s;
    }

    if (dis15 < 8) s='T'; //turn
    return s;
}


/* secondary structure assignment for protein:
 * 1->coil, 2->helix, 3->turn, 4->strand */
void make_sec(double **x, int len, char *sec)
{
    int j1, j2, j3, j4, j5;
    double d13, d14, d15, d24, d25, d35;
    for(int i=0; i<len; i++)
    {     
        sec[i]='C';
        j1=i-2;
        j2=i-1;
        j3=i;
        j4=i+1;
        j5=i+2;        
        
        if(j1>=0 && j5<len)
        {
            d13=sqrt(dist(x[j1], x[j3]));
            d14=sqrt(dist(x[j1], x[j4]));
            d15=sqrt(dist(x[j1], x[j5]));
            d24=sqrt(dist(x[j2], x[j4]));
            d25=sqrt(dist(x[j2], x[j5]));
            d35=sqrt(dist(x[j3], x[j5]));
            sec[i]=sec_str(d13, d14, d15, d24, d25, d35);            
        }    
    } 
    sec[len]=0;
}

/* a c d b: a paired to b, c paired to d */
bool overlap(const int a1,const int b1,const int c1,const int d1,
             const int a2,const int b2,const int c2,const int d2)
{
    return (a2>=a1&&a2<=c1)||(c2>=a1&&c2<=c1)||
           (d2>=a1&&d2<=c1)||(b2>=a1&&b2<=c1)||
           (a2>=d1&&a2<=b1)||(c2>=d1&&c2<=b1)||
           (d2>=d1&&d2<=b1)||(b2>=d1&&b2<=b1);
}

/* find base pairing stacks in RNA*/
void sec_str(int len,char *seq, const vector<vector<bool> >&bp, 
    int a, int b,int &c, int &d)
{
    int i;
    
    for (i=0;i<len;i++)
    {
        if (a+i<len-3 && b-i>0)
        {
            if (a+i<b-i && bp[a+i][b-i]) continue;
            break;
        }
    }
    c=a+i-1;d=b-i+1;
}

/* secondary structure assignment for RNA:
 * 1->unpair, 2->paired with upstream, 3->paired with downstream */
void make_sec(char *seq, double **x, int len, char *sec,const string atom_opt)
{
    int ii,jj,i,j;

    float lb=12.5; // lower bound for " C3'"
    float ub=15.0; // upper bound for " C3'"
    if     (atom_opt==" C4'") {lb=14.0;ub=16.0;}
    else if(atom_opt==" C5'") {lb=16.0;ub=18.0;}
    else if(atom_opt==" O3'") {lb=13.5;ub=16.5;}
    else if(atom_opt==" O5'") {lb=15.5;ub=18.5;}
    else if(atom_opt==" P  ") {lb=16.5;ub=21.0;}

    float dis;
    vector<bool> bp_tmp(len,false);
    vector<vector<bool> > bp(len,bp_tmp);
    bp_tmp.clear();
    for (i=0; i<len; i++)
    {
        sec[i]='.';
        for (j=i+1; j<len; j++)
        {
            if (((seq[i]=='u'||seq[i]=='t')&&(seq[j]=='a'             ))||
                ((seq[i]=='a'             )&&(seq[j]=='u'||seq[j]=='t'))||
                ((seq[i]=='g'             )&&(seq[j]=='c'||seq[j]=='u'))||
                ((seq[i]=='c'||seq[i]=='u')&&(seq[j]=='g'             )))
            {
                dis=sqrt(dist(x[i], x[j]));
                bp[j][i]=bp[i][j]=(dis>lb && dis<ub);
            }
        }
    }
    
    // From 5' to 3': A0 C0 D0 B0: A0 paired to B0, C0 paired to D0
    vector<int> A0,B0,C0,D0;
    for (i=0; i<len-2; i++)
    {
        for (j=i+3; j<len; j++)
        {
            if (!bp[i][j]) continue;
            if (i>0 && j+1<len && bp[i-1][j+1]) continue;
            if (!bp[i+1][j-1]) continue;
            sec_str(len,seq, bp, i,j,ii,jj);
            if (jj<i || j<ii)
            {
                ii=i;
                jj=j;
            }
            A0.push_back(i);
            B0.push_back(j);
            C0.push_back(ii);
            D0.push_back(jj);
        }
    }
    
    //int sign;
    for (i=0;i<A0.size();i++)
    {
        /*
        sign=0;
        if(C0[i]-A0[i]<=1)
        {
            for(j=0;j<A0.size();j++)
            {
                if(i==j) continue;

                if((A0[j]>=A0[i]&&A0[j]<=C0[i])||
                   (C0[j]>=A0[i]&&C0[j]<=C0[i])||
                   (D0[j]>=A0[i]&&D0[j]<=C0[i])||
                   (B0[j]>=A0[i]&&B0[j]<=C0[i])||
                   (A0[j]>=D0[i]&&A0[j]<=B0[i])||
                   (C0[j]>=D0[i]&&C0[j]<=B0[i])||
                   (D0[j]>=D0[i]&&D0[j]<=B0[i])||
                   (B0[j]>=D0[i]&&B0[j]<=B0[i]))
                {
                    sign=-1;
                    break;
                }
            }
        }
        if(sign!=0) continue;
        */

        for (j=0;;j++)
        {
            if(A0[i]+j>C0[i]) break;
            sec[A0[i]+j]='<';
            sec[D0[i]+j]='>';
        }
    }
    sec[len]=0;

    /* clean up */
    A0.clear();
    B0.clear();
    C0.clear();
    D0.clear();
    bp.clear();
}

//get initial alignment from secondary structure alignment
//input: x, y, xlen, ylen
//output: y2x stores the best alignment: e.g., 
//y2x[j]=i means:
//the jth element in y is aligned to the ith element in x if i>=0 
//the jth element in y is aligned to a gap in x if i==-1
void get_initial_ss(bool **path, double **val,
    const char *secx, const char *secy, int xlen, int ylen, int *y2x)
{
    double gap_open=-1.0;
    NWDP_TM(path, val, secx, secy, xlen, ylen, gap_open, y2x);
}


// get_initial5 in TMalign fortran, get_initial_local in TMalign c by yangji
//get initial alignment of local structure superposition
//input: x, y, xlen, ylen
//output: y2x stores the best alignment: e.g., 
//y2x[j]=i means:
//the jth element in y is aligned to the ith element in x if i>=0 
//the jth element in y is aligned to a gap in x if i==-1
bool get_initial5( double **r1, double **r2, double **xtm, double **ytm,
    bool **path, double **val,
    double **x, double **y, int xlen, int ylen, int *y2x,
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
                NWDP_TM(path, val, x, y, xlen, ylen,
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

void score_matrix_rmsd_sec( double **r1, double **r2, double **score,
    const char *secx, const char *secy, double **x, double **y,
    int xlen, int ylen, int *y2x, const double D0_MIN, double d0)
{
    double t[3], u[3][3];
    double rmsd, dij;
    double d01=d0+1.5;
    if(d01 < D0_MIN) d01=D0_MIN;
    double d02=d01*d01;

    double xx[3];
    int i, k=0;
    for(int j=0; j<ylen; j++)
    {
        i=y2x[j];
        if(i>=0)
        {
            r1[k][0]=x[i][0];  
            r1[k][1]=x[i][1]; 
            r1[k][2]=x[i][2];   
            
            r2[k][0]=y[j][0];  
            r2[k][1]=y[j][1]; 
            r2[k][2]=y[j][2];
            
            k++;
        }
    }
    Kabsch(r1, r2, k, 1, &rmsd, t, u);

    
    for(int ii=0; ii<xlen; ii++)
    {        
        transform(t, u, &x[ii][0], xx);
        for(int jj=0; jj<ylen; jj++)
        {
            dij=dist(xx, &y[jj][0]); 
            if (secx[ii]==secy[jj])
                score[ii+1][jj+1] = 1.0/(1+dij/d02) + 0.5;
            else
                score[ii+1][jj+1] = 1.0/(1+dij/d02);
        }
    }
}


//get initial alignment from secondary structure and previous alignments
//input: x, y, xlen, ylen
//output: y2x stores the best alignment: e.g., 
//y2x[j]=i means:
//the jth element in y is aligned to the ith element in x if i>=0 
//the jth element in y is aligned to a gap in x if i==-1
void get_initial_ssplus(double **r1, double **r2, double **score, bool **path,
    double **val, const char *secx, const char *secy, double **x, double **y,
    int xlen, int ylen, int *y2x0, int *y2x, const double D0_MIN, double d0)
{
    //create score matrix for DP
    score_matrix_rmsd_sec(r1, r2, score, secx, secy, x, y, xlen, ylen,
        y2x0, D0_MIN,d0);
    
    double gap_open=-1.0;
    NWDP_TM(score, path, val, xlen, ylen, gap_open, y2x);
}


void find_max_frag(double **x, int len, int *start_max,
    int *end_max, double dcu0, const bool fast_opt)
{
    int r_min, fra_min=4;           //minimum fragment for search
    if (fast_opt) fra_min=8;
    int start;
    int Lfr_max=0;

    r_min= (int) (len*1.0/3.0); //minimum fragment, in case too small protein
    if(r_min > fra_min) r_min=fra_min;
    
    int inc=0;
    double dcu0_cut=dcu0*dcu0;;
    double dcu_cut=dcu0_cut;

    while(Lfr_max < r_min)
    {        
        Lfr_max=0;            
        int j=1;    //number of residues at nf-fragment
        start=0;
        for(int i=1; i<len; i++)
        {
            if(dist(x[i-1], x[i]) < dcu_cut)
            {
                j++;

                if(i==(len-1))
                {
                    if(j > Lfr_max) 
                    {
                        Lfr_max=j;
                        *start_max=start;
                        *end_max=i;                        
                    }
                    j=1;
                }
            }
            else
            {
                if(j>Lfr_max) 
                {
                    Lfr_max=j;
                    *start_max=start;
                    *end_max=i-1;                                        
                }

                j=1;
                start=i;
            }
        }// for i;
        
        if(Lfr_max < r_min)
        {
            inc++;
            double dinc=pow(1.1, (double) inc) * dcu0;
            dcu_cut= dinc*dinc;
        }
    }//while <;    
}

//perform fragment gapless threading to find the best initial alignment
//input: x, y, xlen, ylen
//output: y2x0 stores the best alignment: e.g., 
//y2x0[j]=i means:
//the jth element in y is aligned to the ith element in x if i>=0 
//the jth element in y is aligned to a gap in x if i==-1
double get_initial_fgt(double **r1, double **r2, double **xtm, double **ytm,
    double **x, double **y, int xlen, int ylen, 
    int *y2x, double d0, double d0_search,
    double dcu0, const bool fast_opt, double t[3], double u[3][3])
{
    int fra_min=4;           //minimum fragment for search
    if (fast_opt) fra_min=8;
    int fra_min1=fra_min-1;  //cutoff for shift, save time

    int xstart=0, ystart=0, xend=0, yend=0;

    find_max_frag(x, xlen, &xstart, &xend, dcu0, fast_opt);
    find_max_frag(y, ylen, &ystart, &yend, dcu0, fast_opt);


    int Lx = xend-xstart+1;
    int Ly = yend-ystart+1;
    int *ifr, *y2x_;
    int L_fr=getmin(Lx, Ly);
    ifr= new int[L_fr];
    y2x_= new int[ylen+1];

    //select what piece will be used. The original implement may cause 
    //asymetry, but only when xlen==ylen and Lx==Ly
    //if L1=Lfr1 and L2=Lfr2 (normal proteins), it will be the same as initial1

    if(Lx<Ly || (Lx==Ly && xlen<ylen))
    {        
        for(int i=0; i<L_fr; i++) ifr[i]=xstart+i;
    }
    else if(Lx>Ly || (Lx==Ly && xlen>ylen))
    {        
        for(int i=0; i<L_fr; i++) ifr[i]=ystart+i;
    }
    else // solve asymetric for 1x5gA vs 2q7nA5
    {
        /* In this case, L0==xlen==ylen; L_fr==Lx==Ly */
        int L0=xlen;
        double tmscore, tmscore_max=-1;
        int i, j, k;
        int n1, n2;
        int min_len;
        int min_ali;

        /* part 1, normalized by xlen */
        for(i=0; i<L_fr; i++) ifr[i]=xstart+i;

        if(L_fr==L0)
        {
            n1= (int)(L0*0.1); //my index starts from 0
            n2= (int)(L0*0.89);
            j=0;
            for(i=n1; i<= n2; i++)
            {
                ifr[j]=ifr[i];
                j++;
            }
            L_fr=j;
        }

        int L1=L_fr;
        min_len=getmin(L1, ylen);    
        min_ali= (int) (min_len/2.5); //minimum size of considered fragment 
        if(min_ali<=fra_min1)  min_ali=fra_min1;    
        n1 = -ylen+min_ali; 
        n2 = L1-min_ali;

        for(k=n1; k<=n2; k+=(fast_opt)?3:1)
        {
            //get the map
            for(j=0; j<ylen; j++)
            {
                i=j+k;
                if(i>=0 && i<L1) y2x_[j]=ifr[i];
                else             y2x_[j]=-1;
            }

            //evaluate the map quickly in three iterations
            tmscore=get_score_fast(r1, r2, xtm, ytm, x, y, xlen, ylen, y2x_,
                d0, d0_search, t, u);

            if(tmscore>=tmscore_max)
            {
                tmscore_max=tmscore;
                for(j=0; j<ylen; j++) y2x[j]=y2x_[j];
            }
        }

        /* part 2, normalized by ylen */
        L_fr=Ly;
        for(i=0; i<L_fr; i++) ifr[i]=ystart+i;

        if (L_fr==L0)
        {
            n1= (int)(L0*0.1); //my index starts from 0
            n2= (int)(L0*0.89);

            j=0;
            for(i=n1; i<= n2; i++)
            {
                ifr[j]=ifr[i];
                j++;
            }
            L_fr=j;
        }

        int L2=L_fr;
        min_len=getmin(xlen, L2);    
        min_ali= (int) (min_len/2.5); //minimum size of considered fragment 
        if(min_ali<=fra_min1)  min_ali=fra_min1;    
        n1 = -L2+min_ali; 
        n2 = xlen-min_ali;

        for(k=n1; k<=n2; k++)
        {
            //get the map
            for(j=0; j<ylen; j++) y2x_[j]=-1;

            for(j=0; j<L2; j++)
            {
                i=j+k;
                if(i>=0 && i<xlen) y2x_[ifr[j]]=i;
            }
        
            //evaluate the map quickly in three iterations
            tmscore=get_score_fast(r1, r2, xtm, ytm,
                x, y, xlen, ylen, y2x_, d0,d0_search, t, u);
            if(tmscore>=tmscore_max)
            {
                tmscore_max=tmscore;
                for(j=0; j<ylen; j++) y2x[j]=y2x_[j];
            }
        }

        delete [] ifr;
        delete [] y2x_;
        return tmscore_max;
    }

    
    int L0=getmin(xlen, ylen); //non-redundant to get_initial1
    if(L_fr==L0)
    {
        int n1= (int)(L0*0.1); //my index starts from 0
        int n2= (int)(L0*0.89);

        int j=0;
        for(int i=n1; i<= n2; i++)
        {
            ifr[j]=ifr[i];
            j++;
        }
        L_fr=j;
    }


    //gapless threading for the extracted fragment
    double tmscore, tmscore_max=-1;

    if(Lx<Ly || (Lx==Ly && xlen<=ylen))
    {
        int L1=L_fr;
        int min_len=getmin(L1, ylen);    
        int min_ali= (int) (min_len/2.5);              //minimum size of considered fragment 
        if(min_ali<=fra_min1)  min_ali=fra_min1;    
        int n1, n2;
        n1 = -ylen+min_ali; 
        n2 = L1-min_ali;

        int i, j, k;
        for(k=n1; k<=n2; k+=(fast_opt)?3:1)
        {
            //get the map
            for(j=0; j<ylen; j++)
            {
                i=j+k;
                if(i>=0 && i<L1) y2x_[j]=ifr[i];
                else             y2x_[j]=-1;
            }

            //evaluate the map quickly in three iterations
            tmscore=get_score_fast(r1, r2, xtm, ytm, x, y, xlen, ylen, y2x_,
                d0, d0_search, t, u);

            if(tmscore>=tmscore_max)
            {
                tmscore_max=tmscore;
                for(j=0; j<ylen; j++) y2x[j]=y2x_[j];
            }
        }
    }
    else
    {
        int L2=L_fr;
        int min_len=getmin(xlen, L2);    
        int min_ali= (int) (min_len/2.5);              //minimum size of considered fragment 
        if(min_ali<=fra_min1)  min_ali=fra_min1;    
        int n1, n2;
        n1 = -L2+min_ali; 
        n2 = xlen-min_ali;

        int i, j, k;    

        for(k=n1; k<=n2; k++)
        {
            //get the map
            for(j=0; j<ylen; j++) y2x_[j]=-1;

            for(j=0; j<L2; j++)
            {
                i=j+k;
                if(i>=0 && i<xlen) y2x_[ifr[j]]=i;
            }
        
            //evaluate the map quickly in three iterations
            tmscore=get_score_fast(r1, r2, xtm, ytm,
                x, y, xlen, ylen, y2x_, d0,d0_search, t, u);
            if(tmscore>=tmscore_max)
            {
                tmscore_max=tmscore;
                for(j=0; j<ylen; j++) y2x[j]=y2x_[j];
            }
        }
    }    


    delete [] ifr;
    delete [] y2x_;
    return tmscore_max;
}

//heuristic run of dynamic programing iteratively to find the best alignment
//input: initial rotation matrix t, u
//       vectors x and y, d0
//output: best alignment that maximizes the TMscore, will be stored in invmap
double DP_iter(double **r1, double **r2, double **xtm, double **ytm,
    double **xt, bool **path, double **val, double **x, double **y,
    int xlen, int ylen, double t[3], double u[3][3], int invmap0[],
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
            NWDP_TM(path, val, x, y, xlen, ylen,
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


void output_pymol(const string xname, const string yname,
    const string fname_super, double t[3], double u[3][3], const int ter_opt, 
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
    int i;
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
}

void output_rasmol(const string xname, const string yname,
    const string fname_super, double t[3], double u[3][3], const int ter_opt,
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
    int i;
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
    string resi; // 4-character residue sequence number
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

/* extract rotation matrix based on TMscore8 */
void output_rotation_matrix(const char* fname_matrix,
    const double t[3], const double u[3][3])
{
    stringstream ss;
    ss << "------ The rotation matrix to rotate Structure_1 to Structure_2 ------\n";
    char dest[1000];
    sprintf(dest, "m %18s %14s %14s %14s\n", "t[m]", "u[m][0]", "u[m][1]", "u[m][2]");
    ss << string(dest);
    for (int k = 0; k < 3; k++)
    {
        sprintf(dest, "%d %18.10f %14.10f %14.10f %14.10f\n", k, t[k], u[k][0], u[k][1], u[k][2]);
        ss << string(dest);
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

//output the final results
void output_results(const string xname, const string yname,
    const string chainID1, const string chainID2,
    const int xlen, const int ylen, double t[3], double u[3][3],
    const double TM1, const double TM2,
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
        printf("\n(\":\" denotes residue pairs of d <%4.1f Angstrom, ", d0_out);
        printf("\".\" denotes other aligned residues)\n");
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

    if (strlen(fname_matrix)) output_rotation_matrix(fname_matrix, t, u);

    if (o_opt==1)
        output_pymol(xname, yname, fname_super, t, u, ter_opt,
            mm_opt, split_opt, mirror_opt, seqM, seqxA, seqyA,
            resi_vec1, resi_vec2, chainID1, chainID2);
    else if (o_opt==2)
        output_rasmol(xname, yname, fname_super, t, u, ter_opt,
            mm_opt, split_opt, mirror_opt, seqM, seqxA, seqyA,
            resi_vec1, resi_vec2, chainID1, chainID2,
            xlen, ylen, d0A, n_ali8, rmsd, TM1, Liden);
}

void output_mTMalign_results(const string xname, const string yname,
    const string chainID1, const string chainID2,
    const int xlen, const int ylen, double t[3], double u[3][3],
    const double TM1, const double TM2,
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
        printf("Average aligned length= %d, RMSD= %6.2f, Seq_ID=n_identical/n_aligned= %4.3f\n", n_ali8, rmsd, (n_ali8>0)?Liden/n_ali8:0);
        printf("Average TM-score= %6.5f (normalized by length of shorter structure: L=%d, d0=%.2f)\n", TM2, xlen, d0B);
        printf("Average TM-score= %6.5f (normalized by length of longer structure: L=%d, d0=%.2f)\n", TM1, ylen, d0A);

        if (a_opt==1)
            printf("Average TM-score= %6.5f (if normalized by average length of two structures: L=%.1f, d0=%.2f)\n", TM3, (xlen+ylen)*0.5, d0a);
        if (u_opt)
            printf("Average TM-score= %6.5f (normalized by average L=%.2f and d0=%.2f)\n", TM4, Lnorm_ass, d0u);
        if (d_opt)
            printf("Average TM-score= %6.5f (scaled by user-specified d0=%.2f, and L=%d)\n", TM5, d0_scale, ylen);
    
        //output alignment
        printf("In the following, seqID=n_identical/L.\n\n%s\n", seqM);
    }
    else if (outfmt_opt==1)
    {
        printf("%s\n", seqM);

        printf("# Lali=%d\tRMSD=%.2f\tseqID_ali=%.3f\n",
            n_ali8, rmsd, (n_ali8>0)?Liden/n_ali8:0);

        if (i_opt)
            printf("# User-specified initial alignment: TM=%.5lf\tLali=%4d\trmsd=%.3lf\n", TM_ali, L_ali, rmsd_ali);

        if(a_opt)
            printf("# TM-score=%.5f (normalized by average length of two structures: L=%.1f\td0=%.2f)\n", TM3, (xlen+ylen)*0.5, d0a);

        if(u_opt)
            printf("# TM-score=%.5f (normalized by average L=%.2f\td0=%.2f)\n", TM4, Lnorm_ass, d0u);

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

    if (strlen(fname_matrix)) output_rotation_matrix(fname_matrix, t, u);

    if (o_opt==1)
        output_pymol(xname, yname, fname_super, t, u, ter_opt,
            mm_opt, split_opt, mirror_opt, seqM, seqxA, seqyA,
            resi_vec1, resi_vec2, chainID1, chainID2);
    else if (o_opt==2)
        output_rasmol(xname, yname, fname_super, t, u, ter_opt,
            mm_opt, split_opt, mirror_opt, seqM, seqxA, seqyA,
            resi_vec1, resi_vec2, chainID1, chainID2,
            xlen, ylen, d0A, n_ali8, rmsd, TM1, Liden);
}

double standard_TMscore(double **r1, double **r2, double **xtm, double **ytm,
    double **xt, double **x, double **y, int xlen, int ylen, int invmap[],
    int& L_ali, double& RMSD, double D0_MIN, double Lnorm, double d0,
    double d0_search, double score_d8, double t[3], double u[3][3],
    const int mol_type)
{
    D0_MIN = 0.5;
    Lnorm = ylen;
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

    double tmscore;// collected alined residues from invmap
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
        else if (i != -1) PrintErrorAndQuit("Wrong map!\n");
    }
    L_ali = n_al;

    Kabsch(r1, r2, n_al, 0, &RMSD, t, u);
    RMSD = sqrt( RMSD/(1.0*n_al) );
    
    int temp_simplify_step = 1;
    int temp_score_sum_method = 0;
    d0_search = d0_input;
    double rms = 0.0;
    tmscore = TMscore8_search_standard(r1, r2, xtm, ytm, xt, n_al, t, u,
        temp_simplify_step, temp_score_sum_method, &rms, d0_input,
        score_d8, d0);
    tmscore = tmscore * n_al / (1.0*Lnorm);

    return tmscore;
}

/* copy the value of t and u into t0,u0 */
void copy_t_u(double t[3], double u[3][3], double t0[3], double u0[3][3])
{
    int i,j;
    for (i=0;i<3;i++)
    {
        t0[i]=t[i];
        for (j=0;j<3;j++) u0[i][j]=u[i][j];
    }
}

/* calculate approximate TM-score given rotation matrix */
double approx_TM(const int xlen, const int ylen, const int a_opt,
    double **xa, double **ya, double t[3], double u[3][3],
    const int invmap0[], const int mol_type)
{
    double Lnorm_0=ylen; // normalized by the second protein
    if (a_opt==-2 && xlen>ylen) Lnorm_0=xlen;      // longer
    else if (a_opt==-1 && xlen<ylen) Lnorm_0=xlen; // shorter
    else if (a_opt==1) Lnorm_0=(xlen+ylen)/2.;     // average
    
    double D0_MIN;
    double Lnorm;
    double d0;
    double d0_search;
    parameter_set4final(Lnorm_0, D0_MIN, Lnorm, d0, d0_search, mol_type);
    double TMtmp=0;
    double d;
    double xtmp[3]={0,0,0};

    for(int i=0,j=0; j<ylen; j++)
    {
        i=invmap0[j];
        if(i>=0)//aligned
        {
            transform(t, u, &xa[i][0], &xtmp[0]);
            d=sqrt(dist(&xtmp[0], &ya[j][0]));
            TMtmp+=1/(1+(d/d0)*(d/d0));
            //if (d <= score_d8) TMtmp+=1/(1+(d/d0)*(d/d0));
        }
    }
    TMtmp/=Lnorm_0;
    return TMtmp;
}

void clean_up_after_approx_TM(int *invmap0, int *invmap,
    double **score, bool **path, double **val, double **xtm, double **ytm,
    double **xt, double **r1, double **r2, const int xlen, const int minlen)
{
    delete [] invmap0;
    delete [] invmap;
    DeleteArray(&score, xlen+1);
    DeleteArray(&path, xlen+1);
    DeleteArray(&val, xlen+1);
    DeleteArray(&xtm, minlen);
    DeleteArray(&ytm, minlen);
    DeleteArray(&xt, xlen);
    DeleteArray(&r1, minlen);
    DeleteArray(&r2, minlen);
    return;
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
    if (Lnorm <= 40) ddcc=0.1;   //Lnorm was setted in parameter_set4search
    double local_d0_search = d0_search;

    //************************************************//
    //    get initial alignment from user's input:    //
    //    Stick to the initial alignment              //
    //************************************************//
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
    }

    /******************************************************/
    /*    get initial alignment with gapless threading    */
    /******************************************************/
    if (i_opt<=1)
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

    /* duplicate structure */
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

bool output_cp(const string&xname, const string&yname,
    const string &seqxA, const string &seqyA, const int outfmt_opt,
    int &left_num, int &right_num, int &left_aln_num, int &right_aln_num)
{
    int r;
    bool after_cp=false;
    for (r=0;r<seqxA.size();r++)
    {
        if (seqxA[r]=='*') after_cp=true;
        else 
        {
            if (after_cp)
            {
                right_aln_num++;
                right_num+=(seqxA[r]!='-');
            }
            else
            {
                left_aln_num++;
                left_num+=(seqxA[r]!='-');
            }
        }
    }
    if (after_cp==false)
    {
        if (outfmt_opt<=0) cout<<"No CP"<<endl;
        else if (outfmt_opt==1) cout<<"#No CP"<<endl;
        else if (outfmt_opt==2) cout<<"@"<<xname<<'\t'<<yname<<'\t'<<"No CP"<<endl;
    }
    else
    {
        if (outfmt_opt<=0) cout<<"CP point in structure_1 alignment: "<<left_aln_num<<'/'<<right_aln_num<<'\n'
            <<"CP point in structure_1: "<<left_num<<'/'<<right_num<<endl;
        else if (outfmt_opt==1) 
            cout<<"#CP_in_aln="<<left_aln_num<<'/'<<right_aln_num
               <<"\tCP_in_seq="<<left_num<<'/'<<right_num<<endl;
        else if (outfmt_opt==2) cout<<"@"<<xname<<'\t'<<yname<<'\t'<<left_aln_num
            <<'/'<<right_aln_num<<'\t'<<left_num<<'/'<<right_num<<endl;
    }
    return after_cp;
}
#endif
