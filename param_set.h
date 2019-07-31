/* These functions implement d0 normalization. The d0 for final TM-score
 * output is implemented by parameter_set4final. For both RNA alignment
 * and protein alignment, using d0 set by parameter_set4search yields
 * slightly better results during initial alignment-superposition iteration.
 */
#include <math.h>
#include "basic_fun.h"

void parameter_set4search(const int xlen, const int ylen,
    double &D0_MIN, double &Lnorm,
    double &score_d8, double &d0, double &d0_search, double &dcu0)
{
    //parameter initialization for searching: D0_MIN, Lnorm, d0, d0_search, score_d8
    D0_MIN=0.5; 
    dcu0=4.25;                       //update 3.85-->4.25
 
    Lnorm=getmin(xlen, ylen);        //normalize TMscore by this in searching
    if (Lnorm<=19)                    //update 15-->19
        d0=0.168;                   //update 0.5-->0.168
    else d0=(1.24*pow((Lnorm*1.0-15), 1.0/3)-1.8);
    D0_MIN=d0+0.8;              //this should be moved to above
    d0=D0_MIN;                  //update: best for search    

    d0_search=d0;
    if (d0_search>8)   d0_search=8;
    if (d0_search<4.5) d0_search=4.5;

    score_d8=1.5*pow(Lnorm*1.0, 0.3)+3.5; //remove pairs with dis>d8 during search & final
}

void parameter_set4final_C3prime(const double len, double &D0_MIN,
    double &Lnorm, double &d0, double &d0_search)
{
    D0_MIN=0.3; 
 
    Lnorm=len;            //normalize TMscore by this in searching
    if(Lnorm<=11) d0=0.3;
    else if(Lnorm>11&&Lnorm<=15) d0=0.4;
    else if(Lnorm>15&&Lnorm<=19) d0=0.5;
    else if(Lnorm>19&&Lnorm<=23) d0=0.6;
    else if(Lnorm>23&&Lnorm<30)  d0=0.7;
    else d0=(0.6*pow((Lnorm*1.0-0.5), 1.0/2)-2.5);

    d0_search=d0;
    if (d0_search>8)   d0_search=8;
    if (d0_search<4.5) d0_search=4.5;
}

void parameter_set4final(const double len, double &D0_MIN, double &Lnorm,
    double &d0, double &d0_search, const int mol_type)
{
    if (mol_type>0) // RNA
    {
        parameter_set4final_C3prime(len, D0_MIN, Lnorm,
            d0, d0_search);
        return;
    }
    D0_MIN=0.5; 
 
    Lnorm=len;            //normalize TMscore by this in searching
    if (Lnorm<=21) d0=0.5;          
    else d0=(1.24*pow((Lnorm*1.0-15), 1.0/3)-1.8);
    if (d0<D0_MIN) d0=D0_MIN;   
    d0_search=d0;
    if (d0_search>8)   d0_search=8;
    if (d0_search<4.5) d0_search=4.5;
}

void parameter_set4scale(const int len, const double d_s, double &Lnorm,
    double &d0, double &d0_search)
{
    d0=d_s;          
    Lnorm=len;            //normalize TMscore by this in searching
    d0_search=d0;
    if (d0_search>8)   d0_search=8;
    if (d0_search<4.5) d0_search=4.5;  
}
