//global variables
double **xa, **ya;      //for input vectors xa[0...xlen-1][0..2] and
                        //ya[0...ylen-1][0..2], in general,
                        //ya is regarded as native structure 
                        //--> superpose xa onto ya
double **xtm, **ytm;    //for TMscore search engine
double **xt;            //for saving the superposed version of r_1 or xtm
double **r1, **r2;      //for Kabsch rotation
