//global variables
double score_d8,d0,d0_search,dcu0;//for TMscore search
double **score;                   //Input score table for dynamic programming
bool   **path;                    //for dynamic programming  
double **val;                     //for dynamic programming  
int    xlen, ylen, minlen;        //length of proteins
double **xa, **ya;      //for input vectors xa[0...xlen-1][0..2] and
                        //ya[0...ylen-1][0..2], in general,
                        //ya is regarded as native structure 
                        //--> superpose xa onto ya
int    *xresno, *yresno;//residue numbers, used in fragment gapless threading 
double **xtm, **ytm;    //for TMscore search engine
double **xt;            //for saving the superposed version of r_1 or xtm
char   *seqx, *seqy;    //for the protein sequence 
int    *secx, *secy;    //for the secondary structure 
double **r1, **r2;      //for Kabsch rotation 
double t[3], u[3][3];   //Kabsch translation vector and rotation matrix
