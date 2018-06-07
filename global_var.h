const char *TMalign_version="20180604";   //version 
 
 
//global variables
double D0_MIN;                    //for d0
double Lnorm;                     //normalization length
double score_d8,d0,d0_search,dcu0;//for TMscore search
double **score;            		  //Input score table for dynamic programming
bool   **path;                    //for dynamic programming  
double **val;                     //for dynamic programming  
int    xlen, ylen, minlen;        //length of proteins
int tempxlen, tempylen;
double **xa, **ya;      //for input vectors xa[0...xlen-1][0..2], ya[0...ylen-1][0..2]
                        //in general, ya is regarded as native structure --> superpose xa onto ya
int    *xresno, *yresno;//residue numbers, used in fragment gapless threading 
double **xtm, **ytm;    //for TMscore search engine
double **xt;            //for saving the superposed version of r_1 or xtm
char   *seqx, *seqy;    //for the protein sequence 
int    *secx, *secy;    //for the secondary structure 
double **r1, **r2;      //for Kabsch rotation 
double t[3], u[3][3];   //Kabsch translation vector and rotation matrix

char sequence[10][MAXLEN];// get value from alignment file
double TM_ali, rmsd_ali;  // TMscore and rmsd from standard_TMscore func, 
int L_ali;                // Aligned length from standard_TMscore func, 

//argument variables
char out_reg[MAXLEN];
double Lnorm_ass, Lnorm_d0, d0_scale, d0A, d0B, d0u, d0a;
bool o_opt, a_opt, u_opt, d_opt;
bool i_opt;// flags for -i, with user given initial alignment file
bool m_opt;// flags for -m, output rotation matrix
bool I_opt;// flags for -I, stick to user given initial alignment file
bool fast_opt; // flags for -fast, fast but inaccurate alignment

double TM3, TM4, TM5;
