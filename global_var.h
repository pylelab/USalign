const char *version="20160521";   //version 
 
 
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

int atomxlen, atomylen; // length of atoms
int *ia1, *ia2;         // ATOM indices, used for display
char  **aa1, **aa2;     // "N", or "CA", or "C" 
double **xyza1, **xyza2;// for input vectors xa[0...xlen-1][0..2], ya[0...ylen-1][0..2], just for display
char **ra1, **ra2;      // for the protein sequence 
int  *ir1, *ir2;        // residue numbers, used in fragment gapless threading 

char sequence[10][MAXLEN];// get value from alignment file
double TM_ali, rmsd_ali;  // TMscore and rmsd from standard_TMscore func, 
int L_ali;                // Aligned length from standard_TMscore func, 

char *ins1, *ins2, *ains1, *ains2;// flag characters for data read from PDB file, which begins with "ATOM", and locates at s(27), usually are spaces
int **nres1, **nres2;// number of atoms, nres1(i,j): the number of atoms for ith residue, j usually is 32 for a space


//argument variables
char out_reg[MAXLEN];
double Lnorm_ass, Lnorm_d0, d0_scale, d0A, d0B, d0u, d0a;
bool o_opt, a_opt, u_opt, d_opt, v_opt;
bool i_opt;// flags for -i, with user given initial alignment file
bool m_opt;// flags for -m, output rotation matrix
bool I_opt;// flags for -I, stick to user given initial alignment file

int  fast_level; // 0 - default, 1 - faster
bool f_opt; // flags for -f, fast but inaccurate alignment

double TM3, TM4, TM5;
