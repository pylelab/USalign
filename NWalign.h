/* header for Needleman-Wunsch global sequence alignment */
#ifndef NWalign_H
#define NWalign_H 1

#include "basic_fun.h"

#define MAX(A,B) ((A)>(B)?(A):(B))

using namespace std;

const int gapopen_blosum62=-11;
const int gapext_blosum62=-1;

const int BLOSUM62[24][24]={
//A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
{ 4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0,-2,-1, 0,-4},//A
{-1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3,-1, 0,-1,-4},//R
{-2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3, 3, 0,-1,-4},//N
{-2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3, 4, 1,-1,-4},//D
{ 0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,-3,-3,-2,-4},//C
{-1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2, 0, 3,-1,-4},//Q
{-1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2, 1, 4,-1,-4},//E
{ 0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3,-1,-2,-1,-4},//G
{-2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3, 0, 0,-1,-4},//H
{-1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3,-3,-3,-1,-4},//I
{-1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1,-4,-3,-1,-4},//L
{-1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2, 0, 1,-1,-4},//K
{-1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1,-3,-1,-1,-4},//M
{-2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1,-3,-3,-1,-4},//F
{-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2,-2,-1,-2,-4},//P
{ 1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2, 0, 0, 0,-4},//S
{ 0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0,-1,-1, 0,-4},//T
{-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3,-4,-3,-2,-4},//W
{-2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1,-3,-2,-1,-4},//Y
{ 0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4,-3,-2,-1,-4},//V
{-2,-1, 3, 4,-3, 0, 1,-1, 0,-3,-4, 0,-3,-3,-2, 0,-1,-4,-3,-3, 4, 1,-1,-4},//B
{-1, 0, 0, 1,-3, 3, 4,-2, 0,-3,-3, 1,-1,-3,-1, 0,-1,-3,-2,-2, 1, 4,-1,-4},//Z
{ 0,-1,-1,-1,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-2, 0, 0,-2,-1,-1,-1,-1,-1,-4},//X
{-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4, 1},//*
};

const int gapopen_blastn=-5;
const int gapext_blastn=-2;

const int BLASTN[5][5]={
//a  c  g  t  u
{ 2,-3,-3,-3,-3},//a
{-3, 2,-3,-3,-3},//c
{-3,-3, 2,-3,-3},//g
{-3,-3,-3, 2, 2},//t
{-3,-3,-3, 2, 2},//u
};

const int aa_ascii2int[128]={
    23, /*0   '\0'*/ 23, /*1   SOH */ 23, /*2   STX */ 23, /*3   ETX */
    23, /*4   EOT */ 23, /*5   ENQ */ 23, /*6   ACK */ 23, /*7   '\a'*/
    23, /*8   '\b'*/ 23, /*9   '\t'*/ 23, /*10  '\n'*/ 23, /*11  '\v'*/
    23, /*12  '\f'*/ 23, /*13  '\r'*/ 23, /*14  SO  */ 23, /*15  SI  */
    23, /*16  DLE */ 23, /*17  DC1 */ 23, /*18  DC2 */ 23, /*19  DC3 */
    23, /*20  DC4 */ 23, /*21  NAK */ 23, /*22  SYN */ 23, /*23  ETB */
    23, /*24  CAN */ 23, /*25  EM  */ 23, /*26  SUB */ 23, /*27  ESC */
    23, /*28  FS  */ 23, /*29  GS  */ 23, /*30  RS  */ 23, /*31  US  */
    23, /*32  ' ' */ 23, /*33  !   */ 23, /*34  "   */ 23, /*35  #   */
    23, /*36  $   */ 23, /*37  %   */ 23, /*38  &   */ 23, /*39  '   */
    23, /*40  (   */ 23, /*41  )   */ 23, /*42  *   */ 23, /*43  +   */
    23, /*44  ,   */ 23, /*45  -   */ 23, /*46  .   */ 23, /*47  /   */
    23, /*48  0   */ 23, /*49  1   */ 23, /*50  2   */ 23, /*51  3   */
    23, /*52  4   */ 23, /*53  5   */ 23, /*54  6   */ 23, /*55  7   */
    23, /*56  8   */ 23, /*57  9   */ 23, /*58  :   */ 23, /*59  ;   */
    23, /*60  <   */ 23, /*61  =   */ 23, /*62  >   */ 23, /*63  ?   */
    23, /*64  @   */  0, /*65  A   */ 20, /*66  B   */  4, /*67  C   */
     3, /*68  D   */  6, /*69  E   */ 13, /*70  F   */  7, /*71  G   */
     8, /*72  H   */  9, /*73  I   */ 23, /*74  J   */ 11, /*75  K   */
    10, /*76  L   */ 12, /*77  M   */  2, /*78  N   */ 23, /*79  O   */
    14, /*80  P   */  5, /*81  Q   */  1, /*82  R   */ 15, /*83  S   */
    16, /*84  T   */ 23, /*85  U   */ 19, /*86  V   */ 17, /*87  W   */
    22, /*88  X   */ 18, /*89  Y   */ 21, /*90  Z   */ 23, /*91  [   */
    23, /*92  \   */ 23, /*93  ]   */ 23, /*94  ^   */ 23, /*95  _   */
    23, /*96  `   */  0, /*97  a   */ 23, /*98  b   */  1, /*99  c   */
    23, /*100 d   */ 23, /*101 e   */ 23, /*102 f   */  2, /*103 g   */
    23, /*104 h   */ 23, /*105 i   */ 23, /*106 j   */ 23, /*107 k   */
    23, /*108 l   */ 23, /*109 m   */ 23, /*110 n   */ 23, /*111 o   */
    23, /*112 p   */ 23, /*113 q   */ 23, /*114 r   */ 23, /*115 s   */
     3, /*116 t   */  4, /*117 u   */ 23, /*118 v   */ 23, /*119 w   */
    23, /*120 x   */ 23, /*121 y   */ 23, /*122 z   */ 23, /*123 {   */
    23, /*124 |   */ 23, /*125 }   */ 23, /*126 ~   */ 23, /*127 DEL */
};

/* initialize matrix in gotoh algorithm */
void init_gotoh_mat(int **S, int **JumpH, int **JumpV, int **P,
    int **H, int **V, const int xlen, const int ylen, const int gapopen,
    const int gapext, const int glocal=0, const int alt_init=1)
{
    // fill first row/colum of JumpH,jumpV and path matrix P
    int i,j;
    for (i=0;i<xlen+1;i++)
        for (j=0;j<ylen+1;j++)
            H[i][j]=V[i][j]=P[i][j]=JumpH[i][j]=JumpV[i][j]=0;
    for (i=0;i<xlen+1;i++)
    {
        if (glocal<2) P[i][0]=4; // -
        JumpV[i][0]=i;
    }
    for (j=0;j<ylen+1;j++)
    {
        if (glocal<1) P[0][j]=2; // |
        JumpH[0][j]=j;
    }
    if (glocal<2) for (i=1;i<xlen+1;i++) S[i][0]=gapopen+gapext*(i-1);
    if (glocal<1) for (j=1;j<ylen+1;j++) S[0][j]=gapopen+gapext*(j-1);
    if (alt_init==0)
    {
        for (i=1;i<xlen+1;i++) H[i][0]=gapopen+gapext*(i-1);
        for (j=1;j<ylen+1;j++) V[0][j]=gapopen+gapext*(j-1);
    }
    else
    {
        if (glocal<2) for (i=1;i<xlen+1;i++) V[i][0]=gapopen+gapext*(i-1);
        if (glocal<1) for (j=1;j<ylen+1;j++) H[0][j]=gapopen+gapext*(j-1);
        for (i=0;i<xlen+1;i++) H[i][0]=-99999; // INT_MIN cause bug on ubuntu
        for (j=0;j<ylen+1;j++) V[0][j]=-99999; // INT_MIN;
    }
}

/* locate the cell with highest alignment score. reset path after
 * the cell to zero */
void find_highest_align_score( int **S, int **P,
    int &aln_score, const int xlen,const int ylen)
{
    // locate the cell with highest alignment score
    int max_aln_i=xlen;
    int max_aln_j=ylen;
    int i,j;
    for (i=0;i<xlen+1;i++)
    {
        for (j=0;j<ylen+1;j++)
        {
            if (S[i][j]>=aln_score)
            {
                max_aln_i=i;
                max_aln_j=j;
                aln_score=S[i][j];
            }
        }
    }

    // reset all path after [max_aln_i][max_aln_j]
    for (i=max_aln_i+1;i<xlen+1;i++) for (j=0;j<ylen+1;j++) P[i][j]=0;
    for (i=0;i<xlen+1;i++) for (j=max_aln_j+1;j<ylen+1;j++) P[i][j]=0;
}

/* calculate dynamic programming matrix using gotoh algorithm
 * S     - cumulative scorefor each cell
 * P     - string representation for path
 *         0 :   uninitialized, for gaps at N- & C- termini when glocal>0
 *         1 : \ match-mismatch
 *         2 : | vertical gap (insertion)
 *         4 : - horizontal gap (deletion)
 * JumpH - horizontal long gap number.
 * JumpV - vertical long gap number.
 * all matrices are in the size of [len(seqx)+1]*[len(seqy)+1]
 *
 * glocal - global or local alignment
 *         0 : global alignment (Needleman-Wunsch dynamic programming)
 *         1 : glocal-query alignment
 *         2 : glocal-both alignment
 *         3 : local alignment (Smith-Waterman dynamic programming)
 *
 * alt_init - whether to adopt alternative matrix initialization
 *         1 : use wei zheng's matrix initialization
 *         0 : use yang zhang's matrix initialization, does NOT work
 *             for glocal alignment
 */
int calculate_score_gotoh(const int xlen,const int ylen, int **S,
    int** JumpH, int** JumpV, int **P, const int gapopen,const int gapext,
    const int glocal=0, const int alt_init=1)
{
    int **H;
    int **V;
    NewArray(&H,xlen+1,ylen+1); // penalty score for horizontal long gap
    NewArray(&V,xlen+1,ylen+1); // penalty score for vertical long gap
    
    // fill first row/colum of JumpH,jumpV and path matrix P
    int i,j;
    init_gotoh_mat(S, JumpH, JumpV, P, H, V, xlen, ylen,
        gapopen, gapext, glocal, alt_init);

    // fill S and P
    int diag_score,left_score,up_score;
    for (i=1;i<xlen+1;i++)
    {
        for (j=1;j<ylen+1;j++)
        {
            // penalty of consective deletion
            if (glocal<1 || i<xlen || glocal>=3)
            {
                H[i][j]=MAX(S[i][j-1]+gapopen,H[i][j-1]+gapext);
                JumpH[i][j]=(H[i][j]==H[i][j-1]+gapext)?(JumpH[i][j-1]+1):1;
            }
            else
            {
                H[i][j]=MAX(S[i][j-1],H[i][j-1]);
                JumpH[i][j]=(H[i][j]==H[i][j-1])?(JumpH[i][j-1]+1):1;
            }
            // penalty of consective insertion
            if (glocal<2 || j<ylen || glocal>=3)
            {
                V[i][j]=MAX(S[i-1][j]+gapopen,V[i-1][j]+gapext);
                JumpV[i][j]=(V[i][j]==V[i-1][j]+gapext)?(JumpV[i-1][j]+1):1;
            }
            else
            {
                V[i][j]=MAX(S[i-1][j],V[i-1][j]);
                JumpV[i][j]=(V[i][j]==V[i-1][j])?(JumpV[i-1][j]+1):1;
            }

            diag_score=S[i-1][j-1]+S[i][j]; // match-mismatch '\'
            left_score=H[i][j];     // deletion       '-'
            up_score  =V[i][j];     // insertion      '|'

            if (diag_score>=left_score && diag_score>=up_score)
            {
                S[i][j]=diag_score;
                P[i][j]+=1;
            }
            if (up_score>=diag_score && up_score>=left_score)
            {
                S[i][j]=up_score;
                P[i][j]+=2;
            }
            if (left_score>=diag_score && left_score>=up_score)
            {
                S[i][j]=left_score;
                P[i][j]+=4;
            }
            if (glocal>=3 && S[i][j]<0)
            {
                S[i][j]=0;
                P[i][j]=0;
                H[i][j]=0;
                V[i][j]=0;
                JumpH[i][j]=0;
                JumpV[i][j]=0;
            }
        }
    }
    int aln_score=S[xlen][ylen];

    // re-fill first row/column of path matrix P for back-tracing
    for (i=1;i<xlen+1;i++) if (glocal<3 || P[i][0]>0) P[i][0]=2; // |
    for (j=1;j<ylen+1;j++) if (glocal<3 || P[0][j]>0) P[0][j]=4; // -

    // calculate alignment score and alignment path for swalign
    if (glocal>=3)
        find_highest_align_score(S,P,aln_score,xlen,ylen);

    // release memory
    DeleteArray(&H,xlen+1);
    DeleteArray(&V,xlen+1);
    return aln_score; // final alignment score
}

/* trace back dynamic programming path to diciper pairwise alignment */
void trace_back_gotoh(const char *seqx, const char *seqy,
    int ** JumpH, int ** JumpV, int ** P,
    string& seqxA, string& seqyA, const int xlen, const int ylen)
{
    int i=xlen;
    int j=ylen;
    int gaplen,p;
    char *buf=new char [MAX(xlen,ylen)+1];

    while(i+j)
    {
        gaplen=0;
        if (P[i][j]>=4)
        {
            gaplen=JumpH[i][j];
            j-=gaplen;
            strncpy(buf,seqy+j,gaplen);
            buf[gaplen]=0;
            seqyA=buf+seqyA;

            for (p=0;p<gaplen;p++) buf[p]='-';
            seqxA=buf+seqxA;
        }
        else if (P[i][j] % 4 >= 2)
        {
            gaplen=JumpV[i][j];
            i-=gaplen;
            strncpy(buf,seqx+i,gaplen);
            buf[gaplen]=0;
            seqxA=buf+seqxA;

            for (p=0;p<gaplen;p++) buf[p]='-';
            seqyA=buf+seqyA;
        }
        else
        {
            if (i==0 && j!=0) // only in glocal alignment
            {
                strncpy(buf,seqy,j);
                buf[j]=0;
                seqyA=buf+seqyA;
                for (p=0;p<j;p++) buf[p]='-';
                seqxA=buf+seqxA;
                break;
            }
            if (i!=0 && j==0) // only in glocal alignment
            {
                strncpy(buf,seqx,i);
                buf[i]=0;
                seqxA=buf+seqxA;
                for (p=0;p<i;p++) buf[p]='-';
                seqyA=buf+seqyA;
                break;
            }
            i--;
            j--;
            seqxA=seqx[i]+seqxA;
            seqyA=seqy[j]+seqyA;
        }
    }   
    delete [] buf;
}


/* trace back Smith-Waterman dynamic programming path to diciper 
 * pairwise local alignment */
void trace_back_sw(const char *seqx, const char *seqy,
    int **JumpH, int **JumpV, int **P,
    string& seqxA, string& seqyA, const int xlen, const int ylen)
{
    int i=xlen;
    int j=ylen;
    int gaplen,p;
    char *buf=new char [xlen+ylen+1];

    // find the first non-zero cell in P
    bool found_start_cell=false;
    for (i=xlen;i>=0;i--)
    {
        for (j=ylen;j>=0;j--)
        {
            if (P[i][j]!=0)
            {
                found_start_cell=true;
                break;
            }
        }
        if (found_start_cell) break;
    }

    /* copy C terminal sequence */
    for (p=0;p<ylen-j;p++) buf[p]='-';
    buf[ylen-j]=0;
    seqxA=buf;
    strncpy(buf,seqx+i,xlen-i);
    buf[xlen-i]=0;
    seqxA+=buf;

    strncpy(buf,seqy+j,ylen-j);
    buf[ylen-j]=0;
    seqyA+=buf;
    for (p=0;p<xlen-i;p++) buf[p]='-';
    buf[xlen-i]=0;
    seqyA+=buf;

    if (i<0||j<0)
    {
        delete [] buf;
        return;
    }

    /* traceback aligned sequences */
    while(P[i][j]!=0)
    {
        gaplen=0;
        if (P[i][j]>=4)
        {
            gaplen=JumpH[i][j];
            j-=gaplen;
            strncpy(buf,seqy+j,gaplen);
            buf[gaplen]=0;
            seqyA=buf+seqyA;

            for (p=0;p<gaplen;p++) buf[p]='-';
            seqxA=buf+seqxA;
        }
        else if (P[i][j] % 4 >= 2)
        {
            gaplen=JumpV[i][j];
            i-=gaplen;
            strncpy(buf,seqx+i,gaplen);
            buf[gaplen]=0;
            seqxA=buf+seqxA;

            for (p=0;p<gaplen;p++) buf[p]='-';
            seqyA=buf+seqyA;
        }
        else
        {
            i--;
            j--;
            seqxA=seqx[i]+seqxA;
            seqyA=seqy[j]+seqyA;
        }
    }
    /* copy N terminal sequence */
    for (p=0;p<j;p++) buf[p]='-';
    strncpy(buf+j,seqx,i);
    buf[i+j]=0;
    seqxA=buf+seqxA;

    strncpy(buf,seqy,j);
    for (p=j;p<j+i;p++) buf[p]='-';
    buf[i+j]=0;
    seqyA=buf+seqyA;
    delete [] buf;
}

/* entry function for NWalign */
int NWalign(const char *seqx, const char *seqy, const int xlen, const int ylen,
    string & seqxA,string & seqyA, const int mol_type, const int glocal=0)
{
    int **JumpH;
    int **JumpV;
    int **P;
    int **S;
    NewArray(&JumpH,xlen+1,ylen+1);
    NewArray(&JumpV,xlen+1,ylen+1);
    NewArray(&P,xlen+1,ylen+1);
    NewArray(&S,xlen+1,ylen+1);
    
    int aln_score;
    int gapopen=gapopen_blosum62;
    int gapext =gapext_blosum62;
    int i,j;
    if (mol_type>0) // RNA or DNA
    {
        gapopen=gapopen_blastn;
        gapext =gapext_blastn;
    }
    int aa1=23;
    int aa2=23;
    for (i=0;i<xlen+1;i++)
    {
        if (i) aa1=aa_ascii2int[seqx[i-1]];
        for (j=0;j<ylen+1;j++)
        {
            S[i][j]=0;
            if (i*j==0) continue;
            aa2=aa_ascii2int[seqy[j-1]];
            if (mol_type>0 && aa1<5 && aa2<5) S[i][j]=BLASTN[aa1][aa2];
            else if (mol_type<=0)           S[i][j]=BLOSUM62[aa1][aa2];
        }
    }

    aln_score=calculate_score_gotoh(xlen, ylen, S, JumpH, JumpV, P,
        gapopen, gapext, glocal);

    if (glocal<3) trace_back_gotoh(seqx,seqy,JumpH,JumpV,P,seqxA,seqyA,xlen,ylen);
    else trace_back_sw(seqx,seqy,JumpH,JumpV,P,seqxA,seqyA,xlen,ylen);
    
    DeleteArray(&JumpH, xlen+1);
    DeleteArray(&JumpV, xlen+1);
    DeleteArray(&P, xlen+1);
    DeleteArray(&S, xlen+1);
    return aln_score; // aligment score
}

double get_seqID(const string& seqxA, const string& seqyA,
    string &seqM,double &Liden,int &L_ali)
{
    Liden=0;
    L_ali=0;
    for (int i=0;i<seqxA.size();i++)
    {
        if (seqxA[i]==seqyA[i] && seqxA[i]!='-')
        {
            Liden++;
            seqM+=':';
        }
        else seqM+=' ';
        L_ali+=(seqxA[i]!='-' && seqyA[i]!='-');
    }
    return 1.*Liden/L_ali;
}


void output_NWalign_results(
    const string xname, const string yname,
    const char *chainID1, const char *chainID2,
    const int xlen, const int ylen, const char *seqM, 
    const char *seqxA, const char *seqyA, const double Liden,
    const int L_ali, const int aln_score, const int outfmt_opt)
{
    if (outfmt_opt<=0)
    {
        printf("\nName of Chain_1: %s%s\n", xname.c_str(), chainID1);
        printf("Name of Chain_2: %s%s\n", yname.c_str(), chainID2);
        printf("Length of Chain_1: %d residues\n", xlen);
        printf("Length of Chain_2: %d residues\n\n", ylen);

        printf("Aligned length= %d, Alignment score= %d, Seq_ID=n_identical/n_aligned= %4.3f\n",
            L_ali, aln_score, Liden/L_ali);
        printf("Seq_ID= %6.5f (if normalized by length of Chain_1\n", Liden/xlen);
        printf("Seq_ID= %6.5f (if normalized by length of Chain_2\n", Liden/ylen);
        printf("(You should use Seq_ID normalized by length of the reference structure)\n");
    
        //output alignment
        printf("\n(\":\" denotes pairs with identical residue type)\n");
        printf("%s\n", seqxA);
        printf("%s\n", seqM);
        printf("%s\n", seqyA);
    }
    else if (outfmt_opt==1)
    {
        printf(">%s%s\tL=%d\tseqID=%.3f\n",
            xname.c_str(), chainID1, xlen, Liden/xlen);
        printf("%s\n", seqxA);
        printf(">%s%s\tL=%d\tseqID=%.3f\n",
            yname.c_str(), chainID2, ylen, Liden/ylen);
        printf("%s\n", seqyA);
        printf("# Lali=%d\tseqID_ali=%.3f\n", L_ali, Liden/L_ali);
        printf("$$$$\n");
    }
    else if (outfmt_opt==2)
    {
        printf("%s%s\t%s%s\t%4.3f\t%4.3f\t%4.3f\t%d\t%d\t%d",
            xname.c_str(), chainID1, yname.c_str(), chainID2, 
            Liden/xlen, Liden/ylen, Liden/L_ali,
            xlen, ylen, L_ali);
    }
    cout << endl;
}

#endif
