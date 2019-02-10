/* header for Needleman-Wunsch global sequence alignment */
#ifndef NWalign_H
#define NWalign_H 1

#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <string.h>
#include <iomanip>
#include <stdlib.h>
#include <ctype.h>

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

const string aa_list="ARNDCQEGHILKMFPSTWYVBZX*";

const int gapopen_blastn=-5;
const int gapext_blastn=-2;

const int BLASTN[24][24]={
//a  c  g  t  u  *
{ 2,-3,-3,-3,-3,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},//a
{-3, 2,-3,-3,-3,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},//c
{-3,-3, 2,-3,-3,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},//g
{-3,-3,-3, 2, 2,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},//t
{-3,-3,-3, 2, 2,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},//u
{-3,-3,-3,-3,-3, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},//*
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
};

const string na_list="acgtu*";

/* convert amino acid to int */
inline int aa2int(char aa)
{
    for (int i=0;i<aa_list.length();i++) if (aa_list[i]==aa) return i;
    if (aa!=toupper(aa)) return aa2int(toupper(aa));
    return aa_list.length();
}

inline int na2int(char na)
{
    for (int i=0;i<na_list.length();i++) if (na_list[i]==na) return i;
    if (na!=tolower(na)) return na2int(tolower(na));
    return na_list.length();
}

void aa2int(const char *sequence, const int xlen, vector<int>&seqyint,
    const int mol_type)
{
    seqyint.assign(xlen,0);
    if (mol_type>0) // RNA
        for (int l=0;l<xlen;l++) seqyint[l]=na2int(sequence[l]);
    else // protein
        for (int l=0;l<xlen;l++) seqyint[l]=aa2int(sequence[l]);
    return;
}

/* initialize matrix in gotoh algorithm */
void init_gotoh_mat(vector<vector<int> >&JumpH, vector<vector<int> >&JumpV,
    vector<vector<int> >& P,vector<vector<int> >& S, 
    vector<vector<int> >& H, vector<vector<int> >& V,
    const int len1, const int len2, const int gapopen,const int gapext,
    const int glocal=0, const int alt_init=1)
{
    // fill first row/colum of JumpH,jumpV and path matrix P
    int i,j;
    for (i=0;i<len1+1;i++)
    {
        if (glocal<2) P[i][0]=4; // -
        JumpV[i][0]=i;
    }
    for (j=0;j<len2+1;j++)
    {
        if (glocal<1) P[0][j]=2; // |
        JumpH[0][j]=j;
    }
    if (glocal<2) for (i=1;i<len1+1;i++) S[i][0]=gapopen+gapext*(i-1);
    if (glocal<1) for (j=1;j<len2+1;j++) S[0][j]=gapopen+gapext*(j-1);
    if (alt_init==0)
    {
        for (i=1;i<len1+1;i++) H[i][0]=gapopen+gapext*(i-1);
        for (j=1;j<len2+1;j++) V[0][j]=gapopen+gapext*(j-1);
    }
    else
    {
        if (glocal<2) for (i=1;i<len1+1;i++) V[i][0]=gapopen+gapext*(i-1);
        if (glocal<1) for (j=1;j<len2+1;j++) H[0][j]=gapopen+gapext*(j-1);
        for (i=0;i<len1+1;i++) H[i][0]=-99999; // INT_MIN cause bug on ubuntu
        for (j=0;j<len2+1;j++) V[0][j]=-99999; // INT_MIN;
    }
}

/* locate the cell with highest alignment score. reset path after
 * the cell to zero */
void find_highest_align_score(
    const vector<vector<int> >& S, vector<vector<int> >& P,
    int &aln_score, const int len1,const int len2)
{
    // locate the cell with highest alignment score
    int max_aln_i=len1;
    int max_aln_j=len2;
    int i,j;
    for (i=0;i<len1+1;i++)
    {
        for (j=0;j<len2+1;j++)
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
    for (i=max_aln_i+1;i<len1+1;i++) for (j=0;j<len2+1;j++) P[i][j]=0;
    for (i=0;i<len1+1;i++) for (j=max_aln_j+1;j<len2+1;j++) P[i][j]=0;
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
int calculate_score_gotoh(
    const vector<int>& seqyint1, const vector<int>& seqyint2,
    vector<vector<int> >& JumpH, vector<vector<int> >& JumpV,
    vector<vector<int> >& P,const int ScoringMatrix[24][24],
    const int gapopen,const int gapext,const int glocal=0,
    const int alt_init=1)
{
    int len1=seqyint1.size();
    int len2=seqyint2.size();

    vector<int> temp_int(len2+1,0);
    vector<vector<int> > S(len1+1,temp_int);
    // penalty score for horizontal long gap
    vector<vector<int> > H(len1+1,temp_int);
    // penalty score for vertical long gap
    vector<vector<int> > V(len1+1,temp_int);
    
    // fill first row/colum of JumpH,jumpV and path matrix P
    int i,j;
    init_gotoh_mat(JumpH, JumpV, P, S, H, V, len1, len2,
        gapopen, gapext, glocal, alt_init);

    // fill S and P
    int diag_score,left_score,up_score;
    for (i=1;i<len1+1;i++)
    {
        for (j=1;j<len2+1;j++)
        {
            // penalty of consective deletion
            if (glocal<1 || i<len1 || glocal>=3)
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
            if (glocal<2 || j<len2 || glocal>=3)
            {
                V[i][j]=MAX(S[i-1][j]+gapopen,V[i-1][j]+gapext);
                JumpV[i][j]=(V[i][j]==V[i-1][j]+gapext)?(JumpV[i-1][j]+1):1;
            }
            else
            {
                V[i][j]=MAX(S[i-1][j],V[i-1][j]);
                JumpV[i][j]=(V[i][j]==V[i-1][j])?(JumpV[i-1][j]+1):1;
            }

            diag_score=S[i-1][j-1]; // match-mismatch '\'
            if (seqyint1[i-1]<24 && seqyint2[j-1]<24)
                diag_score+=ScoringMatrix[seqyint1[i-1]][seqyint2[j-1]];
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
    int aln_score=S[len1][len2];

    // re-fill first row/column of path matrix P for back-tracing
    for (i=1;i<len1+1;i++) if (glocal<3 || P[i][0]>0) P[i][0]=2; // |
    for (j=1;j<len2+1;j++) if (glocal<3 || P[0][j]>0) P[0][j]=4; // -

    // calculate alignment score and alignment path for swalign
    if (glocal>=3)
        find_highest_align_score(S,P,aln_score,len1,len2);

    // release memory
    S.clear();
    H.clear();
    V.clear();
    return aln_score; // final alignment score
}

/* trace back dynamic programming path to diciper pairwise alignment */
void trace_back_gotoh(string seqx, string seqy,
    const vector<vector<int> >& JumpH, const vector<vector<int> >& JumpV,
    const vector<vector<int> >& P, string& seqxA, string& seqyA)
{
    int len1=seqx.length();
    int len2=seqy.length();
    
    int i=len1;
    int j=len2;
    int gaplen,p;

    while(i+j)
    {
        gaplen=0;
        if (P[i][j]>=4)
        {
            gaplen=JumpH[i][j];
            for (p=0;p<gaplen;p++) seqxA='-'+seqxA;
            seqyA=seqy.substr(seqy.length()-gaplen,gaplen)+seqyA;
            seqy=seqy.substr(0,seqy.length()-gaplen);
            j-=gaplen;
        }
        else if (P[i][j] % 4 >= 2)
        {
            gaplen=JumpV[i][j];
            seqxA=seqx.substr(seqx.length()-gaplen,gaplen)+seqxA;
            for (p=0;p<gaplen;p++) seqyA='-'+seqyA;
            seqx=seqx.substr(0,seqx.length()-gaplen);
            i-=gaplen;
        }
        else
        {
            if (i==0 && j!=0) // only in glocal alignment
            {
                seqyA=seqy+seqyA;
                for (p=0;p<seqy.length();p++) seqxA='-'+seqxA;
                break;
            }
            if (i!=0 && j==0) // only in glocal alignment
            {
                seqxA=seqx+seqxA;
                for (p=0;p<seqx.length();p++) seqyA='-'+seqyA;
                break;
            }
            seqxA=seqx[seqx.length()-1]+seqxA;
            seqyA=seqy[seqy.length()-1]+seqyA;
            seqx=seqx.substr(0,seqx.length()-1);
            seqy=seqy.substr(0,seqy.length()-1);
            i--;
            j--;
        }
    }   
}


/* trace back Smith-Waterman dynamic programming path to diciper 
 * pairwise local alignment */
void trace_back_sw(string seqx, string seqy,
    const vector<vector<int> >& JumpH, const vector<vector<int> >& JumpV,
    const vector<vector<int> >& P, string& seqxA, string& seqyA)
{
    int len1=seqx.length();
    int len2=seqy.length();
    
    int i=len1;
    int j=len2;
    // find the first non-zero cell in P
    bool found_start_cell=false;
    for (i=len1;i>=0;i--)
    {
        for (j=len2;j>=0;j--)
        {
            if (P[i][j]!=0)
            {
                found_start_cell=true;
                break;
            }
        }
        if (found_start_cell) break;
    }
    if (i<0||j<0) return;
    seqx=seqx.substr(0,i);
    seqy=seqy.substr(0,j);

    int gaplen,p;
    while(P[i][j]!=0)
    {
        gaplen=0;
        if (P[i][j]>=4)
        {
            gaplen=JumpH[i][j];
            for (p=0;p<gaplen;p++) seqxA='-'+seqxA;
            seqyA=seqy.substr(seqy.length()-gaplen,gaplen)+seqyA;
            seqy=seqy.substr(0,seqy.length()-gaplen);
            j-=gaplen;
        }
        else if (P[i][j] % 4 >= 2)
        {
            gaplen=JumpV[i][j];
            seqxA=seqx.substr(seqx.length()-gaplen,gaplen)+seqxA;
            for (p=0;p<gaplen;p++) seqyA='-'+seqyA;
            seqx=seqx.substr(0,seqx.length()-gaplen);
            i-=gaplen;
        }
        else
        {
            seqxA=seqx[seqx.length()-1]+seqxA;
            seqyA=seqy[seqy.length()-1]+seqyA;
            seqx=seqx.substr(0,seqx.length()-1);
            seqy=seqy.substr(0,seqy.length()-1);
            i--;
            j--;
        }
    }
}

/* entry function for NWalign */
int NWalign(const char *seqx, const char *seqy, 
    const vector<int>& seqyint1, const vector<int>& seqyint2, // aa2int
    string & seqxA,string & seqyA, const int mol_type, const int glocal=0)
{
    int len1=seqyint1.size();
    int len2=seqyint2.size();
    vector<int> temp_int(len2+1,0);
    vector<vector<int> > JumpH(len1+1,temp_int);
    vector<vector<int> > JumpV(len1+1,temp_int);
    vector<vector<int> > P(len1+1,temp_int);
    
    int aln_score;
    if (mol_type>0) aln_score=calculate_score_gotoh(seqyint1,seqyint2,
            JumpH,JumpV,P, BLASTN,gapopen_blastn,gapext_blastn,glocal);
    else            aln_score=calculate_score_gotoh(seqyint1,seqyint2,
            JumpH,JumpV,P, BLOSUM62,gapopen_blosum62,gapext_blosum62,glocal);

    if (glocal<3) trace_back_gotoh(seqx,seqy,JumpH,JumpV,P,seqxA,seqyA);
    else trace_back_sw(seqx,seqy,JumpH,JumpV,P,seqxA,seqyA);

    JumpH.clear();
    JumpV.clear();
    P.clear();
    return aln_score; // aligment score
}

double get_seqID(const string& seqxA, const string& seqyA,
    string &seqM,double &Liden,int &L_ali)
{
    Liden=0;
    L_ali=0;
    for (int i=0;i<seqxA.length();i++)
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
    const int xlen, const int ylen,
    const char *seqM, const char *seqxA, const char *seqyA,
    const double Liden, const int L_ali, const int outfmt_opt)
{
    if (outfmt_opt<=0)
    {
        printf("\nName of Chain_1: %s%s\n", xname.c_str(), chainID1);
        printf("Name of Chain_2: %s%s\n", yname.c_str(), chainID2);
        printf("Length of Chain_1: %d residues\n", xlen);
        printf("Length of Chain_2: %d residues\n\n", ylen);

        printf("Aligned length= %d, Seq_ID=n_identical/n_aligned= %4.3f\n", L_ali, Liden/L_ali);
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
