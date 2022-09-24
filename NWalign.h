/* header for Needleman-Wunsch global sequence alignment */
#ifndef NWalign_H
#define NWalign_H 1

#include "basic_fun.h"
#include "BLOSUM.h"

#define MAX(A,B) ((A)>(B)?(A):(B))

using namespace std;

const int gapopen_blosum62=-11;
const int gapext_blosum62=-1;

const int gapopen_blastn=-15; //-5;
const int gapext_blastn =-4;  //-2;

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
    int ** JumpH, int ** JumpV, int ** P, string& seqxA, string& seqyA,
    const int xlen, const int ylen, int *invmap, const int invmap_only=1)
{
    int i,j;
    int gaplen,p;
    char *buf=NULL;

    if (invmap_only) for (j = 0; j < ylen; j++) invmap[j] = -1;
    if (invmap_only!=1) buf=new char [MAX(xlen,ylen)+1];

    i=xlen;
    j=ylen;
    while(i+j)
    {
        gaplen=0;
        if (P[i][j]>=4)
        {
            gaplen=JumpH[i][j];
            j-=gaplen;
            if (invmap_only==1) continue;
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
            if (invmap_only==1) continue;
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
            if (invmap_only) invmap[j]=i;
            if (invmap_only!=1)
            {
                seqxA=seqx[i]+seqxA;
                seqyA=seqy[j]+seqyA;
            }
        }
    }
    delete [] buf;
}


/* trace back Smith-Waterman dynamic programming path to diciper 
 * pairwise local alignment */
void trace_back_sw(const char *seqx, const char *seqy,
    int **JumpH, int **JumpV, int **P, string& seqxA, string& seqyA,
    const int xlen, const int ylen, int *invmap, const int invmap_only=1)
{
    int i;
    int j;
    int gaplen,p;
    bool found_start_cell=false; // find the first non-zero cell in P
    char *buf=NULL;

    if (invmap_only) for (j = 0; j < ylen; j++) invmap[j] = -1;
    if (invmap_only!=1) buf=new char [MAX(xlen,ylen)+1];

    i=xlen;
    j=ylen;
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
    if (invmap_only!=1)
    {
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
    }

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
            if (invmap_only==1) continue;
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
            if (invmap_only==1) continue;
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
            if (invmap_only) invmap[j]=i;
            if (invmap_only!=1)
            {
                seqxA=seqx[i]+seqxA;
                seqyA=seqy[j]+seqyA;
            }
        }
    }
    /* copy N terminal sequence */
    if (invmap_only!=1)
    {
        for (p=0;p<j;p++) buf[p]='-';
        strncpy(buf+j,seqx,i);
        buf[i+j]=0;
        seqxA=buf+seqxA;

        strncpy(buf,seqy,j);
        for (p=j;p<j+i;p++) buf[p]='-';
        buf[i+j]=0;
        seqyA=buf+seqyA;
    }
    delete [] buf;
}

/* entry function for NWalign
 * invmap_only - whether to return seqxA and seqyA or to return invmap
 *               0: only return seqxA and seqyA
 *               1: only return invmap
 *               2: return seqxA, seqyA and invmap */
int NWalign_main(const char *seqx, const char *seqy, const int xlen,
    const int ylen, string & seqxA, string & seqyA, const int mol_type,
    int *invmap, const int invmap_only=0, const int glocal=0)
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
        if (glocal==3)
        {
            gapopen=-5;
            gapext =-2;
        }
    }

    for (i=0;i<xlen+1;i++)
    {
        for (j=0;j<ylen+1;j++)
        {
            if (i*j==0) S[i][j]=0;
            else S[i][j]=BLOSUM[seqx[i-1]][seqy[j-1]];
        }
    }

    aln_score=calculate_score_gotoh(xlen, ylen, S, JumpH, JumpV, P,
        gapopen, gapext, glocal);

    seqxA.clear();
    seqyA.clear();

    if (glocal<3) trace_back_gotoh(seqx, seqy, JumpH, JumpV, P,
            seqxA, seqyA, xlen, ylen, invmap, invmap_only);
    else trace_back_sw(seqx, seqy, JumpH, JumpV, P, seqxA, seqyA,
            xlen, ylen, invmap, invmap_only);

    DeleteArray(&JumpH, xlen+1);
    DeleteArray(&JumpV, xlen+1);
    DeleteArray(&P, xlen+1);
    DeleteArray(&S, xlen+1);
    return aln_score; // aligment score
}

void get_seqID(int *invmap, const char *seqx, const char *seqy, 
    const int ylen, double &Liden,int &L_ali)
{
    Liden=0;
    L_ali=0;
    int i,j;
    for (j=0;j<ylen;j++)
    {
        i=invmap[j];
        if (i<0) continue;
        L_ali+=1;
        Liden+=(seqx[i]==seqy[j]);
    }
    //return L_ali?1.*Liden/L_ali:0;
}

void get_seqID(const string& seqxA, const string& seqyA,
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
    //return L_ali?1.*Liden/L_ali:0;
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
        printf("#score=%d\tLali=%d\tseqID_ali=%.3f\n", aln_score, L_ali, Liden/L_ali);
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

/* extract pairwise sequence alignment from residue index vectors,
 * assuming that "sequence" contains two empty strings.
 * return length of alignment, including gap. */
int extract_aln_from_resi(vector<string> &sequence, char *seqx, char *seqy,
    const vector<string> resi_vec1, const vector<string> resi_vec2,
    const int byresi_opt)
{
    sequence.clear();
    sequence.push_back("");
    sequence.push_back("");

    int i1=0; // positions in resi_vec1
    int i2=0; // positions in resi_vec2
    int xlen=resi_vec1.size();
    int ylen=resi_vec2.size();
    if (byresi_opt==4 || byresi_opt==5) // global or glocal sequence alignment
    {
        int *invmap;
        int glocal=0;
        if (byresi_opt==5) glocal=2;
        int mol_type=0;
        for (i1=0;i1<xlen;i1++)
            if ('a'<seqx[i1] && seqx[i1]<'z') mol_type++;
            else mol_type--;
        for (i2=0;i2<ylen;i2++)
            if ('a'<seqx[i2] && seqx[i2]<'z') mol_type++;
            else mol_type--;
        NWalign_main(seqx, seqy, xlen, ylen, sequence[0],sequence[1],
            mol_type, invmap, 0, glocal);
    }


    map<string,string> chainID_map1;
    map<string,string> chainID_map2;
    if (byresi_opt==3)
    {
        vector<string> chainID_vec;
        string chainID;
        stringstream ss;
        int i;
        for (i=0;i<xlen;i++)
        {
            chainID=resi_vec1[i].substr(5);
            if (!chainID_vec.size()|| chainID_vec.back()!=chainID)
            {
                chainID_vec.push_back(chainID);
                ss<<chainID_vec.size();
                chainID_map1[chainID]=ss.str();
                ss.str("");
            }
        }
        chainID_vec.clear();
        for (i=0;i<ylen;i++)
        {
            chainID=resi_vec2[i].substr(5);
            if (!chainID_vec.size()|| chainID_vec.back()!=chainID)
            {
                chainID_vec.push_back(chainID);
                ss<<chainID_vec.size();
                chainID_map2[chainID]=ss.str();
                ss.str("");
            }
        }
        vector<string>().swap(chainID_vec);
    }
    string chainID1="";
    string chainID2="";
    string chainID1_prev="";
    string chainID2_prev="";
    while(i1<xlen && i2<ylen)
    {
        if (byresi_opt==2)
        {
            chainID1=resi_vec1[i1].substr(5);
            chainID2=resi_vec2[i2].substr(5);
        }
        else if (byresi_opt==3)
        {
            chainID1=chainID_map1[resi_vec1[i1].substr(5)];
            chainID2=chainID_map2[resi_vec2[i2].substr(5)];
        }

        if (chainID1==chainID2)
        {
            if (atoi(resi_vec1[i1].substr(0,4).c_str())<
                atoi(resi_vec2[i2].substr(0,4).c_str()))
            {
                sequence[0]+=seqx[i1++];
                sequence[1]+='-';
            }
            else if (atoi(resi_vec1[i1].substr(0,4).c_str())>
                     atoi(resi_vec2[i2].substr(0,4).c_str()))
            {
                sequence[0]+='-';
                sequence[1]+=seqy[i2++];
            }
            else
            {
                sequence[0]+=seqx[i1++];
                sequence[1]+=seqy[i2++];
            }
            chainID1_prev=chainID1;
            chainID2_prev=chainID2;
        }
        else
        {
            if (chainID1_prev==chainID1 && chainID2_prev!=chainID2)
            {
                sequence[0]+=seqx[i1++];
                sequence[1]+='-';
                chainID1_prev=chainID1;
            }
            else if (chainID1_prev!=chainID1 && chainID2_prev==chainID2)
            {
                sequence[0]+='-';
                sequence[1]+=seqy[i2++];
                chainID2_prev=chainID2;
            }
            else
            {
                sequence[0]+=seqx[i1++];
                sequence[1]+=seqy[i2++];
                chainID1_prev=chainID1;
                chainID2_prev=chainID2;
            }
        }
        
    }
    map<string,string>().swap(chainID_map1);
    map<string,string>().swap(chainID_map2);
    chainID1.clear();
    chainID2.clear();
    chainID1_prev.clear();
    chainID2_prev.clear();
    return sequence[0].size();
}

/* extract pairwise sequence alignment from residue index vectors,
 * return length of alignment, including gap. */
int extract_aln_from_resi(vector<string> &sequence, char *seqx, char *seqy,
    const vector<string> resi_vec1, const vector<string> resi_vec2,
    const vector<int> xlen_vec, const vector<int> ylen_vec,
    const int chain_i, const int chain_j)
{
    sequence.clear();
    sequence.push_back("");
    sequence.push_back("");

    int i1=0; // positions in resi_vec1
    int i2=0; // positions in resi_vec2
    int xlen=xlen_vec[chain_i];
    int ylen=ylen_vec[chain_j];
    int i,j;
    for (i=0;i<chain_i;i++) i1+=xlen_vec[i];
    for (j=0;j<chain_j;j++) i2+=ylen_vec[j];

    i=j=0;
    while(i<xlen && j<ylen)
    {
        if (atoi(resi_vec1[i+i1].substr(0,4).c_str())<
            atoi(resi_vec2[j+i2].substr(0,4).c_str()))
        {
            sequence[0]+=seqx[i++];
            sequence[1]+='-';
        }
        else if (atoi(resi_vec1[i+i1].substr(0,4).c_str())>
                 atoi(resi_vec2[j+i2].substr(0,4).c_str()))
        {
            sequence[0]+='-';
            sequence[1]+=seqy[j++];
        }
        else
        {
            sequence[0]+=seqx[i++];
            sequence[1]+=seqy[j++];
        }
    }
    if (i<xlen && j==ylen)
    {
        for (i;i<xlen;i++)
        {
            sequence[0]+=seqx[i];
            sequence[1]+='-';
        }
    }
    else if (i==xlen && j<ylen)
    {
        for (j;j<ylen;j++)
        {
            sequence[0]+='-';
            sequence[1]+=seqy[j];
        }
    }
    return sequence[0].size();
}

#endif
