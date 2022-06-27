/* Partial implementation of Needleman-Wunsch (NW) dynamic programming for
 * global alignment. The three NWDP_TM functions below are not complete
 * implementation of NW algorithm because gap jumping in the standard Gotoh
 * algorithm is not considered. Since the gap opening and gap extension is
 * the same, this is not a problem. This code was exploited in TM-align
 * because it is about 1.5 times faster than a complete NW implementation.
 * Nevertheless, if gap opening != gap extension shall be implemented in
 * the future, the Gotoh algorithm must be implemented. In rare scenarios,
 * it is also possible to have asymmetric alignment (i.e. 
 * TMalign A.pdb B.pdb and TMalign B.pdb A.pdb have different TM_A and TM_B
 * values) caused by the NWPD_TM implement.
 */

/* Input: score[1:len1, 1:len2], and gap_open
 * Output: j2i[1:len2] \in {1:len1} U {-1}
 * path[0:len1, 0:len2]=1,2,3, from diagonal, horizontal, vertical */
void NWDP_TM(double **score, bool **path, double **val,
    int len1, int len2, double gap_open, int j2i[])
{

    int i, j;
    double h, v, d;

    //initialization
    for(i=0; i<=len1; i++)
    {
        val[i][0]=0;
        //val[i][0]=i*gap_open;
        path[i][0]=false; //not from diagonal
    }

    for(j=0; j<=len2; j++)
    {
        val[0][j]=0;
        //val[0][j]=j*gap_open;
        path[0][j]=false; //not from diagonal
        j2i[j]=-1;    //all are not aligned, only use j2i[1:len2]
    }      


    //decide matrix and path
    for(i=1; i<=len1; i++)
    {
        for(j=1; j<=len2; j++)
        {
            d=val[i-1][j-1]+score[i][j]; //diagonal

            //symbol insertion in horizontal (= a gap in vertical)
            h=val[i-1][j];
            if(path[i-1][j]) h += gap_open; //aligned in last position

            //symbol insertion in vertical
            v=val[i][j-1];
            if(path[i][j-1]) v += gap_open; //aligned in last position


            if(d>=h && d>=v)
            {
                path[i][j]=true; //from diagonal
                val[i][j]=d;
            }
            else 
            {
                path[i][j]=false; //from horizontal
                if(v>=h) val[i][j]=v;
                else val[i][j]=h;
            }
        } //for i
    } //for j

    //trace back to extract the alignment
    i=len1;
    j=len2;
    while(i>0 && j>0)
    {
        if(path[i][j]) //from diagonal
        {
            j2i[j-1]=i-1;
            i--;
            j--;
        }
        else 
        {
            h=val[i-1][j];
            if(path[i-1][j]) h +=gap_open;

            v=val[i][j-1];
            if(path[i][j-1]) v +=gap_open;

            if(v>=h) j--;
            else i--;
        }
    }
}

/* Input: vectors x, y, rotation matrix t, u, scale factor d02, and gap_open
 * Output: j2i[1:len2] \in {1:len1} U {-1}
 * path[0:len1, 0:len2]=1,2,3, from diagonal, horizontal, vertical */
void NWDP_TM(bool **path, double **val, double **x, double **y,
    int len1, int len2, double t[3], double u[3][3],
    double d02, double gap_open, int j2i[])
{
    int i, j;
    double h, v, d;

    //initialization. use old val[i][0] and val[0][j] initialization
    //to minimize difference from TMalign fortran version
    for(i=0; i<=len1; i++)
    {
        val[i][0]=0;
        //val[i][0]=i*gap_open;
        path[i][0]=false; //not from diagonal
    }

    for(j=0; j<=len2; j++)
    {
        val[0][j]=0;
        //val[0][j]=j*gap_open;
        path[0][j]=false; //not from diagonal
        j2i[j]=-1;    //all are not aligned, only use j2i[1:len2]
    }      
    double xx[3], dij;


    //decide matrix and path
    for(i=1; i<=len1; i++)
    {
        transform(t, u, &x[i-1][0], xx);
        for(j=1; j<=len2; j++)
        {
            dij=dist(xx, &y[j-1][0]);    
            d=val[i-1][j-1] +  1.0/(1+dij/d02);

            //symbol insertion in horizontal (= a gap in vertical)
            h=val[i-1][j];
            if(path[i-1][j]) h += gap_open; //aligned in last position

            //symbol insertion in vertical
            v=val[i][j-1];
            if(path[i][j-1]) v += gap_open; //aligned in last position


            if(d>=h && d>=v)
            {
                path[i][j]=true; //from diagonal
                val[i][j]=d;
            }
            else 
            {
                path[i][j]=false; //from horizontal
                if(v>=h) val[i][j]=v;
                else val[i][j]=h;
            }
        } //for i
    } //for j

    //trace back to extract the alignment
    i=len1;
    j=len2;
    while(i>0 && j>0)
    {
        if(path[i][j]) //from diagonal
        {
            j2i[j-1]=i-1;
            i--;
            j--;
        }
        else 
        {
            h=val[i-1][j];
            if(path[i-1][j]) h +=gap_open;

            v=val[i][j-1];
            if(path[i][j-1]) v +=gap_open;

            if(v>=h) j--;
            else i--;
        }
    }
}

/* This is the same as the previous NWDP_TM, except for the lack of rotation
 * Input: vectors x, y, scale factor d02, and gap_open
 * Output: j2i[1:len2] \in {1:len1} U {-1}
 * path[0:len1, 0:len2]=1,2,3, from diagonal, horizontal, vertical */
void NWDP_SE(bool **path, double **val, double **x, double **y,
    int len1, int len2, double d02, double gap_open, int j2i[])
{
    int i, j;
    double h, v, d;

    for(i=0; i<=len1; i++)
    {
        val[i][0]=0;
        path[i][0]=false; //not from diagonal
    }

    for(j=0; j<=len2; j++)
    {
        val[0][j]=0;
        path[0][j]=false; //not from diagonal
        j2i[j]=-1;    //all are not aligned, only use j2i[1:len2]
    }      
    double dij;

    //decide matrix and path
    for(i=1; i<=len1; i++)
    {
        for(j=1; j<=len2; j++)
        {
            dij=dist(&x[i-1][0], &y[j-1][0]);    
            d=val[i-1][j-1] +  1.0/(1+dij/d02);

            //symbol insertion in horizontal (= a gap in vertical)
            h=val[i-1][j];
            if(path[i-1][j]) h += gap_open; //aligned in last position

            //symbol insertion in vertical
            v=val[i][j-1];
            if(path[i][j-1]) v += gap_open; //aligned in last position


            if(d>=h && d>=v)
            {
                path[i][j]=true; //from diagonal
                val[i][j]=d;
            }
            else 
            {
                path[i][j]=false; //from horizontal
                if(v>=h) val[i][j]=v;
                else val[i][j]=h;
            }
        } //for i
    } //for j

    //trace back to extract the alignment
    i=len1;
    j=len2;
    while(i>0 && j>0)
    {
        if(path[i][j]) //from diagonal
        {
            j2i[j-1]=i-1;
            i--;
            j--;
        }
        else 
        {
            h=val[i-1][j];
            if(path[i-1][j]) h +=gap_open;

            v=val[i][j-1];
            if(path[i][j-1]) v +=gap_open;

            if(v>=h) j--;
            else i--;
        }
    }
}

void NWDP_SE(bool **path, double **val, double **x, double **y,
    int len1, int len2, double d02, double gap_open, int j2i[],
    const int hinge)
{
    if (hinge==0)
    {
        NWDP_SE(path, val, x, y, len1, len2, d02, gap_open, j2i);
        return;
    }
    int i, j;
    double h, v, d;

    int L=(len2>len1)?len2:len1;
    int int_min=L*(gap_open-1);

    for (i=0; i<=len1; i++)
    {
        for (j=0; j<=len2; j++)
        {
            val[i][j]=0;
            path[i][j]=false;
        }
    }

    /* fill in old j2i */
    int k=0;
    for (j=0; j<len2; j++)
    {
        i=j2i[j];
        if (i<0) continue;
        path[i+1][j+1]=true;
        val[i+1][j+1]=0;
    }

    double dij;

    //decide matrix and path
    for(i=1; i<=len1; i++)
    {
        for(j=1; j<=len2; j++)
        {
            dij=0;
            if (path[i][j]==false) dij=dist(&x[i-1][0], &y[j-1][0]);    
            d=val[i-1][j-1] +  1.0/(1+dij/d02);

            //symbol insertion in horizontal (= a gap in vertical)
            h=val[i-1][j];
            if(path[i-1][j]) h += gap_open; //aligned in last position

            //symbol insertion in vertical
            v=val[i][j-1];
            if(path[i][j-1]) v += gap_open; //aligned in last position


            if(d>=h && d>=v && val[i][j]==0)
            {
                path[i][j]=true; //from diagonal
                val[i][j]=d;
            }
            else 
            {
                path[i][j]=false; //from horizontal
                if(v>=h) val[i][j]=v;
                else val[i][j]=h;
            }
        } //for i
    } //for j

    //trace back to extract the alignment
    for (j=0;j<=len2;j++) j2i[j]=-1;
    i=len1;
    j=len2;
    while(i>0 && j>0)
    {
        if(path[i][j]) //from diagonal
        {
            j2i[j-1]=i-1;
            i--;
            j--;
        }
        else 
        {
            h=val[i-1][j];
            if(path[i-1][j]) h +=gap_open;

            v=val[i][j-1];
            if(path[i][j-1]) v +=gap_open;

            if(v>=h) j--;
            else i--;
        }
    }
}

/* +ss
 * Input: secondary structure secx, secy, and gap_open
 * Output: j2i[1:len2] \in {1:len1} U {-1}
 * path[0:len1, 0:len2]=1,2,3, from diagonal, horizontal, vertical */
void NWDP_TM(bool **path, double **val, const char *secx, const char *secy,
    const int len1, const int len2, const double gap_open, int j2i[])
{

    int i, j;
    double h, v, d;

    //initialization
    for(i=0; i<=len1; i++)
    {
        val[i][0]=0;
        //val[i][0]=i*gap_open;
        path[i][0]=false; //not from diagonal
    }

    for(j=0; j<=len2; j++)
    {
        val[0][j]=0;
        //val[0][j]=j*gap_open;
        path[0][j]=false; //not from diagonal
        j2i[j]=-1;    //all are not aligned, only use j2i[1:len2]
    }      

    //decide matrix and path
    for(i=1; i<=len1; i++)
    {
        for(j=1; j<=len2; j++)
        {
            d=val[i-1][j-1] + 1.0*(secx[i-1]==secy[j-1]);

            //symbol insertion in horizontal (= a gap in vertical)
            h=val[i-1][j];
            if(path[i-1][j]) h += gap_open; //aligned in last position

            //symbol insertion in vertical
            v=val[i][j-1];
            if(path[i][j-1]) v += gap_open; //aligned in last position

            if(d>=h && d>=v)
            {
                path[i][j]=true; //from diagonal
                val[i][j]=d;
            }
            else 
            {
                path[i][j]=false; //from horizontal
                if(v>=h) val[i][j]=v;
                else val[i][j]=h;
            }
        } //for i
    } //for j

    //trace back to extract the alignment
    i=len1;
    j=len2;
    while(i>0 && j>0)
    {
        if(path[i][j]) //from diagonal
        {
            j2i[j-1]=i-1;
            i--;
            j--;
        }
        else 
        {
            h=val[i-1][j];
            if(path[i-1][j]) h +=gap_open;

            v=val[i][j-1];
            if(path[i][j-1]) v +=gap_open;

            if(v>=h) j--;
            else i--;
        }
    }
}
