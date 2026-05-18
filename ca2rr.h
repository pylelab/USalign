#ifndef TMalign_ca2ct_h
#define TMalign_ca2ct_h 1

#include "TMalign.h"
#include "GeometryTools.h"

using namespace std;

const double CACB=1.53;

void ca2cb(double **xa, char *seqx, char *secx, int xlen, double **ya)
{
    int r;
    for (r=0;r<xlen;r++)
    {
        ya[r][0]=xa[r][0];
        ya[r][1]=xa[r][1];
        ya[r][2]=xa[r][2];
    }
    if (xlen<=2) return;

    double v1[3],v2[3],v3[3];
    for (r=0;r<xlen;r++)
    {
        if (seqx[r]=='G') continue; // no CB for glycine

        if (r==0) subtract(xa[2], xa[1],   v1);
        else      subtract(xa[r], xa[r-1], v1);
        if (!norm(v1)) continue;
        if (r+1==xlen) subtract(xa[r-2], xa[r-1], v2);
        else           subtract(xa[r],   xa[r+1], v2);
        if (!norm(v2)) continue;
        vectorsum(v1,v2,v3);
        if (!norm(v3)) continue;
        ya[r][0]+=v3[0]*CACB;
        ya[r][1]+=v3[1]*CACB;
        ya[r][2]+=v3[2]*CACB;

        /* rotate by 20.9 ~ 23.2 ~ 51.2 */
        if (r+1==xlen) subtract(xa[xlen-1], xa[xlen-3], v2);
        else if (r==0) subtract(xa[2],      xa[0],      v2);
        else           subtract(xa[r+1],    xa[r-1],    v2);
        vectorsum(xa[r],v2,v3);
            
        double angle=35;
        switch (secx[r])
        {
            case 'H': angle=51.2; break; // helix
            case 'E': angle=23.2; break; // parallel/antiparallel 20.9/23.2
            case 'S': angle=23.2; break; // parallel/antiparallel 20.9/23.2
        }

        CoordinateRotation(ya[r], v3,xa[r], angle,v1);
        ya[r][0]=v1[0];
        ya[r][1]=v1[1];
        ya[r][2]=v1[2];
    }
    return;
}

int calFUscore(bool **ct,int d_start, int d_end,double &minFUscore,
    const int mdl_opt, const bool verbose=false)
{
    int minl=0;
    int    l;
    long   N1,N2,N12;
    double FUscore;
    int    r1,r2;
    if (verbose) cout<<"#(domain1,domain2)\tN1\tN2\tN12\tFUscore"<<endl;
    l=1+d_start;
    N1=1+ct[d_start][d_start];
    N2=1;
    N12=1;
    for (r2=1+d_start;r2<d_end;r2++) N12+=ct[d_start][r2];
    for (r1=l;r1<d_end;r1++)
        for (r2=l;r2<d_end;r2++)
            N2+=ct[r1][r2];
    for (r2=1+d_start;r2<d_end;r2++) N12+=ct[d_start][r2];
    FUscore=2*N12*(1./N1+1./N2);
    if (verbose) cout<<'('<<1+d_start<<','<<1+d_start
        <<")("<<2+d_start<<','<<d_end<<")\t"<<N1<<'\t'<<N2<<'\t'
        <<N12<<'\t'<<FUscore<<endl;
    minFUscore=(N1+N2+2*N12);
    for (l=2+d_start;l<d_end;l++)
    {
        N1+=ct[l-1][l-1];
        for (r1=d_start;r1<l-1;r1++) N1+=2*ct[r1][l-1];
        N2-=ct[l-1][l-1];
        for (r2=l;r2<d_end;r2++) N2-=2*ct[l-1][r2];
        for (r1=d_start;r1<l-1;r1++)  N12-=ct[r1][l-1];
        for (r2=l;r2<d_end;r2++) N12+=ct[l-1][r2];
        FUscore=2*N12*(1./N1+1./N2);
        if (verbose) cout<<'('<<1+d_start<<','<<l<<")("<<l+1<<','
            <<d_end<<")\t"<<N1<<'\t'<<N2<<'\t'
            <<N12<<'\t'<<FUscore<<endl;
        if ((d_end-l)<mdl_opt) break;
        if ((l-d_start)>=mdl_opt && FUscore<=minFUscore)
        {
            minl=l;
            minFUscore=FUscore;
        }
    }
    //cout<<"minl="<<minl<<" minFUscore="<<minFUscore<<endl;
    return minl;
}


void iterative_calFUscore(bool **ct, vector<int> &l_vec, const int xlen,
    const int hinge_opt, const int mdl_opt, const bool verbose=true)
{
    int    d_start   =0;
    int    d_end     =xlen;
    double minFUscore=0;
    int minl=calFUscore(ct,d_start,d_end,minFUscore,mdl_opt);
    if (minl<=0) return;
    if (verbose) cout<<"#partition 1: "<<'('<<1+d_start<<','
        <<minl<<")("<<minl+1<<',' <<d_end<<") FUscore="<<minFUscore<<endl;
    l_vec.push_back(minl);
    if (hinge_opt<=2) return;
    vector<pair<double,int> > FU_l_vec;
    for (int iter=2;iter<hinge_opt;iter++)
    {
        for (int r=0;r<=l_vec.size();r++)
        {
            d_start=(r==0)?0:l_vec[r-1];
            d_end  =(r==l_vec.size())?xlen:l_vec[r];
            if (d_end-d_start<2*mdl_opt) continue;
            minl=calFUscore(ct,d_start,d_end,minFUscore,mdl_opt);
            if (minl<=0) continue;
            FU_l_vec.push_back(make_pair(minFUscore,minl));
        }
        if (FU_l_vec.size()==0) break;

        sort(FU_l_vec.begin(),FU_l_vec.end());
        l_vec.push_back(FU_l_vec[0].second);
        sort(l_vec.begin(),l_vec.end());
        if (verbose) cout<<"#partition "<<iter<<": ";
        for (int r=0;r<=l_vec.size();r++)
        {
            d_start=(r==0)?0:l_vec[r-1];
            d_end  =(r==l_vec.size())?xlen:l_vec[r];
            if (verbose) cout<<'('<<1+d_start<<','<<d_end<<')';
        }
        if (verbose) cout<<" FUscore="<<FU_l_vec[0].first<<endl;
        FU_l_vec.clear();
    }
    return;
}

#endif
