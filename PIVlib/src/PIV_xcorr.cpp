#include "PIV_xcorr.h"

void PIV::xxcorr(const std::vector<double>& mat0,
           const std::vector<double>& mat1,
           const int wNx,
           const int wNy,
           int nRows,
           int nCols,
           std::vector<double>& xdisps,
           std::vector<double>& ydisps,
           std::vector<double>& xcorrvals)
{
    int n = static_cast<int>(mat0.size());
    xdisps.resize(n);
    ydisps.resize(n);
    xcorrvals.resize(n);


    // UGLY PARAMETER TODO
    double lim = 10.0;

    double xdisp = 0.0;
    double ydisp = 0.0;
    double xcorrval = 0.0;
    for (int i=0; i<nRows; ++i)
    {
        int imin = ((i-wNy < 0)? 0:i-wNy);
        int imax = ((i+wNy >= nRows)? nRows-1:i+wNy);
        for (int j=0; j<nCols; ++j)
        {
            int jmin = ((j-wNx < 0)? 0:j-wNx);
            int jmax = ((j+wNx >= nCols)? nCols-1:j+wNx);

            int nnCols = jmax-jmin+1;
            int nnRows = imax-imin+1;
            //std::vector<double> txdisps(nnCols*nnRows, 0.0);
            //std::vector<double> tydisps(nnCols*nnRows, 0.0);
            //std::vector<double> txcorrvals(nnCols*nnRows, 0.0);
            xdisp = 0.0;
            ydisp = 0.0;
            xcorrval = 0.0;
            double maxxcorrval = -1.0e10;
            for (int ii=imin; ii<=imax; ++ii)
            {
                for (int jj=jmin; jj<=jmax; ++jj)
                {
                    if (abs(mat1[ij(ii,jj,nRows, nCols)]) > lim) continue; 
                    xcorrval += mat0[ij(i,j,nRows,nCols)]*
                                mat1[ij(ii, jj, nRows, nCols)];
                    if (xcorrval > maxxcorrval)
                    {
                        maxxcorrval = xcorrval;
                        xdisp = jj;
                        ydisp = ii;
                    }
                }
            }

            xdisps[ij(i,j,nRows,nCols)] = j-xdisp;
            ydisps[ij(i,j,nRows,nCols)] = i-ydisp;
            xcorrvals[ij(i,j,nRows,nCols)] = maxxcorrval;
        }
    }

}
