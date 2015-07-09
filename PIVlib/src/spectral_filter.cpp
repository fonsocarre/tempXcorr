#include "spectral_filter.h"

void spectral_filter::initialise(double cutLength)
{
    this->cutLength = cutLength;
}

void spectral_filter::filter(std::vector<double> in,
                             std::vector<double>& out,
                             int nRows, int nCols,
                             double dx, int padding)
{
    const int n = nRows*nCols;
    int nnCols = nCols + 2*padding;
    int nnRows = nRows + 2*padding;

    fftw_complex* tmatrix;
    tmatrix = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*n);

    // initialisation
    for (int i=0; i<n; ++i)
    {
        tmatrix[i][0] = in[i];
        tmatrix[i][1] = (double)0;
    }
    // padding
    fftw_complex* matrix;
    matrix = zeroPadding(tmatrix, nRows, nCols, padding);
    fftw_free(tmatrix);

    fftw_plan p;
    p = fftw_plan_dft_2d(nRows+2*padding, nCols+2*padding, matrix, matrix,
                                                       FFTW_FORWARD,
                                                       FFTW_ESTIMATE);


    //std::cout << "first exec= " << __LINE__<< std::flush << std::endl;
    fftw_execute(p);

    // apply filter
    //fftw_complex* window;
    //window = spectralHanningWindow(nRows, nCols);
    //for (int i=0; i<n; ++i)
    //{
        //matrix[i][0] *= window[i][0];
        //matrix[i][1] *= window[i][1];

        //matrix[i][0] = window[i][0];
        //matrix[i][1] = window[i][1];
    //}
    //fftw_free(window);

    // fx is the frequencies in a row
    double xlength = nnCols*dx;
    double xds = xlength;
    std::vector<double> fx(nnCols/2 + 1, 0.0);
    for (int i=0; i<nnCols/2+1; ++i)
    {
        fx[i] = xds/2. * static_cast<double>(i)/(nnCols/2.0);
    }

    // fy is the frequencies in a column
    double ylength = nnRows*dx;
    double yds = ylength;
    std::vector<double> fy(nnRows/2 + 1, 0.0);
    for (int i=0; i<nnRows/2+1; ++i)
    {
        fy[i] = yds/2. * static_cast<double>(i)/(nnRows/2.0);
    }

    // TODO change filtering, it is wrong!
    for (int i=0; i<nnRows/2 + 1; ++i)
    {
        int ii = i+padding;
        for (int j=0; j<nnCols/2 + 1; ++j)
        {
            int jj = j+padding;
            if (fx[j] >= xlength/this->cutLength || 
                fy[i] >= ylength/this->cutLength)
            {
                //std::cout << i << "  " << j << std::endl << std::flush;
                assert(j+i*nnCols < nnCols*nnRows);
                matrix[j + i*nnCols][0] = 0.0;
                matrix[j + i*nnCols][1] = 0.0;

                assert(nnCols-1-j + (nnCols-1)*(nnRows-1-i) < nnCols*nnRows);
                matrix[nnCols-1-j + (nnCols-1)*(nnRows-1-i)][0] = 0.0;
                matrix[nnCols-1-j + (nnCols-1)*(nnRows-1-i)][1] = 0.0;

                assert(j + (nnCols-1)*(nnRows-i-1) < nnCols*nnRows);
                matrix[j + (nnCols-1)*(nnRows-i-1)][0] = 0.0;
                matrix[j + (nnCols-1)*(nnRows-i-1)][1] = 0.0;

                assert(nnCols-1-j + i*nnCols < nnCols*nnRows);
                matrix[nnCols-1-j + i*nnCols][0] = 0.0;
                matrix[nnCols-1-j + i*nnCols][1] = 0.0;
            }
        }
    }
    
    //std::cout << "= " << __LINE__<< std::flush << std::endl;
    p = fftw_plan_dft_2d(nRows+2*padding, nCols+2*padding, matrix, matrix,
                                                       FFTW_BACKWARD,
                                                       FFTW_ESTIMATE);

    fftw_execute(p);
    //std::cout << "all fft executed" << std::endl;
    
    fftw_complex* nonPadMatrix;
    nonPadMatrix = removePadding(matrix, nRows, nCols, padding);

    out.resize(n);
    for (int i=0; i<n; ++i)
    {
        out[i] = nonPadMatrix[i][0]/(nnRows*nnCols);
    }
    fftw_free(nonPadMatrix);
    fftw_free(matrix);
}

