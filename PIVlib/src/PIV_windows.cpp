#include "PIV_windows.h"

void getHanningWindow(std::vector<double> initial,
                        std::vector<double>& modified,
                        int nRows, int nCols)
{
    const double alpha = 25.0/(double)46.0;
    const double beta  = 21.0/(double)46.0;

    for (int i=0; i<nRows; ++i)
    {
        for (int j=0; j<nCols; ++j)
        {
            modified[i*nCols + j] = initial[i*nCols + j] *
                (alpha - beta*cos(2.*M_PI*i/(nRows-1.))) *
                (alpha - beta*cos(2.*M_PI*j/(nCols-1.)));

        }
    }
}

fftw_complex* spectralHanningWindow(int nRows, int nCols)
{
    int width = 40;
    int n = nRows*nCols;
    //std::vector<double> output(n, 0.0);
    fftw_complex* output;
    fftw_complex* input;

    input =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
    output = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);

    const double alpha = 25.0/(double)46.0;
    const double beta  = 21.0/(double)46.0;

    for (int i=0; i<nRows; ++i)
    {
        for (int j=0; j<nCols; ++j)
        {
            input[i*nCols + j][0] = (alpha - beta*cos(2.*M_PI*i/(nRows-1.))) *
                                     (alpha - beta*cos(2.*M_PI*j/(nCols-1.)));
            input[i*nCols + j][1] = 0.0;
        }
    }
    
    fftw_plan p;
    p = fftw_plan_dft_2d(nRows, nCols, input, output,
                            FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_free(input);
    //for (int i=0; i<n; ++i)
    //{
        //output[i][0] /= sqrt(n);
        //output[i][1] /= sqrt(n);
    //}

    fftw_destroy_plan(p);
    return output;
}

fftw_complex* zeroPadding(fftw_complex* input, int nRows, int nCols,
                                int padding)
{
    int nnRows = nRows + 2*padding;
    int nnCols = nCols + 2*padding;

    fftw_complex* output;
    output = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nnRows*nnCols);

    for (int i=0; i<nnRows*nnCols; ++i)
    {
        output[i][0] = (double)0.0;
        output[i][1] = (double)0.0;
    }
    
    for (int i=0; i<nRows; ++i)
    {
        int ni = i+padding;
        for (int j=0; j<nCols; ++j)
        {
            int nj = j+padding;
            output[ni*nnCols + nj][0] = input[i*nCols + j][0];
            output[ni*nnCols + nj][1] = input[i*nCols + j][1];
        }
    }
    return output;
}

fftw_complex* removePadding(fftw_complex* input, int nRows, int nCols,
                                int padding)
{
    int nnRows = nRows + 2*padding;
    int nnCols = nCols + 2*padding;

    fftw_complex* output;
    output = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nCols*nRows);

    for (int i=0; i<nRows; ++i)
    {
        int ni = i+padding;
        for (int j=0; j<nCols; ++j)
        {
            int nj = j+padding; 
            output[i*nCols + j][0] = input[ni*nnCols + nj][0];
            output[i*nCols + j][1] = input[ni*nnCols + nj][1];
        }
    }
    return output;
}










