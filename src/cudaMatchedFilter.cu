#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cufft.h>

typedef float2  CComplex;
typedef double2 ZComplex;
template<typename T>
static __device__ __host__ inline T operator+(const T a, const T b);
template<typename T>
static __device__ __host__ inline T operator*(const T a, const T b);

// Complex addition
template<typename T> static __device__ __host__ 
inline T operator+(const T a, const T b)
{
    T c;
    c.x = a.x + b.x;
    c.y = a.y + b.y;
    return c;
}
// Complex multiplication
template<typename T> static __device__ __host__
inline T operator*(const T a, const T b)
{
    T c;
    c.x = a.x*b.x - a.y*b.y;
    c.y = a.x*b.y + a.y*b.x;
    return c;
}

/// @brief Computes: \f$ X \leftarrow B X  \f$.
/// @param[in,out] x    On input this is the spectra of the waveform.
///                     On exit, this is the spectra multiplied with the
///                     spectra of the filter - i.e., the convolution.
///                     This is an array whose dimension is [nw].
/// @param[in] b        The spectra of the filter coefficients.
///                     This is an array whose dimension is [nw].
/// @param[in] nw       The number of frequencies. 
static __global__
void multiplySpectra(CComplex *x, const CComplex *b, const int nw)
{
    const int numThreads = blockDim.x*gridDim.x;
    const int threadID   = blockIdx.x*blockDim.x + threadIdx.x;
    for (int i=threadID; i<nw; i=i+numThreads)
    {
        x[i] = x[i]*b[i];
    }
}


struct cufft32z_struct
{
    cufftHandle mPlan;
};

extern "C"
void clear(struct cufft32z_struct *cuft);

extern "C"
void initialize(struct cufft32z_struct *cuft);


void clear(struct cufft32z_struct *cuft)
{
    cufftDestroy(cuft->mPlan);
}

void initialize(struct cufft32z_struct *cuft)
{
    //checkCudaErrors(cufftMakePlan1d(plan_input, new_size, CUFFT_C2C, 1, worksize));  
}

