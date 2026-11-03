/* Micro-benchmark of Intel MKL kernels matching crates/dft-numerics/examples/lib_bench.rs
 * (original code written for this project; MKL is only linked, not copied).
 * Build (sequential MKL; MKLROOT from the QE environment):
 *   gcc -O2 -I$MKLROOT/include tools/lib_bench/mkl_bench.c -o mkl_bench \
 *       -L$MKLROOT/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
 * Run with one thread: MKL_NUM_THREADS=1 taskset -c 6 ./mkl_bench
 */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "mkl.h"
#include "mkl_dfti.h"
#include "mkl_lapacke.h"

static unsigned long long s = 12345ULL;
static double rnd(void) { s = s * 6364136223846793005ULL + 1442695040888963407ULL; return (double)(s >> 11) / 9007199254740992.0 - 0.5; }
static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec * 1e6 + t.tv_nsec * 1e-3; }

int main(void) {
    printf("kind size microseconds_per_op   (op = one complex 3D transform, forward or inverse, in place; one gemm; one eigendecomposition)\n");
    int shapes[5] = {24, 25, 27, 36, 48};
    for (int q = 0; q < 5; q++) {
        MKL_LONG n = shapes[q], dims[3] = {n, n, n};
        size_t N = (size_t)n * n * n;
        MKL_Complex16 *a = malloc(N * sizeof *a);
        for (size_t i = 0; i < N; i++) { a[i].real = rnd(); a[i].imag = rnd(); }
        DFTI_DESCRIPTOR_HANDLE h;
        DftiCreateDescriptor(&h, DFTI_DOUBLE, DFTI_COMPLEX, 3, dims);
        DftiSetValue(h, DFTI_PLACEMENT, DFTI_INPLACE);
        DftiCommitDescriptor(h);
        int reps = 3000000 / (int)N; if (reps < 20) reps = 20;
        for (int i = 0; i < 3; i++) { DftiComputeForward(h, a); DftiComputeBackward(h, a); }
        double t0 = now();
        for (int i = 0; i < reps; i++) { DftiComputeForward(h, a); DftiComputeBackward(h, a); }
        printf("fft3d %ld^3 %.2f\n", (long)n, (now() - t0) / (2.0 * reps));
        DftiFreeDescriptor(&h); free(a);
    }
    int npw = 750, ms[4] = {8, 16, 32, 48};
    for (int q = 0; q < 4; q++) {
        int m = ms[q], k = 8;
        MKL_Complex16 *A = malloc((size_t)npw * m * sizeof *A), *B = malloc((size_t)npw * k * sizeof *B), *C = malloc((size_t)m * k * sizeof *C);
        for (size_t i = 0; i < (size_t)npw * m; i++) { A[i].real = rnd(); A[i].imag = rnd(); }
        for (size_t i = 0; i < (size_t)npw * k; i++) { B[i].real = rnd(); B[i].imag = rnd(); }
        MKL_Complex16 one = {1, 0}, zero = {0, 0};
        int reps = 2000;
        for (int i = 0; i < 10; i++) cblas_zgemm(CblasColMajor, CblasConjTrans, CblasNoTrans, m, k, npw, &one, A, npw, B, npw, &zero, C, m);
        double t0 = now();
        for (int i = 0; i < reps; i++) cblas_zgemm(CblasColMajor, CblasConjTrans, CblasNoTrans, m, k, npw, &one, A, npw, B, npw, &zero, C, m);
        printf("zgemm_CN %dx%d^H * %dx8 %.2f\n", npw, m, npw, (now() - t0) / reps);
        free(A); free(B); free(C);
    }
    int es[3] = {16, 32, 48};
    for (int q = 0; q < 3; q++) {
        int m = es[q];
        MKL_Complex16 *H = malloc((size_t)m * m * sizeof *H), *W = malloc((size_t)m * m * sizeof *W);
        double *ev = malloc(m * sizeof *ev);
        for (int i = 0; i < m; i++) for (int j = 0; j <= i; j++) {
            double re = rnd(), im = (i == j) ? 0.0 : rnd();
            H[i + (size_t)j * m].real = re; H[i + (size_t)j * m].imag = im;
            H[j + (size_t)i * m].real = re; H[j + (size_t)i * m].imag = -im;
        }
        int reps = 300;
        for (int i = 0; i < 3; i++) { for (size_t t = 0; t < (size_t)m * m; t++) W[t] = H[t]; LAPACKE_zheev(LAPACK_COL_MAJOR, 'V', 'U', m, W, m, ev); }
        double t0 = now();
        for (int i = 0; i < reps; i++) { for (size_t t = 0; t < (size_t)m * m; t++) W[t] = H[t]; LAPACKE_zheev(LAPACK_COL_MAJOR, 'V', 'U', m, W, m, ev); }
        printf("zheev %dx%d %.2f\n", m, m, (now() - t0) / reps);
        free(H); free(W); free(ev);
    }
    return 0;
}
