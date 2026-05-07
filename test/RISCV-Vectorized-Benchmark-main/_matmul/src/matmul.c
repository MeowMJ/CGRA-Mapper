/*************************************************************************
* Vectorized matrixmul Kernel
*************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdbool.h>

#define DATA_TYPE
typedef double data_t;
#define TILE_SIZE_K 32   // k维度分块大小

#ifdef USE_RISCV_VECTOR
#include <riscv_vector.h>
#include "../../common/vector_defines.h"

void matrixmul_intrinsics(data_t *a, data_t *b, data_t *c, int n, int m, int p) {

    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < n; j++) {
            size_t gvl = _MMR_VSETVL_E64M1(p);
            vfloat64m1_t vprod = _MM_SET_f64(0, gvl);
            vfloat64m1_t vsum  = _MM_SET_f64(0, gvl);

            for (size_t k = 0; k < p; k += gvl){
                gvl = _MMR_VSETVL_E64M1(p - k);

                // Matrix A row
                vfloat64m1_t va  = _MM_LOAD_f64(&a[i*p+k], gvl);
                // Matrix B column
                vfloat64m1_t vb = _MM_LOAD_STRIDE_f64(&b[k*n+j], n * sizeof(data_t), gvl);

                // A[0]*B[0], A[1]*B[1],... A[n]*B[n]
                vprod  = _MM_MACC_f64(vprod,va, vb, gvl);

            }//k
            gvl = _MMR_VSETVL_E64M1(p);
            vsum   = _MM_REDSUM_f64(vprod,vsum, gvl);
            c[i*n+j] = _MM_VGETFIRST_f64(vsum);
        }//j
    }//i
}

#elif defined(USE_LOOP_TILING)  // 循环分块版本（基于RISC-V向量指令）
#include <riscv_vector.h>
#include "../../common/vector_defines.h"

// matmul 最内部循环执行完毕一次会 load 矩阵 a 和 b 的 p 个数，即 16*p Byte。当 p <= 64 时，SPM <= 1KB.
void matrixmul_tiling(data_t *a, data_t *b, data_t *c, int n, int m, int p) {

    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < n; j++) {
            size_t gvl = _MMR_VSETVL_E64M1(p);
            vfloat64m1_t vprod = _MM_SET_f64(0, gvl);
            vfloat64m1_t vsum  = _MM_SET_f64(0, gvl);
            data_t tmp_c[(TILE_SIZE_K / gvl)];  // FIXME: TILE_SIZE_K / gvl
            for (int kk = 0; kk < p; kk += TILE_SIZE_K) {
                for (size_t k = kk; k < min(p, kk + TILE_SIZE_K); k += gvl){
                    gvl = _MMR_VSETVL_E64M1(p - k);

                    // Tile Load
                    // Matrix A row
                    vfloat64m1_t va  = _MM_LOAD_f64(&a[i*p+k], gvl);
                    // Matrix B column
                    vfloat64m1_t vb = _MM_LOAD_STRIDE_f64(&b[k*n+j], n * sizeof(data_t), gvl);

                    // Execution of A[0]*B[0], A[1]*B[1],... A[n]*B[n]
                    // vprod  = _MM_MACC_f64(vprod,va, vb, gvl);

                }//k
                gvl = _MMR_VSETVL_E64M1(p);
                // A[0]*B[0] + A[1]*B[1],... + A[n]*B[n]
                vsum   = _MM_REDSUM_f64(vprod,vsum, gvl);
                tmp_c[k%kk] = vsum;
            }//kk
            // final C, we need tile accumulation
            for (int t = 0; t < (TILE_SIZE_K / gvl); t++) {
                c[i*n+j] += tmp_c[t];
            }
        }//j
    }//i
}

#elif defined(USE_RISCV_VECTOR_SYSTOLIC)
#include <riscv_vector.h>
#include "../../common/vector_defines.h"

// matmul 最内部循环执行完毕一次会 load 矩阵 a 和 b 的 p 个数，即 16*p Byte。当 p <= 64 时，SPM <= 1KB.
void matrixmul_intrinsics_systolic(data_t *a, data_t *b, data_t *c, int n, int m, int p) {
    for (size_t i = 0; i < m; i++) {
        // 让 va 留在 systolic array 中，只当 c 换行时，即 i++ 时，此数据被读取，即每次读取的 i*p+k 都会存在 PE 的 regsiter 里
        // 那么，每一次读取，data reuse 的次数为 (p^2 - p)/2
        vfloat64m1_t va[p];
        for (size_t k = 0; k < p; k += gvl){
            vfloat64m1_t va[k] = _MM_LOAD_f64(&a[i*p+k], gvl);
        }
        for (size_t j = 0; j < n; j++) {
            size_t gvl = _MMR_VSETVL_E64M1(p);
            vfloat64m1_t vprod = _MM_SET_f64(0, gvl);
            vfloat64m1_t vsum  = _MM_SET_f64(0, gvl);
            data_t tmp_c[(TILE_SIZE_K / gvl)];
            for (int kk = 0; kk < p; kk += TILE_SIZE_K) {
                for (size_t k = kk; k < min(p, kk + TILE_SIZE_K); k += gvl){
                    gvl = _MMR_VSETVL_E64M1(p - k);

                    // Tile Load
                    // Matrix A row
                    // 只要第一个 C0j 读取了 va 的数据后，后面的所有 k 次循环都不需要再重新 load 了
                    vfloat64m1_t tmp_va = va[k];
                    // Matrix B column
                    vfloat64m1_t vb = _MM_LOAD_STRIDE_f64(&b[k*n+j], n * sizeof(data_t), gvl);

                    // Execution of A[0]*B[0], A[1]*B[1],... A[n]*B[n]
                    vprod  = _MM_MACC_f64(vprod,va, vb, gvl);

                }//k
                gvl = _MMR_VSETVL_E64M1(p);
                // A[0]*B[0] + A[1]*B[1],... + A[n]*B[n]
                vsum   = _MM_REDSUM_f64(vprod,vsum, gvl);
                tmp_c[k%kk] = vsum;
            }//kk
            // final C, we need tile accumulation
            for (int t = 0; t < (TILE_SIZE_K / gvl); t++) {
                c[i*n+j] += tmp_c[t]; // 按行访问
            }
        }//j
    }//i
}

#elif defined(USE_LOOP_TILING_SYSTOLIC)
#include <riscv_vector.h>
#include "../../common/vector_defines.h"

// matmul 最内部循环执行完毕一次会 load 矩阵 a 和 b 的 p 个数，即 16*p Byte。当 p <= 64 时，SPM <= 1KB.
void matrixmul_tiling_systolic(data_t *a, data_t *b, data_t *c, int n, int m, int p) {
    for (size_t i = 0; i < m; i++) {
        // 让 va 留在 systolic array 中，只当 c 换行时，即 i++ 时，此数据被读取，即每次读取的 i*p+k 都会存在 PE 的 regsiter 里
        // 那么，每一次读取，data reuse 的次数为 (p^2 - p)/2
        vfloat64m1_t va[p];
        for (int kk = 0; kk < p; kk += TILE_SIZE_K) {
            for (size_t k = kk; k < min(p, kk + TILE_SIZE_K); k += gvl){
                vfloat64m1_t va[k] = _MM_LOAD_f64(&a[i*p+k], gvl);
            }
        }
        for (size_t j = 0; j < n; j++) {
            size_t gvl = _MMR_VSETVL_E64M1(p);
            vfloat64m1_t vprod = _MM_SET_f64(0, gvl);
            vfloat64m1_t vsum  = _MM_SET_f64(0, gvl);
            data_t tmp_c[(TILE_SIZE_K / gvl)];
            for (int kk = 0; kk < p; kk += TILE_SIZE_K) {
                for (size_t k = kk; k < min(p, kk + TILE_SIZE_K); k += gvl){
                    gvl = _MMR_VSETVL_E64M1(p - k);

                    // Tile Load
                    // Matrix A row
                    // 只要第一个 C0j 读取了 va 的数据后，后面的所有 k 次循环都不需要再重新 load 了
                    vfloat64m1_t tmp_va = va[k];
                    // Matrix B column
                    vfloat64m1_t vb = _MM_LOAD_STRIDE_f64(&b[k*n+j], n * sizeof(data_t), gvl);

                    // Execution of A[0]*B[0], A[1]*B[1],... A[n]*B[n]
                    vprod  = _MM_MACC_f64(vprod,va, vb, gvl);

                }//k
                gvl = _MMR_VSETVL_E64M1(p);
                // A[0]*B[0] + A[1]*B[1],... + A[n]*B[n]
                vsum   = _MM_REDSUM_f64(vprod,vsum, gvl);
                tmp_c[k%kk] = vsum;
            }//kk
            // final C, we need tile accumulation
            for (int t = 0; t < (TILE_SIZE_K / gvl); t++) {
                c[i*n+j] += tmp_c[t]; // 按行访问
            }
        }//j
    }//i
}


#else // !USE_RISCV_VECTOR

void matmul_serial(data_t *a, data_t *b, data_t *c, int n, int m, int p) {
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < n; ++j) {
            c[i * n + j] = 0;
            for (int k = 0; k < p; ++k) {
                c[i * n + j] += a[i * p + k] * b[k * n + j];
            }
        }
}

#endif


bool compare( size_t dm, size_t dn, data_t *a ,data_t *b) {
    bool result = false;
    for (int i = 0; i < dm; i++) {
        for (int j = 0; j < dn; j++) {
            if(a[i*dn+j] != b[i*dn+j]) {
              result = true;
            }
        }

    }
    return result;
}
