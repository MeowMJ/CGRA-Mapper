clang-12 -O0 -emit-llvm -fno-unroll-loops -fno-vectorize kernel.cpp -S -o kernelO0.ll
clang-12 -O0 -emit-llvm -fno-unroll-loops -fno-vectorize ../kernels/conv/conv.c -S -o convO0.ll
clang-12 -O0 -emit-llvm -fno-unroll-loops -fno-vectorize ../kernels/spmv/spmv.c -S -o spmvO0.ll
clang-12 -O0 -emit-llvm -fno-unroll-loops -fno-vectorize ../kernels/relu/relu.c -S -o reluO0.ll
clang-12 -O3 -emit-llvm -fno-unroll-loops -fno-vectorize kernel.cpp -S -o kernelO3.ll
clang-12 -O3 -emit-llvm -fno-unroll-loops -fno-vectorize ../kernels/conv/conv.c -S -o convO3.ll
clang-12 -O3 -emit-llvm -fno-unroll-loops -fno-vectorize ../kernels/spmv/spmv.c -S -o spmvO3.ll
clang-12 -O3 -emit-llvm -fno-unroll-loops -fno-vectorize ../kernels/relu/relu.c -S -o reluO3.ll