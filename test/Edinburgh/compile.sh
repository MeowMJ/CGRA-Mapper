clang-12 -emit-llvm -O3 -fno-unroll-loops -fno-vectorize -o kernel.bc -c relu.c
clang-12 -emit-llvm -O3 -fno-unroll-loops -fno-vectorize -o conv.bc -c ../kernels/conv/conv.c
clang-12 -emit-llvm -O3 -fno-unroll-loops -fno-vectorize -o spmv.bc -c ../kernels/spmv/spmv.c
clang-12 -emit-llvm -O3 -fno-unroll-loops -fno-vectorize -o relu.bc -c ../kernels/relu/relu.c