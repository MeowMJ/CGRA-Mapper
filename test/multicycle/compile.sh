clang-12 -emit-llvm -O3 -fno-unroll-loops -fno-vectorize -o multicycle_test.bc -c multicycle_test.cpp
opt-12 --loop-unroll --unroll-count=4 multicycle_test.bc -o multicycle_test.bc