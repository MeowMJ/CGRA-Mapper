# clang-12 -emit-llvm -fno-unroll-loops -fno-discard-value-names -O3 -o kernel.bc -c 1113IFELMY.cpp
clang-12 -O3 -emit-llvm -fno-unroll-loops -fno-discard-value-names 1113IFELMY.cpp -S -o kernel.ll
#llvm-dis fir.bc -o fir.ll
