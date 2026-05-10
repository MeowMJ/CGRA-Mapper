<<<<<<< HEAD
# clang-12 -emit-llvm -fno-unroll-loops -fno-discard-value-names -O3 -o kernel.bc -c 1113IFELMY.cpp
clang-12 -O3 -emit-llvm -fno-unroll-loops -fno-discard-value-names 1113IFELMY.cpp -S -o kernel.ll
#llvm-dis fir.bc -o fir.ll
=======
clang-12 -emit-llvm -fno-unroll-loops -O0 -o kernel.bc -c kernel.cpp
llvm-dis-12 kernel.bc -o O0kernel.ll
#clang-12 -emit-llvm -fno-unroll-loops -mllvm -force-vector-width=4 -O3 -o kernel.bc -c ./_matmul/src/matmul.c
>>>>>>> f5eeba755245451367251c9a787615a8c681ef2d
