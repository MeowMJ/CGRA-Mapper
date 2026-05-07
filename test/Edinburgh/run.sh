opt-12 -load ../../build/src/libmapperPass.so -mapperPass kernel.bc | tee trace.log
dot -Tpng _Z6kernelPfS_S_.dot -o kernel.png
opt-12 -load ../../build/src/libmapperPass.so -mapperPass conv.bc
dot -Tpng kernel.dot -o conv.png
rm kernel.dot
opt-12 -load ../../build/src/libmapperPass.so -mapperPass spmv.bc
dot -Tpng kernel.dot -o spmv.png
opt-12 -load ../../build/src/libmapperPass.so -mapperPass relu.bc
dot -Tpng kernel.dot -o relu.png