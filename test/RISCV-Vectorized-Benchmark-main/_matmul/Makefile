#makefile

# Variables for the compiler and options
RISCV_BIN=${RISCV}/bin

start:
	rm -f bin/*.*

serial:
	@echo "serial"
	${RISCV_BIN}/bin/clang -march=rv64gcv -O2 -fno-vectorize -o bin/matmul_serial.exe src/matmul.c src/main.c
	${RISCV_BIN}/bin/llvm-objdump  --mattr=+m,+f,+d,+a,+c -ds  bin/matmul_serial.exe > bin/matmul_serial.dump

vector:
	@echo "vector"
	$(RISCV_BIN)/riscv64-unknown-elf-gcc -march=rv64gcv -O2 -I../common -o bin/matmul_vector.exe src/matmul.c src/main.c
	$(RISCV_BIN)/riscv64-unknown-elf-objdump -ds bin/matmul_vector.exe > bin/matmul_vector.dump

clean:
	rm -f bin/*.dump; \
	rm -f bin/*.exe;
