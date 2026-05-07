#makefile

# Variables for the compiler and options
RISCV_BIN=${RISCV}/bin
PK=${RISCV}/riscv64-unknown-elf/bin/pk

# Variables for source code
SRC=src/spmv.c src/main.c
BIN_DIR=bin
BIN_SERIAL=$(BIN_DIR)/spmv_serial.exe
BIN_VECTOR=$(BIN_DIR)/spmv_vector.exe
DUMP_SERIAL=$(BIN_DIR)/spmv_serial.dump
DUMP_VECTOR=$(BIN_DIR)/spmv_vector.dump

start:
	rm -f bin/*.*

serial:
	@echo "Compiling serial version..."
	$(RISCV_BIN)/riscv64-unknown-elf-gcc -march=rv64gc -O2 -o $(BIN_SERIAL) $(SRC)
	# $(RISCV_BIN)/llvm-objdump --mattr=+m,+f,+d,+a,+c -ds $(BIN_SERIAL) > $(DUMP_SERIAL)

vector:
	@echo "Compiling vector version..."
	$(RISCV_BIN)/riscv64-unknown-elf-gcc -march=rv64gcv -O2 -I../common -o $(BIN_VECTOR) $(SRC)
	$(RISCV_BIN)/riscv64-unknown-elf-objdump -ds $(BIN_VECTOR) > $(DUMP_VECTOR)	

clean:
	rm -f bin/*.dump; \
	rm -f bin/*.exe;
