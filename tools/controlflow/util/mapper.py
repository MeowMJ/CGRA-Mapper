# ----------------------------------------------------------------------------
#   Filename: scheduler.py                                                  /
#   Description: simulate multi-kernel running on multi-CGRA                /
# ----------------------------------------------------------------------------

import heapq
import os
import subprocess
import json
import eventlet    # for time out
import pandas as pd
import math

# ----------------------------------------------------------------------------
#   global variables                                                        /
# ----------------------------------------------------------------------------

DICT_CSV = {'kernels': "", 'DFG nodes': "", 'DFG edges': "", 'recMII': "", 'mappingII': "", 'expandableII': "", 'utilization': ""}  # column names of generated CSV
DICT_COLUMN = len(DICT_CSV)
VECTOR_LANE = 2
JSON_NAME = "./param.json"
TIME_OUT_SET = 180
KERNEL_DIRECTORY = "../../test/EscortBench"
KERNEL_FUSION = False

def init_args(args):
    """init config"""
    global JSON_NAME, TIME_OUT_SET, KERNEL_DIRECTORY, KERNEL_FUSION
    JSON_NAME = args.json_name
    KERNEL_DIRECTORY = args.kernel_directory
    TIME_OUT_SET = args.time_out_set
    KERNEL_FUSION = args.fusion

def update_args(args):
    """init config"""
    global KERNEL_FUSION
    KERNEL_FUSION = args

# ----------------------------------------------------------------------------
#   class defination                                                         /
# ----------------------------------------------------------------------------

class Kernel:
    def __init__(self, kernel_name, kernel_id, arrive_period, unroll_factor, vector_factor, total_iterations, cgra_rows, cgra_columns, hardware_type):
        """
        Initialize an instance of the Kernel class.

        Parameters:
            kernel_name (str): The name of the kernel.
            kernel_id (int): The ID of the kernel.
            arrive_period (int): The period at which the same kernel will arrive again.
            unroll_factor (int): The unroll factor of the kernel.
            vector_factor (int): The vector factor of the kernel.
            total_iterations (int): The total number of iterations of the kernel.
            cgra_rows (int): The number of rows in the CGRA.
            cgra_columns (int): The number of columns in the CGRA.
        """
        self.kernel_name = kernel_name
        self.kernel_id = kernel_id
        self.arrive_period = arrive_period
        self.unroll_factor = unroll_factor
        self.vector_factor = vector_factor
        self.df = pd.DataFrame(DICT_CSV, index=[0])
        self.base_ii = 0  # II when using 1 CGRA, actual II, if fused, base_ii is fused_ii
        self.expandable_ii = 0  # II when using 2 CGRAs, expandable II, if fused, expandable_ii is individual_ii
        self.utilization = 0
        self.total_iterations = math.ceil(total_iterations / (self.unroll_factor*self.vector_factor))
        self.rows = cgra_rows
        self.columns = cgra_columns
        self.hardware_type = hardware_type
        self.load_kernel_data()


    def __lt__(self, other):
        """
        Compare two Kernel by id.
        """
        return self.kernel_id < other.kernel_id

    def load_kernel_data(self):
        prefix = './tmp/t_'
        csv_name = f'{prefix}{self.kernel_name}_{self.rows}x{self.columns}_unroll{self.unroll_factor}_vector{self.vector_factor}.csv'
        self.get_ii(csv_name)

        self.is_valid = bool(self.base_ii)
        print(f"Kernel {self.kernel_name} loaded with unroll_factor {self.unroll_factor}, vector factor {self.vector_factor}, cgra size {self.rows}x{self.columns}, hardware_type {self.hardware_type}")

    def comp_kernel(self):
        """
        This is a func compile a kernel using clang with selected unrolling factor.

        Returns: function name of kernel.
        """
        file_source = (self.kernel_name.split("."))[0]
        # corner case
        if self.kernel_name == "conv.c" and self.unroll_factor == 4:
            self.unroll_factor = 2
        if self.kernel_name == "fft.c" and self.unroll_factor == 2:
            self.unroll_factor = 1
        if self.kernel_name == "relu+histogram.c" and self.unroll_factor == 4 and self.rows == 12:
            self.unroll_factor = 2
        if self.kernel_name == "spmv.c" and self.unroll_factor == 2 and self.rows == 4:
            self.unroll_factor = 1

        if self.unroll_factor == 1 and self.vector_factor == 1:
            compile_command = f"clang-12 -emit-llvm -fno-unroll-loops -fno-vectorize -O3 -o kernel.bc -c {KERNEL_DIRECTORY}/{file_source}/{self.kernel_name}"
        elif self.unroll_factor == 1 and self.vector_factor != 1:
            compile_command = f"clang-12 -emit-llvm -fno-unroll-loops -O3 -mllvm -force-vector-width={self.vector_factor} -o kernel.bc -c {KERNEL_DIRECTORY}/{file_source}/{self.kernel_name}"
        elif self.unroll_factor != 1 and self.vector_factor == 1:
            compile_command = f"clang-12 -emit-llvm -funroll-loops -mllvm -unroll-count={self.unroll_factor} -fno-vectorize -O3 -o kernel.bc -c {KERNEL_DIRECTORY}/{file_source}/{self.kernel_name}"
        else:
            # print("Error, invalid unroll and vector factor combination.")
            return

        compile_proc = subprocess.Popen([compile_command, '-u'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        (compile_out, compile_err) = compile_proc.communicate()

        disassemble_command = f"llvm-dis-12 kernel.bc -o kernel.ll"
        disassemble_proc = subprocess.Popen([disassemble_command, '-u'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        (disassemble_out, disassemble_err) = disassemble_proc.communicate()


        if compile_err:
            print(f"Compile warning message for {self.kernel_name}: {compile_err}")
        if disassemble_err:
            # print(f"Disassemble error message for {self.kernel_name}: {disassemble_err}")
            return

        # collect the potentially targeting kernel/function from kernel.ll
        ir_file = open(f'kernel.ll', 'r')
        ir_lines = ir_file.readlines()

        # strips the newline character
        for line in ir_lines:
            if "define " in line and "{" in line and "@" in line:
                func_name = line.split("@")[1].split("(")[0]
                if "kernel" in func_name:
                    target_kernel = func_name
                    break

        ir_file.close()
        # print(f"Target kernel function for {self.kernel_name}: {target_kernel}")
        return target_kernel

    def map_kernel(self):
        """
        This is a func for mapping a kernel and gain information during mapping.

        Returns: NULL
        """
        get_map_command = f"opt-12 -load ../../build/src/libmapperPass.so -mapperPass kernel.bc"
        gen_map_proc = subprocess.Popen([get_map_command, "-u"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        dataS = []    # for get results from subprocess and output to pandas
        kernels_source = (self.kernel_name.split("."))[0]
        dataS.append(kernels_source)

        try:
            eventlet.monkey_patch()
            with eventlet.Timeout(TIME_OUT_SET, True):
                with gen_map_proc.stdout:
                    gen_map_proc.stdout.flush()
                    for line in iter(gen_map_proc.stdout.readline, b''):
                        output_line = line.decode("ISO-8859-1")
                        if "DFG node count: " in output_line:
                            dataS.append(int(output_line.split("DFG node count: ")[1].split(";")[0]))
                            dataS.append(int(output_line.split("DFG edge count: ")[1].split(";")[0]))
                        if "[RecMII: " in output_line:
                            dataS.append(int(output_line.split("[RecMII: ")[1].split("]")[0]))
                        if "[Mapping II: " in output_line:
                            self.base_ii = int(output_line.split("[Mapping II: ")[1].split("]")[0])
                            dataS.append(self.base_ii)
                        if "[ExpandableII: " in output_line:
                            self.expandable_ii = int(output_line.split("[ExpandableII: ")[1].split("]")[0])
                            dataS.append(self.expandable_ii)
                        if "tile avg fu utilization: " in output_line:
                            self.utilization = min(float(output_line.split("avg overall utilization: ")[1].split("%")[0])/100,1)
                            dataS.append(self.utilization)
                        if "[Mapping Fail]" in output_line:
                            print(f"{self.kernel_name} mapping failed.")
        except eventlet.timeout.Timeout:
            dataS = [0]*(DICT_COLUMN)
            # print("Skipping a specific config for kernel: ", self.kernel_name, "Because it runs more than", TIME_OUT_SET/60 , "minute(s).")

        if len(dataS) != DICT_COLUMN:
            dataS.extend([0]*(DICT_COLUMN-len(dataS)))

        self.df.loc[len(self.df.index)] = dataS

    def map_kernel_skip(self):
        """
        This is a func gain DFG information only without mapping.

        Returns: NULL
        """
        get_map_command = f"opt-12 -load ../../build/src/libmapperPass.so -mapperPass kernel.bc"
        gen_map_proc = subprocess.Popen([get_map_command, "-u"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        # Holds the results from subprocess and output to pandas.
        dataS = []
        kernels_source = (self.kernel_name.split("."))[0]
        dataS.append(kernels_source)
        # The first 4 element of dataS is not empty: kernelsSource, DFG node count, DFG edge count, RecMII.
        k_data_s_head = 4

        try:
            eventlet.monkey_patch()
            with eventlet.Timeout(TIME_OUT_SET, True):
                with gen_map_proc.stdout:
                    gen_map_proc.stdout.flush()
                    for line in iter(gen_map_proc.stdout.readline, b''):
                        output_line = line.decode("ISO-8859-1")
                        if "DFG node count: " in output_line:
                            dataS.append(int(output_line.split("DFG node count: ")[1].split(";")[0]))
                            dataS.append(int(output_line.split("DFG edge count: ")[1].split(";")[0]))
                        if "[RecMII: " in output_line:
                            dataS.append(int(output_line.split("[RecMII: ")[1].split("]")[0]))
                            dataS.extend([0]*(DICT_COLUMN-k_data_s_head))
                            break

        except eventlet.timeout.Timeout:
            dataS = [0]*(DICT_COLUMN)
            # print("Skipping a specific config for kernel: ", self.kernel_name, "Because it runs more than", TIME_OUT_SET/60, "minute(s).")

        self.df.loc[len(self.df.index)] = dataS

    def get_ii(self, csv_name):
        """
        This is a func to compile, run and map kernels under neura_json and store the mapping result in csv

        Returns: name of the csv that collects information of mapped kernels
        """
        # print("Generating", csv_name)
        target_kernel = self.comp_kernel()
        target_fusion_strategy = []
        if self.hardware_type == 0:
            target_fusion_strategy = ["default_heterogeneous"]
        else:
            target_fusion_strategy = []
        neura_json = {
            "kernel": target_kernel,
            "targetFunction": False,
            "targetNested": False,
            "targetLoopsID": [0],
            "doCGRAMapping": True,
            "row": self.rows,
            "column": self.columns,
            "precisionAware": False,
            "fusionStrategy": target_fusion_strategy,
            "isTrimmedDemo": True,
            "heuristicMapping": True,
            "parameterizableCGRA": False,
            "vectorizationMode": "all",
            "diagonalVectorization": False,
            "bypassConstraint": 4,
            "isStaticElasticCGRA": False,
            "ctrlMemConstraint": 10,
            "regConstraint": 8,
            "incrementalMapping"    : False,
            "vectorFactorForIdiv "  : 1,
            "fusionPattern"         : {
                                        "4" : ["phi", "add", "icmp", "br"]
                                        },
            "additionalFunc"        : {
                                        "load" : [0,1,2,3],
                                        "store": [0,1,2,3]
                                        },
            "incrementalMapping"   : False,
            "pathSupportDim"       : 4, # TODO
            "ctrlType"             : self.hardware_type,
            "supportDVFS": False,
            "DVFSIslandDim": 1,
            "DVFSAwareMapping": False,
            "enablePowerGating": False,
            "expandableMapping" : True
        }



        json_object = json.dumps(neura_json, indent=4)

        with open(JSON_NAME, "w") as outfile:
            outfile.write(json_object)
        if True:
            self.map_kernel()
        else:
            self.map_kernel_skip()

        self.df.to_csv(csv_name)
        return csv_name

    def read_ii(self, csv_name):
        """
        This is a func to read from csv generated from get_ii()

        Returns: csv_name
        """
        try:
            df = pd.read_csv(csv_name)
            self.base_ii = int(df['mappingII'].iloc[1])
            self.expandable_ii = int(df['expandableII'].iloc[1])
            if 'utilization' in df.columns:
                self.utilization = min(float(df['utilization'].iloc[1]),1.0)
            else:
                self.get_ii()
                return csv_name
        except FileNotFoundError:
            # print(f"CSV file {csv_name} not found.")
            self.get_ii()
            return csv_name
        except ValueError:
            # print(f"Error extracting II values from {csv_name}.")
            self.get_ii()
            return csv_name

        return csv_name

    def return_ii(self, num_cgras):
        """
        Get the initiation interval (II) based on the number of CGRAs allocated.

        Parameters:
            num_cgras (int): Number of CGRAs allocated.

        Returns:
            int: The initiation interval (II).
        """
        if num_cgras == 1:
            return self.base_ii
        elif num_cgras == 2:
            return self.expandable_ii
        else:
            raise ValueError("Number of CGRAs must be 1 or 2.")

    def return_total_iterations(self):
        """
        Total iterations for the kernel, affected by unroll_factor and vector_factor

        Returns:
            int: Total iterations.
        """
        return self.total_iterations