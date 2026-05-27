# ----------------------------------------------------------------------------
#   Filename: main.py                                                       /
#   Description: heterogenous control flow support on CGRAs                 /
# ----------------------------------------------------------------------------

import argparse
import json
import os
from pathlib import Path
import time
import util.mapper as mapper
import util.visualizer as visualizer

# ----------------------------------------------------------------------------
#   global variables                                                        /
# ----------------------------------------------------------------------------
VISUALIZATION = True
TESTME = False
FUSION = False
RESULT = ''

# Static kernel data (name: (sort_id, total_iterations))
KERNEL_DATA = {
    "adpcm_coder_kernel.cpp": (0, 1024),
    "basicmath.cpp": (0, 2048),
    "CBT_kernel.cpp": (0, 512),
    "consumer.cpp": (0, 1024),
    "CreateBinaryTree_kernel.cpp": (0, 512),
    "nwkernel_kernel.cpp": (0, 1024),
    "office.cpp": (0, 2048),
    "patricia.cpp": (0, 1024),
    "relu.c": (0, 4096),
    "security.cpp": (0, 1024),
    "solver0.cpp": (0, 2048),
    "susan.c": (0, 1024),
    "wtree_per_kernel.cpp": (0, 1024),
}

# Case configuration dictionary (task_id: kernels, cgra_rows, cgra_cols — all same length)
TASK_CONFIGS = {
    # Task 1: small kernels on 4x4 CGRA
    1: {
        'KERNELS': ['adpcm_coder_kernel.cpp', 'basicmath.cpp', 'CBT_kernel.cpp', 'consumer.cpp', 'CreateBinaryTree_kernel.cpp', 'nwkernel_kernel.cpp', 'office.cpp', 'patricia.cpp', 'relu.c', 'security.cpp', 'solver0.cpp', 'susan.c', 'wtree_per_kernel.cpp'],
        'CGRA_ROWS': [4, 4, 4, 4],
        'CGRA_COLS': [4, 4, 4, 4]
    },
    # Task 2: mix of kernels on varied CGRA sizes
    2: {
        'KERNELS': ['fir.cpp', 'latnrm.c', 'dtw.cpp', 'mvt.c', 'relu+histogram.c'],
        'CGRA_ROWS': [4, 4, 6, 6, 8],
        'CGRA_COLS': [4, 4, 6, 6, 8]
    },
    # Task 3: all kernels on 8x8 CGRA
    3: {
        'KERNELS': ['fir.cpp', 'latnrm.c', 'fft.c', 'dtw.cpp', 'spmv.c', 'conv.c', 'mvt.c', 'gemm.c', 'relu+histogram.c'],
        'CGRA_ROWS': [8]*9,
        'CGRA_COLS': [8]*9
    },
    # Task 4: heavy kernels on larger 12x12 CGRA
    4: {
        'KERNELS': ['conv.c', 'gemm.c', 'relu+histogram.c', 'fft.c', 'spmv.c'],
        'CGRA_ROWS': [12, 12, 12, 8, 8],
        'CGRA_COLS': [12, 12, 12, 8, 8]
    }
}


# ----------------------------------------------------------------------------
#   function defination                                                      /
# ----------------------------------------------------------------------------

def str_to_bool(value):
    if isinstance(value, bool):
        return value
    if str(value).lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif str(value).lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    raise argparse.ArgumentTypeError('Invalid boolean value (accepted: 0/1, true/false, yes/no)')


def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Multi-CGRA Task Scheduling Tool'
    )
    # Core application arguments
    parser.add_argument('--test', type=str_to_bool, default=TESTME,
                       help='Run tests in CI/CD [y/n]')
    parser.add_argument('--cgra-config', type=int, default= 4,
                       help='Path to CGRA configuration file')
    parser.add_argument('--json-name', type=str, default= "./param.json",
                       help='JSON configuration file name')
    parser.add_argument('--kernel-directory', type=str, default= "../../test/EscortBench",
                       help='Kernel directory path')
    parser.add_argument('--result-directory', type=str, default= "../../test/EscortBench",
                    help='RESULT directory path')
    parser.add_argument('--time-out-set', type=int, default= 180,
                       help='Timeout setting for operations')
    parser.add_argument('--visualize', type=str_to_bool, default=VISUALIZATION,
                       help='Generate visualization figures [y/n]')
    parser.add_argument('--fusion', type=str_to_bool, default=FUSION,
                       help='Default fusion strategy for kernels [y/n]')

    return parser.parse_args()


def load_configuration():
    """Load and merge configurations from multiple sources with priority:
    1. Command line arguments (highest priority)
    2. Default values (lowest priority)
    """
    # Update global configuration with command line arguments
    global VISUALIZATION, TESTME, FUSION
    # Parse command line arguments
    args = parse_arguments()
    VISUALIZATION = args.visualize
    TESTME = args.test
    FUSION = args.fusion
    mapper.init_args(args)
    print(f"Test in CI/CD: {args.test}")
    print(f"Timeout: {args.time_out_set}")
    print(f"Visualization: {args.visualize}")
    print(f"FUSION: {args.fusion}")


# ========== Task Loading Function ==========

# def load_tasks(task_id, task_type="baseline"):
#     """
#     Load task list based on task_id and CGRA type

#     Args:
#         task_id: Configuration case ID
#         task_type: "baseline" or "task", corresponding to 12x12 and 4x4 CGRA respectively

#     Returns:
#         task_list: List of task objects
#     """
#     global TASK_CONFIGS, KERNEL_DATA
#     if task_id not in TASK_CONFIGS:
#         raise ValueError(f"Task{task_id} configuration does not exist")

#     config = TASK_CONFIGS[task_id]
#     A_P = config['A_P']
#     UNROLL_FACTORS = config['UNROLL_FACTORS']
#     VECTOR_FACTORS = config['VECTOR_FACTORS']

#     # Validate parameter lengths
#     lists = [KERNEL_DATA, A_P, UNROLL_FACTORS, VECTOR_FACTORS]
#     if len(set(len(lst) for lst in lists if lst)) > 1:
#         raise ValueError(f"Task{task_id} parameter length mismatch: {[len(lst) for lst in lists]}")

#     # Set CGRA dimensions
#     if task_type == "baseline":
#         cgra_rows, cgra_columns = 12, 12
#     elif task_type == "task":
#         cgra_rows, cgra_columns = 4, 4
#     else:
#         raise ValueError("task_type must be either 'baseline' or 'task'")

#     # Generate task list
#     task_list = []
#     for i, (kernel_name, (kernel_id, total_iters, _)) in enumerate(KERNEL_DATA.items()):
#         task = mapper.Kernel(
#             kernel_name=kernel_name,
#             kernel_id=kernel_id,
#             arrive_period=A_P[i] if A_P else 0,
#             unroll_factor=UNROLL_FACTORS[i],
#             vector_factor=VECTOR_FACTORS[i],
#             total_iterations=total_iters,
#             cgra_rows=cgra_rows,
#             cgra_columns=cgra_columns
#         )
#         task_list.append(task)

#     return task_list


# def load_tasks_from_file(filename):
#     """
#     Load task list from JSON file

#     Args:
#         filename: Input JSON filename

#     Returns:
#         task_list: List of reconstructed task objects
#     """
#     if not os.path.exists(filename):
#         raise FileNotFoundError(f"Task file {filename} not found")

#     with open(filename, 'r') as f:
#         tasks_data = json.load(f)

#     # Reconstruct task objects from dictionaries
#     task_list = []
#     for task_dict in tasks_data:
#         task = mapper.Kernel(
#             kernel_name=task_dict['kernel_name'],
#             kernel_id=task_dict['kernel_id'],
#             arrive_period=task_dict['arrive_period'],
#             unroll_factor=task_dict['unroll_factor'],
#             vector_factor=task_dict['vector_factor'],
#             total_iterations=task_dict['total_iterations'],
#             cgra_rows=task_dict['cgra_rows'],
#             cgra_columns=task_dict['cgra_columns']
#         )
#         task_list.append(task)

#     print(f"Tasks loaded from {filename}")
#     return task_list


# ========== Run Function ==========

def run_kernels(task_type, hardware_type):
    """
    Run kernels on 12x12 CGRA baseline and collect results.

    Kernel list and CGRA dimensions are read directly from TASK_CONFIGS.
    """
    config = TASK_CONFIGS[task_type]
    kernel_names = config['KERNELS']
    cgra_rows = config['CGRA_ROWS']
    cgra_cols = config['CGRA_COLS']

    for i, name in enumerate(kernel_names):
        _, total_iters = KERNEL_DATA[name]
        mapper.Kernel(
            kernel_name=name,
            kernel_id=i,
            arrive_period=0,
            unroll_factor=1,
            vector_factor=1,
            total_iterations=total_iters,
            cgra_rows=cgra_rows[i],
            cgra_columns=cgra_cols[i],
            hardware_type=hardware_type
        )


def main():
    """Main workflow control function"""
    start = time.time()
    # 1. Load configuration (includes parsing arguments)
    print("=== Control Flow Task Scheduling Tool ===")
    load_configuration()

    # 2. Create output directory
    print(f"Intermediate result in: ./tmp")
    print(f"Final result in: {RESULT}")
    output_dir = Path("./tmp")
    output_dir.mkdir(parents=True, exist_ok=True)

    # 3. Execute scheduling — iterate over each (task_type, hardware_type) pair
    HARDWARE_TYPE = ["baseline", "escort", "4dcgra"] if not TESTME else ["baseline"]
    print(f"[Step 1] Loading tasks and Running tasks on CGRAs...")
    if TESTME:
        run_kernels(task_type=1, hardware_type=2)  # baseline 对应 2
    else:
        hardware_type_map = {"escort": 0, "4dcgra": 1, "baseline": 2}
        for task_type, hw_str in zip(TASK_CONFIGS, HARDWARE_TYPE):
            # task_type 为 TASK_CONFIGS 的 key（1, 2, 3, 4）
            # hw_str 为 ["baseline", "escort", "4dcgra"] 中的字符串
            # 映射为数字：escort→0, 4dcgra→1, baseline→2
            run_kernels(task_type=task_type, hardware_type=hardware_type_map[hw_str])

        # 5. Generate visualization
        if VISUALIZATION:  # Use global variable
            print(f"[Step 3] Generating visualization figures...")
            # Generate Fig9
            genFigs = visualizer.SimulationDataAnalyzer(kernel_data=KERNEL_DATA, output_dir=RESULT)
            genFigs.genFig9("./fig/Fig9.png")
            genFigs.genFig10("./fig/Fig10.png")
            genFigs.genFig11("./fig/Fig11.png")


    print("\n=== Scheduling completed successfully! ===")
    end = time.time()
    execution_time = end - start
    print(f"Time cost: {execution_time/60:.2f} min")


if __name__ == '__main__':
    main()