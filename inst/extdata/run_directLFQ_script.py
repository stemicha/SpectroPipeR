
# import directlfq.normalization as lfqnorm
# import directlfq.protein_intensity_estimation as lfqprot_estimation
# import directlfq.utils as lfqutils
# import pandas as pd
# import numpy as np
import directlfq.lfq_manager as lfqmanager

import warnings
import sys
import multiprocess


def run_pipline(input_file, num_cores):
  df1=lfqmanager.run_lfq(input_file, num_cores = num_cores)
  
if __name__ == '__main__':
  
  
  if len(sys.argv) > 1:
    input_file = sys.argv[1]
  else: 
    print("Lack input!")
    sys.exit(1) ## if miss input, interrupt and return 1
    
  if len(sys.argv) > 2:
    print(len(sys.argv))
    n_cores = int(sys.argv[2])
  else: 
    n_cores = 32
    
  if n_cores > multiprocess.cpu_count():
    n_cores = multiprocess.cpu_count() - 1
  
  print("Input file:"," ",input_file)
  print("Use", n_cores, "cores")
  
  run_pipline(input_file = input_file, num_cores = n_cores)
  sys.exit(0)
