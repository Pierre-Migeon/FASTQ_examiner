
import numpy as np



N_count = np.zeros((500,), dtype=int)
T_count = np.zeros((500,), dtype=int)



N_count[0] = 25
T_count[0] = 4

N_count[0] = N_count[0] / T_count[0]

print(N_count[0])
