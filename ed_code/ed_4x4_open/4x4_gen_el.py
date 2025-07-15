import numpy as np
from scipy.linalg import block_diag
import itertools
import scipy as sc
import scipy.sparse as scs
from joblib import Parallel, delayed
import sys
import copy

def ten2four(n):
    n = n % (4294967296)
    s = np.base_repr(n,4)
    return np.array(list(s.zfill(16)),dtype=int)

def e_count(state):
    el = np.count_nonzero(state)
    return 3*el #6/4*2

def get_H_el(gauge_indices_base4):
    return gauge_indices_base4.astype(bool).sum(1)*3

def main():
    gauge_indices = np.load("gauge_indices_4x4_sort.npz")["arr_0"].astype(int)
    gauge_indices_base4 = np.array([ten2four(i) for i in gauge_indices])
    h_el = get_H_el(gauge_indices_base4)
    h2 = scs.block_diag(h_el)
    scs.save_npz("h_el",h2)



main()
quit()