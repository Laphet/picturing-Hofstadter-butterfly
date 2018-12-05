from scipy.sparse import diags
import numpy as np
from numpy import linalg as LA
import math
import timeit

order = 8
ondiagonal = np.array([2.0*math.cos(2.0 * math.pi * j * 0.34)
                       for j in range(order)])
diagonals = [ondiagonal, np.ones(order-1), np.ones(order-1)]
tri_mat = diags(diagonals, offsets=[0, 1, -1]).toarray()
start = timeit.default_timer()
w = LA.eigvalsh(tri_mat)
end = timeit.default_timer()
print(w)
print('using default QR algorithm in numpy to solve the eigenvalues of a n={:d} tridiagonal matrix consumes time {:f}s\n'.format(
    order, end-start))
