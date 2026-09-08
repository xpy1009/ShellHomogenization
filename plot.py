import numpy as np
import matplotlib.pyplot as plt


n3 = np.loadtxt('../data/N3.txt')
n5 = np.loadtxt('../data/N5.txt')
n7 = np.loadtxt('../data/N7.txt')
h3 = np.loadtxt('../data/H3.txt')
h5 = np.loadtxt('../data/H5.txt')
h7 = np.loadtxt('../data/H7.txt')

plt.plot(n3[:, 0], n3[:, 1] * n3[:, 2], 'C0', linestyle=(0, (2, 2)), linewidth=3)
plt.plot(n5[:, 0], n5[:, 1] * n5[:, 2], 'C0', linestyle=(0, (4, 2)), linewidth=3)
plt.plot(n7[:, 0], n7[:, 1] * n7[:, 2], 'C0', linewidth=3)
plt.plot(h3[:, 0], h3[:, 1] * h3[:, 2], 'C1', linestyle=(0, (3, 3)), linewidth=2)
plt.plot(h5[:, 0], h5[:, 1] * h5[:, 2], 'C1', linestyle=(0, (6, 3)), linewidth=2)
plt.plot(h7[:, 0], h7[:, 1] * h7[:, 2], 'C1', linewidth=2)

plt.grid(True)

plt.show()