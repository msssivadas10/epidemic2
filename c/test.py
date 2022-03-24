import numpy as np
import matplotlib.pyplot as plt


t, s, i, r = np.loadtxt('output.txt', unpack = 1)


plt.figure()
plt.plot(t, s)
plt.plot(t, i)
plt.plot(t, r)
plt.show()