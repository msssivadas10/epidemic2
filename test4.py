from turtle import color
import numpy as np
import numpy.random as rnd
import matplotlib.pyplot as plt
from matplotlib.widgets import Button

npart = 50
scale = 1.0
stop  = False
dt    = 0.005
r0    = scale / npart**0.5

def button_action(x):
    global stop
    stop = True
    plt.close()
    return

x = rnd.uniform(size = (npart, 2)) * scale
v = rnd.uniform(size = (npart, 2)) * 2 - 1

fig, ax = plt.subplots(figsize = (6,6))
axb     = plt.axes([0.45, 0.9, 0.1, 0.075])
button  = Button(axb, 'Stop')
button.on_clicked(button_action)

while not stop:
    ax.cla()
    ax.scatter(x[:,0], x[:,1], s = 7, color = 'black')
    ax.set(xlim = [0, scale], ylim = [0, scale])
    plt.pause(0.01)

    t0 = np.mean(v**2)

    # compute force:
    dx = x[..., None] - x.T
    r = np.sqrt(np.sum(dx**2, 1))
    a  = np.zeros_like(r)
    m  = (r != 0)
    a[m] = (r0 / r[m])**3
    a[m]  = (24 * a[m] * (2 * a[m] - 1) / r[m])
    a = np.sum(dx @ a, -1)


    # update positions:
    x += v * dt / 2.0
    v += a * dt
    x += v * dt / 2.0
    print(a)


    # v *= np.sqrt(t0 / np.mean(v**2))

    mask = (x < 0.0)
    x[mask] = -x[mask]; v[mask] = -v[mask]

    mask = (x > scale)
    x[mask] = 2 * scale - x[mask]; v[mask] = -v[mask]

plt.show()
