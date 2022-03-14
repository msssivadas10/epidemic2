import numpy as np
import matplotlib.pyplot as plt
plt.style.use('ggplot')

n = 500
dt = 0.1

pos = np.random.uniform(0, 100, (n, 2))

theta = np.random.uniform(0, 2*np.pi, n)
vel =  np.array([np.cos(theta), np.sin(theta)]).T

# q = -np.ones(n) * 10


fig, ax = plt.subplots(figsize=[6,6])

while 1:
    ax.cla()

    ax.plot(pos[:,0], pos[:,1], 'o', ms = 3, color = 'black')

    ax.set(xlim=[0,100], ylim=[0,100])

    pos = pos + vel * dt * 0.5

    r = pos - 50.0
    _r = np.sqrt(r[:,:1]**2 + r[:,1:]**2 + 0.1)

    # acc = q[:,None] * 10 * r / (r[:,:1]**2 + r[:,1:]**2 + 0.1)**1.5 

    acc = -10 * np.exp(-(_r/100)**2 / 0.05) * r / _r

    vel = vel + acc * dt

    pos = pos + vel * dt * 0.5

    mask = (pos > 100) | (pos < 0)
    vel[mask] = -vel[mask]

    # vel = vel / (vel[:,:1]**2 + vel[:,1:]**2)**0.5

    plt.pause(0.05)

# x, y = np.mgrid[0:100:51j, 0:100:51j]

# r = ((x - 50)**2 + (y - 50)**2 + 0.1)**0.5


# f = 0.1 * np.exp(-(r/100)**2 / 0.05) 

# fx = (x-50)*f / r
# fy = (y-50)*f / r

# a = ax.pcolor(x, y, np.sqrt(fx**2+fy**2))

# plt.colorbar(a)


plt.show()