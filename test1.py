import numpy as np
import numpy.random as rnd
import matplotlib.pyplot as plt
plt.style.use('ggplot')
import matplotlib.animation as anim
from pprint import pprint


def make_grid(pos, subdiv, scale):
    cs    = scale / subdiv
    cells = (pos // cs).astype('int')
    grid  = {}
    for k, (i,j) in enumerate(cells):
        if (i,j) not in grid.keys():
            grid[i,j] = []
        grid[i,j].append(k)
    
    # remove empty
    _grid = {}
    for key, val in grid.items():
        if len(val):
            _grid[key] = np.array(val)        
    return _grid

def get_nearest(grid, pos, x, r, scale, subdiv):
    cs = scale / subdiv
    cr = -int(-r // cs)
    ci, cj = (np.asarray(x) // cs).astype('int')
    near = []
    for i in range(ci-cr, ci+cr+1):
        for j in range(cj-cr, cj+cr+1):
            if (i,j) not in grid.keys():
                continue
            k = grid[i,j]
            nn = ( np.sum((pos[k,:] - x)**2, -1) < r**2 )
            near = near + list(k[nn])
    return np.array(near)

def a_attr(pos, scale, f0, fpos, w = 0.05):
    r  = pos - fpos
    _r = np.sqrt(r[:,:1]**2 + r[:,1:]**2 + 0.1)

    acc = -f0 * np.exp(-(_r/scale)**2 / w) * r / _r
    return acc

def init(n, scale, subdiv):
    pos   = rnd.uniform(0, 1, (n, 2)) * scale

    theta = rnd.uniform(0, 2*np.pi, n)
    vel   =  np.array([np.cos(theta), np.sin(theta)]).T
    
    grid  = make_grid(pos, subdiv, scale)
    return pos, vel, grid

def update(pos, vel, dt, f0 = 0.0, fpos = None, fvel = None, ):
    pos[...] += vel[...] * dt 

    if fpos is not None:
        acc = a_attr(pos, scale, f0, fpos, 0.05)
        vel[...] += acc * dt

    # pos = pos + vel * dt * 0.5

    mask = (pos < 0)
    pos[mask] = -pos[mask]
    vel[mask] = -vel[mask]

    mask = (pos > scale)
    pos[mask] = 2 * scale - pos[mask]
    vel[mask] = -vel[mask]


    if fvel is not None:
        fpos[...] += fvel[...] * dt

        mask = (fpos < 0)
        fpos[mask] = -fpos[mask]
        fvel[mask] = -fvel[mask]

        mask = (fpos > scale)
        fpos[mask] = 2 * scale - fpos[mask]
        fvel[mask] = -fvel[mask]

    # vel = vel / (vel[:,:1]**2 + vel[:,1:]**2)**0.5

    # vel = np.clip(vel, -1, 1)
    grid = make_grid(pos, subdiv, scale)
    return grid

def init_disease(n, pi, ri, T):
    state = np.zeros(n, 'int')
    t2rec = np.zeros(n)

    # infect random:
    i = rnd.randint(n)
    state[i] = 1
    t2rec[i] = T
    return state, t2rec

def update_diease(pos, grid, state, t2rec, pi, ri, T, scale, subdiv):
    inf = (state == 1)
    t2rec[inf] -= dt

    # transmission:
    for xi in pos[inf,:]:
        nn = get_nearest(grid, pos, xi, ri, scale, subdiv)

        if not len(nn):
            continue
        to_inf        = nn[(rnd.random(size = len(nn)) < pi) & (state[nn] == 0)]
        state[to_inf] = 1
        t2rec[to_inf] = T

    rec = (t2rec < 0.0)
    t2rec[rec] = 0.0
    state[rec] = 2
    return 

def count_stats(state):
    return [np.sum(state==0), np.sum(state==1), np.sum(state==2)]
    
def simulation(scale, n, dt, subdiv, f0, fc, fv, pi, ri, T, tmax = None):
    pos, vel, grid = init(n, scale, subdiv)
    state, t2rec   = init_disease(n, pi, ri, T)

    t, sir  = [], []

    ti = 0
    while (1 if tmax is None else ti < tmax):
        t.append(ti)
        stat = count_stats(state)
        sir.append(stat)

        grid = update(pos, vel, dt, f0, fc, fv)

        update_diease(pos, grid, state, t2rec, pi, ri, T, scale, subdiv)
        ti += 1

        if not stat[1]:
            break
    return t, sir

def mean_peakinf(rep, scale, n, dt, subdiv, f0, fc, fv, pi, T, tmax = None):
    tm, im = 0.0, 0.0
    for j in range(rep):
        t, sir = simulation(scale, n, dt, subdiv, f0, fc, fv, pi, ri, T, tmax)
        
        sir = np.array(sir)
        k   = np.argmax(sir[:,1])
        tm += t[k]
        im += sir[k, 1]
    return tm / rep, im / rep 


scale  = 10
n      = 500
dt     = 0.1
subdiv = 10

f0, fc, fv = 100, np.array([0.5, 0.5]) * scale, np.array([0.0, 0.0])

pi, ri, T  = 0.1, 0.05 * scale, 50 * dt

# pos, vel, grid = init(n, scale, subdiv)

# state, t2rec   = init_disease(n, pi, ri, T)

# t, sir = [], []

# t, sir = simulation(scale, n, dt, subdiv, f0, fc, fv, pi, ri, T)
# sir = np.array(sir)

fig, ax = plt.subplots(figsize=[6,6])


rs = np.linspace(0.0, 0.2, 101)
for f0 in np.linspace(0, 100, 5):
    tm, im = [], []
    for ri in rs:
        ri = ri * scale
        _tm, _im = mean_peakinf(500, scale, n, dt, subdiv, f0, fc, fv, pi, ri, T)
        tm.append(_tm)
        im.append(_im)
    ax.plot(rs, im, 'o-', ms = 2, label = f"{f0:.3f}")
ax.legend()
# colours = np.array(['black', 'C0', 'gray'])


# ax.plot(t, sir[:,0], color = colours[0])
# ax.plot(t, sir[:,1], color = colours[1])
# ax.plot(t, sir[:,2], color = colours[2])

# ax.plot(f0s, im, 'o-')
# ax.plot(f0s, tm, 'o-')


# ax.set(xlim = [0, scale], ylim = [0, scale])

# l0, = ax.plot([], [], 'o', ms = 1.5, color = 'black')
# l1, = ax.plot([], [], 'o', ms = 1.5, color = 'C0')
# l2, = ax.plot([], [], 'o', ms = 1.5, color = 'gray')

# def update_plot(i):
#     l0.set_data(pos[state==0, 0], pos[state==0, 1])
#     l1.set_data(pos[state==1, 0], pos[state==1, 1])
#     l2.set_data(pos[state==2, 0], pos[state==2, 1])

#     grid = update(pos, vel, f0, fc, fv)

#     update_diease(pos, grid, state, t2rec, pi, ri, T, scale, subdiv)

# ani = anim.FuncAnimation(fig, update_plot, frames=1000, interval=200, blit=False, repeat=False)
# ani.save('x.gif')

# ti = 0
# while 1:
#     ax.cla()
#     ax.scatter(pos[:,0], pos[:,1], s = 5, c = colours[state])
#     ax.set(xlim=[0,scale], ylim=[0,scale])
#     plt.pause(0.005)

#     t.append(ti)
#     stat = count_stats(state)
#     sir.append(stat)

#     grid = update(pos, vel, f0, fc, fv)

#     update_diease(pos, grid, state, t2rec, pi, ri, T, scale, subdiv)
#     ti += 1

#     if not stat[1]:
#         break
    

plt.show()