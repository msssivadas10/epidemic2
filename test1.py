import numpy as np
import numpy.random as rnd
import matplotlib.pyplot as plt
plt.style.use('ggplot')
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

def a_attr(pos, scale, f0, w = 0.05):
    r  = pos - 50.0
    _r = np.sqrt(r[:,:1]**2 + r[:,1:]**2 + 0.1)

    acc = -f0 * np.exp(-(_r/scale)**2 / w) * r / _r
    return acc

def update(pos, vel):
    pos = pos + vel * dt 

    acc = a_attr(pos, scale, f0, 0.05)
    vel = vel + acc * dt

    # pos = pos + vel * dt * 0.5

    mask = (pos < 0)
    pos[mask] = -pos[mask]
    vel[mask] = -vel[mask]

    mask = (pos > scale)
    pos[mask] = 2 * scale - pos[mask]
    vel[mask] = -vel[mask]

    # vel = vel / (vel[:,:1]**2 + vel[:,1:]**2)**0.5

    # vel = np.clip(vel, -1, 1)
    grid = make_grid(pos, subdiv, scale)
    return pos, vel, grid

def init(n, scale, subdiv):
    pos   = rnd.uniform(0, 1, (n, 2)) * scale

    theta = rnd.uniform(0, 2*np.pi, n)
    vel   =  np.array([np.cos(theta), np.sin(theta)]).T
    
    grid  = make_grid(pos, subdiv, scale)
    return pos, vel, grid

def init_disease(n, pi, ri, T):
    state = np.zeros(n, 'int')
    t2rec = np.zeros(n)

    # infect random:
    i = rnd.randint(n)
    state[i] = 1
    t2rec[i] = T
    return state, t2rec

def update_diease():
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

    


scale  = 100
n      = 100
dt     = 0.1
f0, fc = 0, np.array([50, 50])
subdiv = 10

pi, ri, T = 0.1, 0.1 * scale, 50 * dt

pos, vel, grid = init(n, scale, subdiv)

state, t2rec = init_disease(n, pi, ri, T)


# nn = get_nearest(grid, pos, [50,50], 10, scale, subdiv)
# print(nn)
colours = np.array(['black', 'C0', 'gray'])
fig, ax = plt.subplots(figsize=[6,6])

while 1:
    ax.cla()
    ax.scatter(pos[:,0], pos[:,1], s = 10, c = colours[state])
    ax.set(xlim=[0,scale], ylim=[0,scale])
    plt.pause(0.005)

    pos, vel, grid = update(pos, vel)

    
    

plt.show()