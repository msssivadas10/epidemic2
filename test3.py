from itertools import product
from typing import Any
import numpy as np
import numpy.random as rnd


class Disease:
    __slots__ = 'p', 'trec', 

    def __init__(self, p: float, trec: float) -> None:
        if p < 0.0 or p > 1.0:
            raise ValueError("p (transmission probability) must be a value in [0, 1]")
        if trec < 0.0:
            raise ValueError("trec must be a positive value")
        self.p, self.trec = p, trec

    def __repr__(self) -> str:
        return f"Disease(p={self.p}, trec={self.trec})"


class ForceField:
    __slots__ = 'x', 'v', 'str', 'rad', 'soft', 'dt', 'scale', 'state', 

    def __init__(self, x: Any, v: Any = [0., 0.], _str: float = 1.0, rad: float = 1.0, soft: float = 1e-3) -> None:
        self.x, self.v = np.asarray(x), np.asarray(v)
        
        # force parameters:
        self.str  = _str  # force strength
        self.rad  = rad   # radius parameter set the extend of force
        self.soft = soft  # softening parameter

        # other settings:
        self.dt    = ... # time step
        self.scale = ... # space scaling
        self.state = True

    def setup(self, dt: float, scale: float) -> None:
        self.dt, self.scale = dt, scale
        
        self.x   *= self.scale
        self.rad *= self.scale

    def a(self, x: Any) -> Any:
        if not self.state:
            return 0.0
        x  = np.asarray(x)
        dx = self.x - x
        r  = np.sqrt(dx[:,0]**2 + dx[:,1]**2 + self.soft**2)
        # f  = self.str * np.exp(-r / self.rad) / r
        srp = (self.rad / r)**3
        f = 24 * self.str * srp * (2 * srp - 1)
        return f[:, None] * dx

    def move(self) -> None:
        if self.dt is ... :
            raise ValueError("force is not setup")
        self.x += self.v * self.dt

    def switch(self, state: bool) -> None:
        self.state = bool(state)

class Particles:
    __slots__ = (
                    'n', 'dt', 'T', 'ri', 'x', 'v', 'status', 't2rec', 
                    'force', 'disease', 'grid', 'inter', 'fc_affect', 'fi_affect'
                )
    scale     = 100.0
    subdiv    = 10

    def __init__(self, n: int, ri: float = 0.1, T: float = 1.0, dt: float = 0.01) -> None:
        self.n  = n
        self.dt = dt
        self.T  = 100.0 * T
        self.ri = ri * self.scale

        self.x = rnd.uniform(size = (self.n, 2)) * self.scale

        self.v      = np.zeros_like(self.x)
        angle       = rnd.uniform(size = self.n) * (2*np.pi)
        self.v[:,0] = np.cos(angle)
        self.v[:,1] = np.sin(angle)
        self.v     -= np.mean(self.v, 0)
        self.v     *= np.sqrt(1.5 * self.T / np.mean(self.v**2))

        self.status = np.zeros(self.n, 'int')
        self.t2rec  = np.zeros(self.n)

        self.fc_affect = np.ones(self.n)
        self.fi_affect = np.ones(self.n)

        self.force   = None
        self.disease = None
        self.inter   = False

        self.createGrid()

    def move(self) -> None:
        def xstep():
            self.x += self.v * self.dt * 0.5
            
            mask    = (self.x < 0.0)
            self.x[mask] = 0.0; self.v[mask] *= -rnd.uniform(0.5, 1.0)

            mask    = (self.x > self.scale)
            self.x[mask] = self.scale; self.v[mask] *= -rnd.uniform(0.5, 1.0)

        xstep()

        # self.v -= 0.01 *self.v * self.dt
        self.v += (self.af + self.ai) * self.dt + self.v * rnd.uniform(-1, 1, self.v.shape) * 0.05
        # f1 = np.sqrt(1.5 * self.T / np.mean(self.v**2 * self.fi_affect[:,None]))
        f1 = np.sqrt(1.5 * self.T / np.mean(self.v**2))
        # self.v[self.fi_affect.astype('bool')]  *= f1
        self.v *= f1


        xstep()

        self.createGrid()

        if self.force:
            self.force.move()

        # disease dynamics:
        _inf  = np.where(self.status == 1)[0]
        self.t2rec[_inf] -= self.dt

        for k in _inf:
            near  = self.nearest(self.x[k,:], self.ri)
            if not len(near):
                continue
            near  = near[self.status[near] == 0]
            _2inf = near[rnd.uniform(size = near.shape) < self.disease.p]
            self.status[_2inf] = 1
            self.t2rec[_2inf]  = self.disease.trec
        
        _2rec = (self.t2rec < 0.0)
        self.t2rec[_2rec] = 0.0; self.status[_2rec] = 2

    @property
    def af(self) -> Any:
        if self.force:
            return self.force.a(self.x) * self.fc_affect[:,None]
        return 0.0

    @property
    def ai(self) -> Any:
        if self.inter:
            s  = self.scale / np.sqrt(self.n)
            rc = 2**(1/6) * s
            a  = np.zeros_like(self.x)
            for k1, xy1 in enumerate(self.x):
                k2       = self.nearest(xy1, rc)
                k2       = k2[k2 > k1]
                if not len(k2):
                    continue
                xy2      = self.x[k2, :]
                dx       = xy2 - xy1
                r        = np.sqrt(dx[:,0]**2 + dx[:,1]**2) + 1e-6
                srp      = (s / r)**3  # (s/r)**p with p = 3
                u        = 24 * srp * (2 * srp - 1) / r 
                u[u>100] = 100
                u        = u[:,None] * self.fi_affect[k2, None] * dx
                a[k1,:] -= np.sum(u, 0) 
                a[k2,:] += u 
            return a #* self.fi_affect[:,None]
        return 0.0

    def addForce(self, f: ForceField) -> None:
        if not isinstance(f, ForceField):
            raise TypeError("f must be a 'ForceField'")
        self.force = f
        self.force.setup(self.dt, self.scale)

    def infectRandom(self, d: Disease) -> None:
        if not isinstance(d, Disease):
            raise TypeError("d must be a 'Disease'")
        self.disease = d 

        # infect a random particle
        i = rnd.randint(self.n)
        self.status[i] = 1
        self.t2rec[i]  = self.disease.trec

    def createGrid(self) -> None:
        cellsize = self.scale / self.subdiv
        cx       = (self.x // cellsize).astype('int')
        grid     = {}
        for k, ij in enumerate(cx):
            ij = tuple(ij)
            if ij not in grid:
                grid[ij] = []
            grid[ij].append(k)
        self.grid = grid

    def nearest(self, xy: Any, dist: float) -> list:
        cellsize = self.scale / self.subdiv
        i, j     = (np.array(xy) // cellsize).astype('int')
        steps    = -int(-dist // cellsize)
        near     = []
        for _ij in product( range(i-steps, i+steps+1), range(j-steps, j+steps+1) ):
            if _ij not in self.grid.keys():
                continue
            k = np.array(self.grid[_ij])
            k = k[(self.x[k,0] - xy[0])**2 + (self.x[k,1] - xy[1])**2 < dist**2]
            near = near + list(k)
        return np.array(near)

    def switchInteractions(self, state: bool, val :float = ..., ) -> None:
        self.inter = bool(state)
        if val is not  ... :
            if val < 0.0 or val > 1.0:
                raise ValueError("'val' must be a number between 0 and 1")
            self.fi_affect = (rnd.uniform(size = self.n) < val).astype('int')

    def switchForce(self, state: bool, val :float = ..., ) -> None:
        self.force.switch(state)
        if val is not ... :
            if val < 0.0 or val > 1.0:
                raise ValueError("'val' must be a number between 0 and 1")
            self.fc_affect = (rnd.uniform(size = self.n) < val).astype('int')

        


g = Particles(100, 0.01)
f = ForceField([0.5, 0.5], _str = 100.0, rad = 0.5)
g.addForce(f)
d = Disease(0.2, 200*g.dt)
# g.infectRandom(d)
# g.switchForce(1)
# print(g.ri)
g.switchInteractions(1, 0.9)



import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize = [6, 6])

col = np.array([[0.0, 0.2, 0.4], [0.9, 0.0, 0.1], [0.0, 0.8, 0.5]])
for i in range(500):
    ax.cla()
    ax.set_title(f"steps = {i}")
    ax.scatter(g.x[:,0], g.x[:,1], s = 5, c = col[g.status,:])
    # for k in np.where(g.status==1)[0]:
    #     ax.add_patch(plt.Circle(g.x[k,:], g.ri, color=col[1,:], fill=0))
    ax.set(xlim = [0, g.scale], ylim = [0, g.scale])
    plt.pause(0.0001)
    g.move()

    # if i == 100:
    #     g.switchInteractions(1, 0.9)

plt.show()

# print(g.grid)

