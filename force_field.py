from typing import Any, Dict
import numpy as np
from numpy.random import uniform as unif

class forcefield:
    """
    A force field.
    """

    def __init__(self, sources: Dict[tuple, float]) -> None:
        if not isinstance(sources, dict):
            raise TypeError("sources must be a 'dict'")
        elif not min( map( lambda o: len(o) == 2, sources.keys() ) ):
            raise TypeError("keys must be iterables of length 2")
        elif not min( map( lambda o: np.ndim(o) == 0, sources.values() ) ):
            raise TypeError("values must be 'float'")
        
        self.src = np.asarray(list(sources.keys()))
        self.w   = np.asarray(list(sources.values()))

    def field(self, x: Any, y: Any) -> Any:
        """ 
        Potential field at a point.
        """
        ux, uy = 0.0, 0.0
        for (sx, sy), sq in zip(self.src, self.w):
            dx, dy = (sx - x), (sy - y)
            r = np.sqrt(dx*dx + dy*dy + 1e-8)
            u = sq * np.exp(-1.0 * r**2)

            ux += u * dx / r
            uy += u * dy / r
        return ux, uy

class agent:
    """
    A random walker.
    """

    def __init__(self, x: float, y: float, vx: float, vy: float, rinf: float = 2.0, lim: float = None) -> None:
        self.x, self.y   = x, y
        self.vx, self.vy = vx, vy
        self.status = 's'
        self.rinf   = rinf
        self.lim    = lim
        self.field  = None

    def setfield(self, field: forcefield) -> None:
        self.field = field

    def move(self, dt: float = 0.01) -> None:
        """ 
        move the agent.
        """
        self.x += self.vx * dt / 2.0
        self.y += self.vy * dt / 2.0

        if self.lim:
            self.x = self.x % self.lim
            self.y = self.y % self.lim

        ax, ay = 0.0, 0.0
        if self.field:
            ax, ay = self.field.field(self.x, self.y)

        self.vx += ax * dt
        self.vy += ay * dt

        self.x += self.vx * dt / 2.0
        self.y += self.vy * dt / 2.0

        if self.lim:
            self.x = self.x % self.lim
            self.y = self.y % self.lim
    

class community:
    """
    A community of random walkers.
    """

    def __init__(self, size: int, limit: float) -> None:
        self.size  = size
        self.limit = limit
        
        self.people = [ 
                            agent(
                                    x = unif(0., limit), y = unif(0, limit),
                                    vx = unif(-5, 5), vy = unif(-5, 5),
                                    lim = self.limit,
                                 ) for i in range(size)
                      ]

    def __getitem__(self, __i: int) -> agent:
        return self.people[__i]

    def move(self, dt: float = 0.01) -> None:
        for p in self.people:
            p.move(dt)

    def setfield(self, field: forcefield) -> None:
        for p in self.people:
            p.setfield(field)




        
src = { 
            (3.0, 3.0): 100.0, 
            (7.0, 7.0): 100.0,
            (3.0, 7.0): 100.0,
            (7.0, 3.0): 100.0,
      }
f = forcefield(src)
c = community(50, 10.0)
c.setfield(f)

import matplotlib.pyplot as plt
plt.style.use('ggplot')

x, y = np.mgrid[0:10:201j, 0:10:201j]
ux, uy = f.field(x, y)
u = np.sqrt(ux**2 + uy**2)

fig, ax = plt.subplots(figsize=[8,8])

for t in range(50):
    ax.cla()

    ax.pcolor(x, y, u, cmap = 'Reds', alpha = 1.0)
    # ax.plot(f.src[:,0], f.src[:,1], 'o', ms = 3, color = 'black')

    for p in c.people:
        ax.plot(p.x, p.y, 'o', ms = 3, color = 'black')

    ax.set_title(f't = {t}')
    ax.set(xlim = [0,10], ylim = [0,10])

    plt.pause(0.001)

    c.move()

plt.show()