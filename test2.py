import numpy as np

class community:
    """ a group of random walking agents """

    def __init__(self, n: int, v0: float = 1.0, lim: float = 100.0, rinf: float = 0.05, fcentre: list = ..., f0: float = 10.0, w: float = 0.05, dt: float = 0.1) -> None:
        if fcentre is ... :
            fcentre = [lim / 2.0, lim / 2.0]
            
        self.n       = n                   # number of agents
        self.lim     = lim                 # space limits
        self.rinf    = self.lim * rinf     # infection radius
        self.dt      = dt                  # time step for integration
        
        self.fcentre    = np.asarray(fcentre) # centre of the force field
        self.f0, self.w = f0, w               # force parameters

        # create random positions
        self.pos = np.random.uniform(0.0, self.lim, (self.n, 2))

        # create random velocities
        angle    = np.random.uniform(0.0, 2*np.pi, self.n)
        self.vel = v0 * np.array([np.cos(angle), np.sin(angle)]).T

        # infection status: all are susceptible now
        self.state = np.zeros(self.n) # 0 -> susceptible, 1 -> infected, 2 -> recovered

        self.soft = 0.1 # softening value

    def move(self) -> None:
        """ move the agents by one step. """
        # using leapfrog integration
        self.pos = self.pos + self.vel * self.dt * 0.5 # postion half-step

        # get the acceleration / force
        if self.f0:
            sep = self.pos - self.fcentre
            r   = np.sqrt(np.sum(sep**2, -1) + self.soft**2)
            acc = -self.f0 * np.exp(-(r / self.lim)**2 / self.w) / r
            acc = acc[:, None] * sep # acceleration = force / mass

            self.vel = self.vel + acc * self.dt # velocity step

        self.pos = self.pos + self.vel * self.dt * 0.5 # postion half-step

        # put reflection boundary conditions
        mask            = (self.pos > self.lim) | (self.pos < 0.0) # crossing the boundary
        self.vel[mask] *= -1




import matplotlib.pyplot as plt
import matplotlib.animation as anim
plt.style.use('ggplot')

com = community(n = 100, )


fig, ax = plt.subplots(figsize = [10,10])
ax.axis([0., com.lim, 0, com.lim])

for i in range(1000):
    ax.cla()
    ax.plot(com.pos[:,0], com.pos[:,1], 'o', ms = 3, color = 'black')
    com.move()
    ax.set(xlim=[0,com.lim], ylim=[0,com.lim])
    plt.pause(0.01)

plt.show()