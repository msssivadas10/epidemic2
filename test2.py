from itertools import product
import numpy as np

class community:
    """ a group of random walking agents """

    def __init__(self, n: int, v0: float = 1.0, lim: float = 100.0, rinf: float = 0.01, dt: float = 0.1) -> None:    
        self.n       = n                   # number of agents
        self.lim     = lim                 # space limits
        self.rinf    = self.lim * rinf     # infection radius
        self.pinf    = 0.0                 # infection probability
        self.trec    = 0.0                 # recovery time
        self.time = 0.0
        self.dt      = dt                  # time step for integration
        
        self.fcentre    = np.array([0.0, 0.0]) # centre of the force field
        self.vcenter    = np.array([0.0, 0.0]) # velocity of the force field
        self.f0, self.w = 0.0, 0.0             # force parameters
        self.soft     = 0.1                    # softening value

        self.e0, self.r0 = 1.0, 0.01 * self.lim # agent force parameters

        self.subdiv   = 10  # grid subdivisions
        self.cellsize = self.lim / self.subdiv

        if self.rinf > self.cellsize:
            raise ValueError("rinf cannot greater than cellsize")

        # create random positions
        self.pos = np.random.uniform(0.0, self.lim, (self.n, 2))

        # create random velocities
        angle    = np.random.uniform(0.0, 2*np.pi, self.n)
        self.vel = v0 * np.array([np.cos(angle), np.sin(angle)]).T

        # infection status: all are susceptible now
        # 0 -> susceptible, 1 -> infected, 2 -> recovered
        self.state = np.zeros(self.n, dtype = 'int') 
        self.t2rec = np.zeros(self.n) # time left to recover

        self.gridify()

        self.cforce_on   = False
        self.aforce_on   = False
        self.has_disease = False

    def config_cforce(self, f0: float = ... , w: float = ..., centre: tuple = ..., vcentre: tuple = ..., ) -> None:
        """ create a central force field """
        if centre is not ... :
            self.fcentre = self.lim * np.asarray(centre)
        if f0 is not ... :
            self.f0      = f0
        if w is not ... :
            self.w       = w
        if vcentre is not  ... :
            self.vcenter = np.asarray(vcentre)

    def config_aforce(self, e0: float = ... , r0: float = ... ) -> None:
        """ configure agent force """
        if e0 is not ... :
            self.e0 = e0
        if r0 is not ... :
            self.r0 = r0 * self.lim

    def config_disease(self, pinf: float = 0.2, trec: float = 10) -> None:
        """ configure the disease """
        self.pinf, self.trec = pinf, trec
        self.has_disease     = True

    def infect_random(self) -> None:
        """ infect a random agent. """
        if not self.has_disease:
            self.config_disease()

        i = np.random.randint(0, self.n)
        self.state[i] = 1
        self.t2rec[i] = self.trec

    def move(self) -> None:
        """ move the agents by one step. """
        # using leapfrog integration
        self.pos = self.pos + self.vel * self.dt * 0.5 # postion half-step

        # get the acceleration / force
        acc = 0.0

        if self.cforce_on: # acceleration due to the central force
            acc += self.central_acc(self.pos)
        if self.aforce_on: # acceleration due to agent forces
            acc += self.agent_acc(self.pos) 

        self.vel = self.vel * (1 + 0.05 * np.random.uniform(-1, 1, size = self.vel.shape)) + acc * self.dt # velocity step

        self.pos = self.pos + self.vel * self.dt * 0.5 # postion half-step

        # put reflection boundary conditions
        mask            = (self.pos  < 0) # crossing the left/bottom boundary
        self.pos[mask]  = 0
        self.vel[mask] *= -0.8

        mask            = (self.pos > self.lim) # crossing the right/top boundary
        self.pos[mask]  = self.lim
        self.vel[mask] *= -0.8

        # move force field
        if sum(self.vcenter**2) > 1e-8:
            self.fcentre = self.fcentre + self.vcenter * self.dt

             # put reflection boundary conditions
            mask            = (self.fcentre  < 0) # crossing the left/bottom boundary
            self.fcentre[mask]  = 0
            self.vcenter[mask] *= -1.0

            mask            = (self.fcentre > self.lim) # crossing the right/top boundary
            self.fcentre[mask]  = self.lim
            self.vcenter[mask] *= -1.0


        self.gridify() # re-grid

        self.update_state() # update state

        self.time += self.dt

    def central_acc(self, pos: float) -> float:
        """ acceleration due to central force """
        sep = pos - self.fcentre
        r   = np.sqrt(np.sum(sep**2, -1) + self.soft**2)
        acc = -self.f0 * np.exp(-(r / self.lim)**2 / self.w) / r
        acc = acc[:, None] * sep # acceleration = force / mass
        return acc

    def agent_acc(self, pos: float) -> float: # TODO
        """ acceleration due to agent force. """
        # return 0.0
        acc = np.zeros_like(pos)
        for p in pos:
            near = self.nearest(p)
            if not len(near):
                continue
            sep = p - self.pos[near,:]
            r   = np.sqrt(np.sum(sep**2, -1))

            # filter non-zero values
            mask         = (r > 1e-8)
            near, r, sep = near[mask], r[mask], sep[mask,:]

            acc_ = -self.e0 * np.exp(-(r / self.r0)**2) / r
            acc[near, :] += acc_[:, None] * sep 
        return acc

    def gridify(self) -> None:
        """ arrange the agents in a grid. """
        subdiv   = self.subdiv
        cellsize = self.cellsize

        cellpos = (self.pos // cellsize).astype(int)

        cells, inv = np.unique(cellpos, axis = 0, return_inverse = True)
        
        # grid = {(i,j): [] for i in range(subdiv) for j in range(subdiv)}
        grid = {}
        for i, (ci, cj) in enumerate(cells):
            grid[ci, cj] = np.where(inv == i)[0]

        self.grid = grid

    def nearest(self, pos: list) -> list:
        """ get the nearest agents to an agent. """
        def clip(x: float, a: float, b: float) -> float:
            """ clip value between a and b """
            return max(a, min(x, b))

        def intersect(cell: tuple, center: list) -> bool:
            """ check if a cell intersect with a cell. """
            r = self.cellsize

            x, y   = cell[0] * r, cell[1] * r
            near_x = clip(center[0], x, x + r)
            near_y = clip(center[1], y, y + r)
            return (center[0] - near_x)**2 + (center[1] - near_y)**2 < r**2

        i, j = (np.asarray(pos) // self.cellsize).astype(int)

        search_cells = product(
                                range(
                                        max(0, i-1),
                                        min(i+2, self.subdiv)
                                     ),
                                range(
                                        max(0, j-1),
                                        min(j+2, self.subdiv)
                                     )
                              )
        neighbours = []
        for cell in search_cells:
            if cell not in self.grid.keys():
                continue
            if not intersect(cell, pos):
                continue
            k     = self.grid[cell]
            pos2  = self.pos[k, :] # others in the cell
            
            
            select, = np.where(np.sqrt(np.sum((pos2 - pos)**2, axis = -1)) <= self.cellsize**2)
            neighbours = neighbours + list(k[select])
        return np.array(neighbours)

    def update_state(self) -> None:
        """ update the status. """
        # infect nearest heighbours
        inf = np.where(self.state == 1)[0] # infected agents
        if not len(inf):
            return

        for i in inf:
            near = self.nearest(self.pos[i,:])

            exp  = near[self.state[near] == 0] # nearset exposed
            prob = np.random.uniform(0, 1, exp.shape)

            to_infect = exp[(prob < self.pinf) & (np.sum((self.pos[i,:] - self.pos[exp,:])**2, -1) < self.rinf**2)]

            self.state[to_infect] = 1
            self.t2rec[to_infect] = self.trec

        # update recovery time
        self.t2rec[inf] -= self.dt 

        # recover some agents
        to_recover, = np.where(self.t2rec < 0.)
        self.t2rec[to_recover] = 0.0
        self.state[to_recover] = 2

    def switch_central_force(self, state: bool) -> None:
        """ switch on/off the central force. """
        self.cforce_on = state

    def switch_agent_force(self, state: bool) -> None:
        """ switch on/off the agent force. """
        self.aforce_on = state







import matplotlib.pyplot as plt
import matplotlib.animation as anim
plt.style.use('ggplot')

<<<<<<< HEAD
com = community(n = 100)
com.config_cforce(f0 = 10.0, w = 0.05, centre = [0.5, 0.5])
com.config_aforce(r0 = 0.1, e0 = 1)
=======
com = community(n = 2500, v0 = 1)
com.config_cforce(f0 = 10.0, w = 0.05, centre = [0.5, 0.5])#, vcentre = [1.0, 1.0])
com.config_aforce(r0 = 0.1, e0 = 1e-2)
com.config_disease(pinf = 0.2, trec = 10)
>>>>>>> 15cb2d8a23f0850748995a6173fcdcd4ee6c730f

com.infect_random()
# com.switch_central_force(1)
com.switch_agent_force(1)


fig, ax = plt.subplots(figsize = [10,10])
ax.axis([-5, com.lim+5, -5, com.lim+5])

colours = np.array(['#000000', '#ff0000', '#0000ff'])
# for i in range(1000):
while len(np.where(com.state == 1)[0]) > 0:
    ax.cla()
    col = colours[com.state]
    ax.scatter(com.pos[:,0], com.pos[:,1], s = 5, c = col)

    # ax.scatter(com.pos[:,0], com.pos[:,1], s = 5, c = 'black')

    # k = com.nearest([50, 50])

    # ax.scatter(com.pos[k,0], com.pos[k,1], s = 5, c = 'red')

    ax.plot([0, com.lim, com.lim, 0, 0], [0, 0, com.lim, com.lim, 0], lw = 2, color = 'gray')

    com.move()
    # ax.set(xlim=[0,com.lim], ylim=[0,com.lim])
    ax.set_title("t = {:.3f}".format(com.time))
    plt.pause(0.01)



# __x = 0
# for i in range(com.subdiv):
#     ax.plot([0., com.lim], [__x, __x], lw = 0.5, color = 'black')
#     ax.plot([__x, __x], [0., com.lim], lw = 0.5, color = 'black')
#     __x += com.cellsize


# for (ci, cj), i in com.grid.items():
#     if not i:
#         continue
#     c = '#0000ff'
#     if ci%2 ^ cj%2:
#         c = '#ff0000'
#     ax.plot(com.pos[i,0], com.pos[i,1], 'o', ms = 3, color = c)#, color = 'black')

# plt.show()