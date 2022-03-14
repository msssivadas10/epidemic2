from typing import Any, Dict
import numpy as np

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

        
src = { tuple(np.random.uniform(0.0, 10.0, 2)): np.random.uniform(0, 1) for i in range(10) }
f = forcefield(src)

import matplotlib.pyplot as plt
plt.style.use('ggplot')

plt.figure()
plt.plot(f.src[:,0], f.src[:,1], 'o', ms = 3, color = 'black')
plt.show()