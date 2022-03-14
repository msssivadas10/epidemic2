from typing import Any, Dict
import numpy as np

class forcefield:
    """
    A force field.
    """

    def __init__(self, sources: Dict[tuple, float]) -> None:
        if not isinstance(sources, dict):
            raise TypeError("sources must be a 'dict'")
        elif min( map( lambda o: len(o) == 2, sources.keys() ) ):
            raise TypeError("keys must be iterables of length 2")
        elif min( map( lambda o: isinstance(o, float), sources.values() ) ):
            raise TypeError("values must be 'float'")
        
        self.src = np.asarray(list(sources.keys()))
        self.w   = np.asarray(list(sources.values()))

        
        
