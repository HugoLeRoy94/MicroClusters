import numpy as np
import BOX
from typing import List as List, Tuple

class Object:
    def __init__(self,position:Tuple[int, int, int]):
        self._position:Tuple[int,int,int] = position

    @property
    def position(self) -> Tuple[int,int,int]:
        return self._position

    @position.setter
    def position(self, value: Tuple[int,int,int]) -> None:
        if not isinstance(value, tuple):
            raise TypeError(f"Expected 'tuple', got {type(value).__name__}")
        self._position = value
    def isempty(self)->int:
        return False
    def Index(self)->int:
        return 0
class Empty(Object):
    def __init__(self,position:Tuple[int, int, int]):
        super().__init__(position)
        return
    def isempty(self)->bool:
        return True
    def Index(self)->int:
        return 0
class RNA(Object):
    def __init__(self,
                 length:int,
                 x:List[int],
                 y:List[int],
                 z:List[int])-> None:
        return
    def Index(self)->int:
        return 2

class DHH1(Object):
    def __init__(self,position:Tuple[int, int, int])->None:
        super().__init__(position)

    def Index(self)->int:
        return 1
    def get_site_to_exchange(self,box:BOX)->Tuple:
        while True:
            x = np.random.randint(0,box.size)
            y = np.random.randint(0,box.size)
            z = np.random.randint(0,box.size)
            if (x,y,z) != self._position:
                break
        return (x,y,z)