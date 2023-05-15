from dataclasses import dataclass
from typing import List

@dataclass
class poly:
    coeffs: List[int]

    def __init__(self, inp: List[int] = [0]*10):
        if len(inp) > 10:
            raise ValueError("Polynomial can't have more than 10 coeffs")
        if len(inp) < 10:
            inp = [0]*(10-len(inp)) + inp
        self.coeffs = inp
