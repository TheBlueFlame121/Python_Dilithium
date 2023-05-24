# Contains elements from params.h 

# from config import *
from typing import List

class Parameters:
    SEEDBYTES = 32
    CRHBYTES = 64
    N = 256
    Q = 8380417
    D = 13
    ROOT_OF_UNITY = 1753

    def __init__(self, mode:int):
        assert mode in [2, 3, 5]
        self.DILITHIUM_MODE = mode
        if self.DILITHIUM_MODE == 2:
            self.K = 4
            self.L = 4
            self.ETA = 2
            self.TAU = 39
            self.BETA = 78
            self.GAMMA1 = (1 << 17)
            self.GAMMA2 = ((self.Q-1)//88)
            self.OMEGA = 80
        elif self.DILITHIUM_MODE == 3:
            self.K = 6
            self.L = 5
            self.ETA = 4
            self.TAU = 49
            self.BETA = 196
            self.GAMMA1 = (1 << 19)
            self.GAMMA2 = ((self.Q-1)//32)
            self.OMEGA = 55
        elif self.DILITHIUM_MODE == 5:
            self.K = 8
            self.L = 7
            self.ETA = 2
            self.TAU = 60
            self.BETA = 120
            self.GAMMA1 = (1 << 19)
            self.GAMMA2 = ((self.Q-1)//32)
            self.OMEGA = 75

        self.POLYT1_PACKEDBYTES = 320
        self.POLYT0_PACKEDBYTES = 416
        self.POLYVECH_PACKEDBYTES = self.OMEGA + self.K

        if self.GAMMA1 == (1 << 17):
            self.POLYZ_PACKEDBYTES = 576
        elif self.GAMMA1 == (1 << 19):
            self.POLYZ_PACKEDBYTES = 640


        if self.GAMMA2 == (self.Q-1)/88:
            self.POLYW1_PACKEDBYTES = 192
        elif self.GAMMA2 == (self.Q-1)/32:
            self.POLYW1_PACKEDBYTES = 128

        if self.ETA == 2:
            self.POLYETA_PACKEDBYTES = 96
        elif self.ETA == 4:
            self.POLYETA_PACKEDBYTES = 128

        self.CRYPTO_PUBLICKEYBYTES = (self.SEEDBYTES + self.K*self.POLYT1_PACKEDBYTES)
        self.CRYPTO_SECRETKEYBYTES = (3*self.SEEDBYTES \
                                       + self.L*self.POLYETA_PACKEDBYTES \
                                       + self.K*self.POLYETA_PACKEDBYTES \
                                       + self.K*self.POLYT0_PACKEDBYTES)
        self.CRYPTO_BYTES = (self.SEEDBYTES + self.L*self.POLYZ_PACKEDBYTES + self.POLYVECH_PACKEDBYTES)

        STREAM256_BLOCKBYTES = 136
        if self.ETA == 2:
            self.POLY_UNIFORM_ETA_NBLOCKS = ((136 + STREAM256_BLOCKBYTES - 1)//STREAM256_BLOCKBYTES)
        elif self.ETA == 4:
            self.POLY_UNIFORM_ETA_NBLOCKS = ((227 + STREAM256_BLOCKBYTES - 1)//STREAM256_BLOCKBYTES)

        self.POLY_UNIFORM_GAMMA1_NBLOCKS = ((self.POLYZ_PACKEDBYTES + STREAM256_BLOCKBYTES - 1)//STREAM256_BLOCKBYTES)

    def set_mode(self, mode:int):
        assert mode in [2, 3, 5]
        self.__init__(mode)

g = Parameters(2)
