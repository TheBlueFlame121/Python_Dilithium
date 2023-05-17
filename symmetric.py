from params import *
from aes256_ctr_drbg import *
from fips202 import *

if DILITHIUM_USE_AES:
    STREAM128_BLOCKBYTES = AES256CTR_BLOCKBYTES
    STREAM256_BLOCKBYTES = AES256CTR_BLOCKBYTES
else:
    STREAM128_BLOCKBYTES = SHAKE128_RATE
    STREAM256_BLOCKBYTES = SHAKE256_RATE
