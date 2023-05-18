from params import *
from fips202 import *


STREAM128_BLOCKBYTES = SHAKE128_RATE
STREAM256_BLOCKBYTES = SHAKE256_RATE

stream128_state = SHAKE128.new
stream256_state = SHAKE256.new


def shake128_stream_init(state: Crypto.Hash.SHAKE128.SHAKE128_XOF, seed: bytes, nonce: int):
    state.update(seed)
    state.update(nonce.to_bytes(2, 'little'))


def shake256_stream_init(state: Crypto.Hash.SHAKE256.SHAKE256_XOF, seed: bytes, nonce: int):
    state.update(seed)
    state.update(nonce.to_bytes(2, 'little'))


stream128_init = shake128_stream_init
stream128_squeezeblocks = shake128_squeezeblocks

stream256_init = shake256_stream_init
stream256_squeezeblocks = shake256_squeezeblocks
