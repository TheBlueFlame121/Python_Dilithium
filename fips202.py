# Contains stuff from fips202.c
# Instead of implementing Keccak from scratch, a library is used

from Crypto.Hash import SHAKE128, SHAKE256, SHA3_256, SHA3_512

SHAKE128_RATE = 168
SHAKE256_RATE = 136
SHA3_256_RATE = 136
SHA3_512_RATE = 72


#################################################
# Name:        shake128_squeezeblocks
#
# Description: Squeeze step of SHAKE128 XOF. Squeezes full blocks of
#              SHAKE256_RATE bytes each. Can be called multiple times
#              to keep squeezing. Assumes next block has not yet been
#              started (state->pos = SHAKE256_RATE).
#
# Arguments:   - state: SHAKE128_XOF object
#              - int nblocks: requested number of blocks
#
# Returns the output of the XOF
##################################################
def shake128_squeezeblocks(nblocks:int, state:Crypto.Hash.SHAKE128.SHAKE128_XOF) -> List[int]:
    return list(state.read(nblocks*SHAKE128_RATE))


#################################################
# Name:        shake128
#
# Description: SHAKE128 XOF with non-incremental API
#
# Arguments:   - bytes input: input bytes for the XOF
#              - int output_len: requested output length
#
# Returns the output of the XOF
##################################################
def shake128(input:bytes, output_len: int) -> bytes:
    temp = SHAKE128.new()
    temp.update(input)
    return temp.read(output_len)


#################################################
# Name:        shake256_squeezeblocks
#
# Description: Squeeze step of SHAKE256 XOF. Squeezes full blocks of
#              SHAKE256_RATE bytes each. Can be called multiple times
#              to keep squeezing. Assumes next block has not yet been
#              started (state->pos = SHAKE256_RATE).
#
# Arguments:   - state: SHAKE256_XOF object
#              - int nblocks: requested number of blocks
#
# Returns the output of the XOF
##################################################
def shake256_squeezeblocks(nblocks:int, state:Crypto.Hash.SHAKE256.SHAKE256_XOF) -> List[int]:
    return list(state.read(nblocks*SHAKE256_RATE))


#################################################
# Name:        shake256
#
# Description: SHAKE256 XOF with non-incremental API
#
# Arguments:   - bytes input: input bytes for the XOF
#              - int output_len: requested output length
#
# Returns the output of the XOF
##################################################
def shake256(input:bytes, output_len: int) -> bytes:
    temp = SHAKE256.new()
    temp.update(input)
    return temp.read(output_len)


#################################################
# Name:        sha3_256
#
# Description: SHA3-256 with non-incremental API
#
# Arguments:   - bytes input: input to be hashed
#
# Returns the output of the hash function
##################################################
def sha3_256(input: bytes) -> bytes:
    temp = SHA3_256.new()
    temp.update(input)
    return temp.digest()


#################################################
# Name:        sha3_512
#
# Description: SHA3-512 with non-incremental API
#
# Arguments:   - bytes input: input to be hashed
#
# Returns the output of the hash function
##################################################
def sha3_512(input:bytes) -> bytes:
    temp = SHA3_512.new()
    temp.update(input)
    return temp.digest()
