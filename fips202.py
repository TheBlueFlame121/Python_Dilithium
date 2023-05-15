# Contains stuff from fips202.c
# Instead of implementing Keccak from scratch, a library is used

from Crypto.Hash import SHAKE128, SHAKE256, SHA3_256, SHA3_512


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
