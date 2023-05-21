# Contains elements from sign.h and sign.py
from params import *
from packing import *
from polyvec import *
from poly import *
from os import urandom
from symmetric import *
from fips202 import *


#################################################
# Name:        crypto_sign_keypair
#
# Description: Generates public and private key.
#
# Arguments:   - List[int] pk: output public key (allocated array
#                              of CRYPTO_PUBLICKEYBYTES bytes)
#              - List[int] sk: output private key (allocated array
#                              of CRYPTO_SECRETKEYBYTES bytes)
#
# Returns 0 (success)
##################################################
def crypto_sign_keypair(pk:List[int], sk:List[int]) -> int:
    # seedbuf = [0]*(2*SEEDBYTES+CRHBYTES)
    # tr = [0]*SEEDBYTES
    mat = [polyvecl() for _ in range(K)]
    s1 = polyvecl()
    s1hat = polyvecl()
    s2 = polyveck()
    t1 = polyveck()
    t0 = polyveck()

    # Get randomness for rho, rhoprime and key
    seedbuf = shake256(urandom(SEEDBYTES), 2*SEEDBYTES + CRHBYTES)
    rho = list(seedbuf[:SEEDBYTES])
    rhoprime = list(seedbuf[SEEDBYTES:SEEDBYTES+CRHBYTES])
    key = list(seedbuf[-SEEDBYTES:])

    # Expand matrix
    polyvec_matrix_expand(mat, rho)

    # Sample short vectors s1 and s2
    polyvecl_uniform_eta(s1, rhoprime, 0)
    polyveck_uniform_eta(s2, rhoprime, L)

    # Matrix-vector multiplication
    for i in range(L):
        for j in range(N):
            s1hat.vec[i].coeffs[j] = s1.vec[i].coeffs[j]
    polyvecl_ntt(s1hat)
    polyvec_matrix_pointwise_montgomery(t1, mat, s1hat)
    polyveck_reduce(t1)
    polyveck_invntt_tomont(t1)

    # Add error vector s2
    polyveck_add(t1, t1, s2)

    # Extract t1 and write public key
    polyveck_caddq(t1)
    polyveck_power2round(t1, t0, t1)
    pack_pk(pk, rho, t1)

    # Compute H(rho, t1) and write secret key
    tr = list(shake256(bytes(pk), SEEDBYTES))
    pack_sk(sk, rho, tr, key, t0, s1, s2)

    return 0
