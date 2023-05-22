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


#################################################
# Name:        crypto_sign_signature
#
# Description: Computes signature.
#
# Arguments:   - List[int] sig:   output signature (of length CRYPTO_BYTES)
#              - int siglen:      output length of signature (UNUSED)
#              - List[int] m:     message to be signed
#              - int mlen:        length of message (UNUSED)
#              - List[int] sk:    bit-packed secret key
#
# Returns 0 (success)
##################################################
def crypto_sign_signature(sig:List[int], siglen:int, m:List[int], mlen:int, sk:List[int]) -> int:
    # seedbuf = [0]*(3*SEEDBYTES + 2*CRHBYTES)
    rho, tr, key = ([0]*SEEDBYTES for _ in range(3))
    # mu, rhoprime = ([0]*CRHBYTES for _ in range(2))

    nonce = 0
    mat = [polyvecl() for _ in range(K)]
    s1 = polyvecl()
    y = polyvecl()
    z = polyvecl()
    t0 = polyveck()
    s2 = polyveck()
    w1 = polyveck()
    w0 = polyveck()
    h = polyveck()
    cp = poly()
    # state = stream256_state()

    unpack_sk(rho, tr, key, t0, s1, s2, sk)

    # Compute CRH(tr, msg)
    state = stream256_state()
    state.update(bytes(tr))
    state.update(bytes(m))
    mu = list(state.read(CRHBYTES))

    rhoprime = list(shake256(bytes(key+mu), CRHBYTES))

    # Expand matrix and transform vectors
    polyvec_matrix_expand(mat, rho)
    polyvecl_ntt(s1)
    polyveck_ntt(s2)
    polyveck_ntt(t0)

    while True:
        # Sample intermediate vector y
        polyvecl_uniform_gamma1(y, rhoprime, nonce)
        nonce += 1

        # Matrix-vector multiplication
        for i in range(L):
            for j in range(N):
                z.vec[i].coeffs[j] = y.vec[i].coeffs[j]
        polyvecl_ntt(z)
        polyvec_matrix_pointwise_montgomery(w1, mat, z)
        polyveck_reduce(w1)
        polyveck_invntt_tomont(w1)

        # Decompose w and call the random oracle
        polyveck_caddq(w1)
        polyveck_decompose(w1, w0, w1)
        polyveck_pack_w1(sig, w1)

        state = stream256_state()
        state.update(bytes(mu))
        state.update(bytes(sig[:K*POLYW1_PACKEDBYTES]))
        temp = list(state.read(SEEDBYTES))
        for i in range(SEEDBYTES):
            sig[i] = temp[i]
        poly_challenge(cp, sig[:SEEDBYTES])
        poly_ntt(cp)

        # Compute z, reject if it reveals secret
        polyvecl_pointwise_poly_montgomery(z, cp, s1)
        polyvecl_invntt_tomont(z)
        polyvecl_add(z, z, y)
        polyvecl_reduce(z)
        if polyvecl_chknorm(z, GAMMA1-BETA):
            continue

        # Check that subtracting cs2 does not change high bits of w and low bits
        # do not reveal secret information
        polyveck_pointwise_poly_montgomery(h, cp, s2)
        polyveck_invntt_tomont(h)
        polyveck_sub(w0, w0, h)
        polyveck_reduce(w0)
        if polyveck_chknorm(w0, GAMMA2 - BETA):
            continue

        # Compute hints for w1
        polyveck_pointwise_poly_montgomery(h, cp, t0)
        polyveck_invntt_tomont(h)
        polyveck_reduce(h)
        if polyveck_chknorm(h, GAMMA2):
            continue

        polyveck_add(w0, w0, h)
        n = polyveck_make_hint(h, w0, w1)
        if n > OMEGA:
            continue

        # Write signature
        pack_sig(sig, sig, z, h)

        break
    return 0


#################################################
# Name:        crypto_sign
#
# Description: Compute signed message.
#
# Arguments:   - List[int] sm: output signed message (allocated
#                              array with CRYPTO_BYTES + mlen bytes),
#                              can be equal to m
#              - int smlen:    output length of signed
#                              message (UNUSED)
#              - List[int] m:  message to be signed
#              - int mlen:     length of message (UNUSED)
#              - List[int] sk: bit-packed secret key
#
# Returns 0 (success)
##################################################
def crypto_sign(sm:List[int], smlen:int, m:List[int], mlen:int, sk:List[int]) -> int:
    for i in range(len(m)):
        sm[CRYPTO_BYTES + len(m) - 1 -i] = m[len(m) - 1 - i]
    crypto_sign_signature(sm, smlen, m, mlen, sk)
    return 0


#################################################
# Name:        crypto_sign_verify
#
# Description: Verifies signature.
#
# Arguments:   - List[int] sig: input signature
#              - int siglen:    length of signature (UNUSED)
#              - List[int] m:   message
#              - int mlen:      length of message (UNUSED)
#              - List[int] pk:  bit-packed public key
#
# Returns 0 if signature could be verified correctly and -1 otherwise
##################################################
def crypto_sign_verify(sig: List[int], siglen:int, m:List[int], mlen:int, pk:List[int]) -> int:
    buf = [0]*(K*POLYW1_PACKEDBYTES)
    rho = [0]*SEEDBYTES
    # mu = [0]*CRHBYTES
    c = [0]*SEEDBYTES
    # c2 = [0]*SEEDBYTES
    cp = poly()
    mat = [polyvecl() for _ in range(K)]
    z = polyvecl()
    t1 = polyveck()
    w1 = polyveck()
    h = polyveck()
    # state = stream256_state()

    if len(sig) != CRYPTO_BYTES:
        return -1

    unpack_pk(rho, t1, pk)
    if unpack_sig(c, z, h, sig):
        return -1
    if polyvecl_chknorm(z, GAMMA1 - BETA):
        return -1

    # Compute CRH(H(rho, t1), msg)
    mu = list(shake256(bytes(pk), SEEDBYTES))
    state = stream256_state()
    state.update(bytes(mu))
    state.update(bytes(m))
    mu = list(state.read(CRHBYTES))

    # Matrix-vector multiplication; compute Az - c2^dt1
    poly_challenge(cp, c)
    polyvec_matrix_expand(mat, rho)

    polyvecl_ntt(z)
    polyvec_matrix_pointwise_montgomery(w1, mat, z)

    poly_ntt(cp)
    polyveck_shiftl(t1)
    polyveck_ntt(t1)
    polyveck_pointwise_poly_montgomery(t1, cp, t1)

    polyveck_sub(w1, w1, t1)
    polyveck_reduce(w1)
    polyveck_invntt_tomont(w1)

    # Reconstruct w1
    polyveck_caddq(w1)
    polyveck_use_hint(w1, w1, h)
    polyveck_pack_w1(buf, w1)

    # Call random oracle and verify challenge
    state = stream256_state()
    state.update(bytes(mu))
    state.update(bytes(buf))
    c2 = state.read(SEEDBYTES)
    for i in range(SEEDBYTES):
        if c[i] != c2[i]:
            return -1
    
    return 0
