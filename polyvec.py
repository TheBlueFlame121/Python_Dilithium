# Contains elements from polyvec.h and polyvec.c
from params import *
from poly import *


class polyvecl:
    vec: List[poly]

    def __init__(self, inp: list[poly] = None):
        if inp is None:
            inp = [poly() for _ in range(L)]
        if len(inp) > L:
            raise ValueError("Polynomial Vector L can't have more than L polys")
        if len(inp) < L:
            inp = inp + [poly() for _ in range(L-len(inp))].copy()
        self.vec = inp


class polyveck:
    vec: List[poly]

    def __init__(self, inp: list[poly] = None):
        if inp is None:
            inp = [poly() for _ in range(K)]
        if len(inp) > K:
            raise ValueError("Polynomial Vector K can't have more than K polys")
        if len(inp) < K:
            inp = inp + [poly() for _ in range(K-len(inp))].copy()
        self.vec = inp

#################################################
# Name:        expand_mat
#
# Description: Implementation of ExpandA. Generates matrix A with uniformly
#              random coefficients a_{i,j} by performing rejection
#              sampling on the output stream of SHAKE128(rho|j|i)
#              or AES256CTR(rho,j|i).
#
# Arguments:   - polyvecl mat[K]: output matrix
#              - List[int] rho: byte array containing seed rho
##################################################
def polyvec_matrix_expand(mat: List[polyvecl], rho: List[int]):
    for i in range(K):
        for j in range(L):
            poly_uniform(mat[i].vec[j], bytes(rho), (i<<8) + j)


def polyvec_matrix_pointwise_montgomery(t:polyveck, mat:List[polyvecl], v:polyvecl):
    for i in range(K):
        polyvecl_pointwise_acc_montgomery(t.vec[i], mat[i], v)


##############################################################
############ Vectors of polynomials of length L ##############
##############################################################

def polyvecl_uniform_eta(v:polyvecl, seed:List[int], nonce:int):
    for i in range(L):
        poly_uniform_eta(v.vec[i], seed, nonce)
        nonce+=1


def polyvecl_uniform_gamma1(v: polyvecl, seed:List[int], nonce:int):
    for i in range(L):
        poly_uniform_gamma1(v.vec[i], seed, L*nonce+i)


def polyvecl_reduce(v:polyvecl):
    for i in range(L):
        poly_reduce(v.vec[i])


#################################################
# Name:        polyvecl_add
#
# Description: Add vectors of polynomials of length L.
#              No modular reduction is performed.
#
# Arguments:   - polyvecl w: output vector
#              - polyvecl u: first summand
#              - polyvecl v: second summand
##################################################
def polyvecl_add(w:polyvecl, u:polyvecl, v:polyvecl):
    for i in range(L):
        poly_add(w.vec[i], u.vec[i], v.vec[i])


#################################################
# Name:        polyvecl_ntt
#
# Description: Forward NTT of all polynomials in vector of length L. Output
#              coefficients can be up to 16*Q larger than input coefficients.
#
# Arguments:   - polyvecl v: input/output vector
##################################################
def polyvecl_ntt(v:polyvecl):
    for i in range(L):
        poly_ntt(v.vec[i])


#################################################
# Name:        polyveck_invntt_tomont
#
# Description: Inverse NTT and multiplication by 2^{32} of polynomials
#              in vector of length L. Input coefficients need to be less
#              than 2*Q.
#
# Arguments:   - polyvecl *v: pointer to input/output vector
##################################################
def polyvecl_invntt_tomont(v:polyvecl):
    for i in range(L):
        poly_invntt_tomont(v.vec[i])


def polyvecl_pointwise_poly_montgomery(r:polyvecl, a:poly, v:polyvecl):
    for i in range(L):
        poly_pointwise_montgomery(r.vec[i], a, v.vec[i])


#################################################
# Name:        polyvecl_pointwise_acc_montgomery
#
# Description: Pointwise multiply vectors of polynomials of length L, multiply
#              resulting vector by 2^{-32} and add (accumulate) polynomials
#              in it. Input/output vectors are in NTT domain representation.
#
# Arguments:   - poly w: output polynomial
#              - polyvecl u: first input vector
#              - polyvecl v: second input vector
##################################################
def polyvecl_pointwise_acc_montgomery(w:poly, u:polyvecl, v:polyvecl):
    t = poly()
    poly_pointwise_montgomery(w, u.vec[0], v.vec[0])
    for i in range(1, L):
        poly_pointwise_montgomery(t, u.vec[0], v.vec[0])
        poly_add(w, w, t)


#################################################
# Name:        polyvecl_chknorm
#
# Description: Check infinity norm of polynomials in vector of length L.
#              Assumes input polyvecl to be reduced by polyvecl_reduce().
#
# Arguments:   - polyvecl v: vector
#              - int B: norm bound
#
# Returns 0 if norm of all polynomials is strictly smaller than B <= (Q-1)/8
# and 1 otherwise.
##################################################
def polyvecl_chknorm(v:polyvecl, bound:int) -> int:
    for i in range(L):
        if poly_chknorm(v.vec[i], bound):
            return 1

    return 0


##############################################################
############ Vectors of polynomials of length K ##############
##############################################################


def polyveck_uniform_eta(v:polyveck, seed:List[int], nonce:int):
    for i in range(K):
        poly_uniform_eta(v.vec[i], seed, nonce)
        nonce+=1


#################################################
# Name:        polyveck_reduce
#
# Description: Reduce coefficients of polynomials in vector of length K
#              to representatives in [-6283009,6283007].
#
# Arguments:   - polyveck v: input/output vector
##################################################
def polyveck_reduce(v:polyveck):
    for i in range(K):
        poly_reduce(v.vec[i])


#################################################
# Name:        polyveck_caddq
#
# Description: For all coefficients of polynomials in vector of length K
#              add Q if coefficient is negative.
#
# Arguments:   - polyveck v: input/output vector
##################################################
def polyveck_caddq(v:polyveck):
    for i in range(K):
        poly_caddq(v.vec[i])


#################################################
# Name:        polyveck_add
#
# Description: Add vectors of polynomials of length K.
#              No modular reduction is performed.
#
# Arguments:   - polyveck w: output vector
#              - polyveck u: first summand
#              - polyveck v: second summand
##################################################
def polyveck_add(w:polyveck, u:polyveck, v:polyveck):
    for i in range(K):
        poly_add(w.vec[i], u.vec[i], v.vec[i])


#################################################
# Name:        polyveck_sub
#
# Description: Subtract vectors of polynomials of length K.
#              No modular reduction is performed.
#
# Arguments:   - polyveck w: output vector
#              - polyveck u: first input vector
#              - polyveck v: second input vector to be
#                                   subtracted from first input vector
##################################################
def polyveck_sub(w:polyveck, u:polyveck, v:polyveck):
    for i in range(K):
        poly_sub(w.vec[i], u.vec[i], v.vec[i])


#################################################
# Name:        polyveck_shiftl
#
# Description: Multiply vector of polynomials of Length K by 2^D without modular
#              reduction. Assumes input coefficients to be less than 2^{31-D}.
#
# Arguments:   - polyveck v: input/output vector
##################################################
def polyveck_shiftl(v:polyveck):
    for i in range(K):
        poly_shiftl(v.vec[i])


#################################################
# Name:        polyveck_ntt
#
# Description: Forward NTT of all polynomials in vector of length K. Output
#              coefficients can be up to 16*Q larger than input coefficients.
#
# Arguments:   - polyveck v: input/output vector
##################################################
def polyveck_ntt(v:polyveck):
    for i in range(K):
        poly_ntt(v.vec[i])


#################################################
# Name:        polyveck_invntt_tomont
#
# Description: Inverse NTT and multiplication by 2^{32} of polynomials
#              in vector of length K. Input coefficients need to be less
#              than 2*Q.
#
# Arguments:   - polyveck v: input/output vector
##################################################
def polyveck_invntt_tomont(v:polyveck):
    for i in range(K):
        poly_invntt_tomont(v.vec[i])


def polyveck_pointwise_poly_montgomery(r:polyveck, a:poly, v:polyveck):
    for i in range(K):
        poly_pointwise_montgomery(r.vec[i], a, v.vec[i])


#################################################
# Name:        polyveck_chknorm
#
# Description: Check infinity norm of polynomials in vector of length K.
#              Assumes input polyveck to be reduced by polyveck_reduce().
#
# Arguments:   - polyveck v: pointer to vector
#              - int B: norm bound
#
# Returns 0 if norm of all polynomials are strictly smaller than B <= (Q-1)/8
# and 1 otherwise.
##################################################
def polyveck_chknorm(v:polyveck, bound:int) -> int:
    for i in range(K):
        if poly_chknorm(v.vec[i], bound):
            return 1

    return 0


#################################################
# Name:        polyveck_power2round
#
# Description: For all coefficients a of polynomials in vector of length K,
#              compute a0, a1 such that a mod^+ Q = a1*2^D + a0
#              with -2^{D-1} < a0 <= 2^{D-1}. Assumes coefficients to be
#              standard representatives.
#
# Arguments:   - polyveck v1: output vector of polynomials with
#                              coefficients a1
#              - polyveck v0: output vector of polynomials with
#                              coefficients a0
#              - polyveck v: input vector
##################################################
def polyveck_power2round(v1:polyveck, v0:polyveck, v:polyveck):
    for i in range(K):
        poly_power2round(v1.vec[i], v0.vec[i], v.vec[i])


#################################################
# Name:        polyveck_decompose
#
# Description: For all coefficients a of polynomials in vector of length K,
#              compute high and low bits a0, a1 such a mod^+ Q = a1*ALPHA + a0
#              with -ALPHA/2 < a0 <= ALPHA/2 except a1 = (Q-1)/ALPHA where we
#              set a1 = 0 and -ALPHA/2 <= a0 = a mod Q - Q < 0.
#              Assumes coefficients to be standard representatives.
#
# Arguments:   - polyveck v1: output vector of polynomials with
#                              coefficients a1
#              - polyveck v0: output vector of polynomials with
#                              coefficients a0
#              - polyveck v: input vector
##################################################
def polyveck_decompose(v1:polyveck, v0:polyveck, v:polyveck):
    for i in range(N):
        poly_decompose(v1.vec[i], v0.vec[i], v.vec[i])


#################################################
# Name:        polyveck_make_hint
#
# Description: Compute hint vector.
#
# Arguments:   - polyveck h: output vector
#              - polyveck v0: low part of input vector
#              - polyveck v1: high part of input vector
#
# Returns number of 1 bits.
##################################################
def polyveck_make_hint(h:polyveck, v0:polyveck, v1:polyveck):
    s = 0
    for i in range(K):
        s += poly_make_hint(h.vec[i], v0.vec[i], v1.vec[i])

    return s


#################################################
# Name:        polyveck_use_hint
#
# Description: Use hint vector to correct the high bits of input vector.
#
# Arguments:   - polyveck w: output vector of polynomials with
#                             corrected high bits
#              - polyveck *u: input vector
#              - polyveck *h: input hint vector
##################################################
def polyveck_use_hint(w:polyveck, u:polyveck, h:polyveck):
    for i in range(K):
        poly_use_hint(w.vec[i], u.vec[i], h.vec[i])


def polyveck_pack_w1(r:List[int], w1:polyveck):
    for i in range(K):
        polyw1_pack(r[i*POLYW1_PACKEDBYTES:], w1.vec[i])
