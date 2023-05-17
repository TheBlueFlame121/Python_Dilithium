from dataclasses import dataclass
from params import *
from ntt import *
from reduce import *
from rounding import *
from symmetric import *


@dataclass
class poly:
    coeffs: List[int]

    def __init__(self, inp: List[int] = [0]*N):
        if len(inp) > N:
            raise ValueError("Polynomial can't have more than 10 coeffs")
        if len(inp) < N:
            inp = [0]*(N-len(inp)) + inp
        self.coeffs = inp


#################################################
# Name:        poly_reduce
#
# Description: Inplace reduction of all coefficients of polynomial to
#              representative in [-6283009,6283007].
#
# Arguments:   - poly a: input/output polynomial
##################################################
def poly_reduce(a: poly):
    for i in range(N):
        a.coeffs[i] = reduce32(a.coeffs[i])


#################################################
# Name:        poly_caddq
#
# Description: For all coefficients of in/out polynomial add Q if
#              coefficient is negative.
#
# Arguments:   - poly a: input/output polynomial
##################################################
def poly_caddq(a: poly):
    for i in range(N):
        a.coeffs[i] = caddq(a.coeffs[i])


#################################################
# Name:        poly_add
#
# Description: Add polynomials. No modular reduction is performed.
#
# Arguments:   - poly c: output polynomial
#              - poly a: first summand
#              - poly b: second summand
##################################################
def poly_add(c: poly, a: poly, b: poly):
    for i in range(N):
        c.coeffs[i] = a.coeffs[i] + b.coeffs[i]


#################################################
# Name:        poly_sub
#
# Description: Subtract polynomials. No modular reduction is
#              performed.
#
# Arguments:   - poly c: output polynomial
#              - poly a: first input polynomial
#              - poly b: second input polynomial to be
#                           subtraced from first input polynomial
##################################################
def poly_sub(c: poly, a: poly, b: poly):
    for i in range(N):
        c.coeffs[i] = a.coeffs[i] - b.coeffs[i]


#################################################
# Name:        poly_shiftl
#
# Description: Multiply polynomial by 2^D without modular reduction. Assumes
#              input coefficients to be less than 2^{31-D} in absolute value.
#
# Arguments:   - poly a: input/output polynomial
##################################################
def poly_shiftl(a: poly):
    for i in range(N):
        a.coeffs[i] <<= D


#################################################
# Name:        poly_ntt
#
# Description: Inplace forward NTT. Coefficients can grow by
#              8*Q in absolute value.
#
# Arguments:   - poly a: input/output polynomial
##################################################
def poly_ntt(a: poly):
    ntt(a.coeffs)


#################################################
# Name:        poly_invntt_tomont
#
# Description: Inplace inverse NTT and multiplication by 2^{32}.
#              Input coefficients need to be less than Q in absolute
#              value and output coefficients are again bounded by Q.
#
# Arguments:   - poly a: input/output polynomial
##################################################
def poly_invntt_tomont(a: poly):
    invntt_tomont(a.coeffs)


#################################################
# Name:        poly_pointwise_montgomery
#
# Description: Pointwise multiplication of polynomials in NTT domain
#              representation and multiplication of resulting polynomial
#              by 2^{-32}.
#
# Arguments:   - poly c: output polynomial
#              - poly a: first input polynomial
#              - poly b: second input polynomial
##################################################
def poly_pointwise_montgomery(c: poly, a: poly, b: poly):
    for i in range(N):
        c.coeffs[i] = montgomery_reduce(a.coeffs[i] * b.coeffs[i])


#################################################
# Name:        poly_power2round
#
# Description: For all coefficients c of the input polynomial,
#              compute c0, c1 such that c mod Q = c1*2^D + c0
#              with -2^{D-1} < c0 <= 2^{D-1}. Assumes coefficients to be
#              standard representatives.
#
# Arguments:   - poly a1: output polynomial with coefficients c1
#              - poly a0: output polynomial with coefficients c0
#              - poly a: input polynomial
##################################################
def poly_power2round(a1: poly, a0: poly, a: poly):
    for i in range(N):
        a0.coeffs[i], a1.coeffs[i] = power2round(a.coeffs[i])


#################################################
# Name:        poly_decompose
#
# Description: For all coefficients c of the input polynomial,
#              compute high and low bits c0, c1 such c mod Q = c1*ALPHA + c0
#              with -ALPHA/2 < c0 <= ALPHA/2 except c1 = (Q-1)/ALPHA where we
#              set c1 = 0 and -ALPHA/2 <= c0 = c mod Q - Q < 0.
#              Assumes coefficients to be standard representatives.
#
# Arguments:   - poly a1: output polynomial with coefficients c1
#              - poly a0: output polynomial with coefficients c0
#              - poly a: input polynomial
##################################################
def poly_decompose(a1: poly, a0: poly, a: poly):
    for i in range(N):
        a0.coeffs[i], a1.coeffs[i] = decompose(a.coeffs[i])


#################################################
# Name:        poly_make_hint
#
# Description: Compute hint polynomial. The coefficients of which indicate
#              whether the low bits of the corresponding coefficient of
#              the input polynomial overflow into the high bits.
#
# Arguments:   - poly h: output hint polynomial
#              - poly a0: low part of input polynomial
#              - poly a1: high part of input polynomial
#
# Returns number of 1 bits.
##################################################
def poly_make_hint(h: poly, a0: poly, a1: poly):
    s = 0
    for i in range(N):
        h.coeffs[i] = int(make_hint(a0.coeffs[i], a1.coeffs[i]))
        s += h.coeffs[i]
    return s


#################################################
# Name:        poly_use_hint
#
# Description: Use hint polynomial to correct the high bits of a polynomial.
#
# Arguments:   - poly b: output polynomial with corrected high bits
#              - poly a: input polynomial
#              - poly h: input hint polynomial
##################################################
def poly_use_hint(b: poly, a: poly, h:poly):
    for i in range(N):
        b.coeffs[i] = use_hint(a.coeffs[i], h.coeffs[i])


#################################################
# Name:        poly_chknorm
#
# Description: Check infinity norm of polynomial against given bound.
#              Assumes input coefficients were reduced by reduce32().
#
# Arguments:   - poly a: polynomial
#              - int B: norm bound
#
# Returns 0 if norm is strictly smaller than B <= (Q-1)/8 and 1 otherwise.
##################################################
def poly_chknorm(a: poly, B: int) -> int:
    if B > (Q-1)//8:
        return 1

    for i in range(N):
        t = a.coeffs[i] >> 31
        t = a.coeffs[i] - (t & 2*a.coeffs[i])

    if t >= B:
        return 1

    return 0


#################################################
# Name:        rej_uniform
#
# Description: Sample uniformly random coefficients in [0, Q-1] by
#              performing rejection sampling on array of random bytes.
#
# Arguments:   - List[int] a: output array (allocated)
#              - int l: number of coefficients to be sampled
#              - int buf: array of random bytes
#              - int buflen: length of array of random bytes
#
# Returns number of sampled coefficients. Can be smaller than len if not enough
# random bytes were given.
##################################################
def rej_uniform(a:List[int], l:int, buf: List[int], buflen:int) -> int:
    ctr = pos = 0
    while ctr<l and pos+3 <= buflen:
        t = buf[pos]
        pos+=1
        t |= buf[pos] << 8
        pos+=1
        t |= buf[pos] << 16
        pos+=1
        t &= 0x7FFFFF

        if t< Q:
            a[ctr] = t
            ctr += 1
    return ctr
