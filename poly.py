# Contains elements from poly.h and poly.c
from dataclasses import dataclass
from params import *
from ntt import *
from reduce import *
from rounding import *
from symmetric import *


class poly:
    coeffs: List[int]

    def __init__(self, inp: list[int] = None):
        if inp is None:
            inp = [0 for _ in range(N)]
        if len(inp) > N:
            raise ValueError("Polynomial can't have more than N coeffs")
        if len(inp) < N:
            inp = inp + [0 for _ in range(N-len(inp))].copy()
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


#################################################
# Name:        poly_uniform
#
# Description: Sample polynomial with uniformly random coefficients
#              in [0,Q-1] by performing rejection sampling on the
#              output stream of SHAKE256(seed|nonce) or AES256CTR(seed,nonce).
#
# Arguments:   - poly a: output polynomial
#              - List[int] seed: byte array with seed of length SEEDBYTES
#              - int nonce: 2-byte nonce
##################################################
POLY_UNIFORM_NBLOCKS = ((768 + STREAM128_BLOCKBYTES - 1)//STREAM128_BLOCKBYTES)
def poly_uniform(a: poly, seed: bytes, nonce: int):
    buflen = POLY_UNIFORM_NBLOCKS*STREAM128_BLOCKBYTES
    buf = [0]*(POLY_UNIFORM_NBLOCKS*STREAM128_BLOCKBYTES + 2)
    state = stream128_state()

    stream128_init(state, seed, nonce)
    buf = stream128_squeezeblocks(POLY_UNIFORM_NBLOCKS, state)

    ctr = rej_uniform(a.coeffs, N, buf, buflen)

    while ctr < N:
        off = buflen % 3
        for i in range(off):
            buf[i] = buf[buflen - off + i]

        buf = buf[:off] + stream128_squeezeblocks(1, state)
        buflen = STREAM128_BLOCKBYTES + off;
        temp = [0]*(N-ctr)
        ctr_temp = rej_uniform(temp, N-ctr, buf, buflen)
        for i in range(N-ctr):
            a.coeffs[ctr+i] = temp[i]
        ctr+=ctr_temp



#################################################
# Name:        rej_eta
#
# Description: Sample uniformly random coefficients in [-ETA, ETA] by
#              performing rejection sampling on array of random bytes.
#
# Arguments:   - List[int] a: output array (allocated)
#              - int l: number of coefficients to be sampled
#              - List[int] buf: array of random bytes
#              - int buflen: length of array of random bytes
#
# Returns number of sampled coefficients. Can be smaller than len if not enough
# random bytes were given.
##################################################
def rej_eta(a:List[int], l:int, buf:List[int], buflen:int) -> int:
    ctr = pos = 0
    while ctr<l and pos<buflen:
        t0 = buf[pos] & 0x0F
        t1 = buf[pos] >> 4
        pos+=1

        if ETA == 2:
            if t0<15:
                t0 = t0 - (205*t0 >> 10)*5
                a[ctr] = 2-t0
                ctr+=1
            if t1<15 and ctr<l:
                t1 = t1 - (205*t1 >> 10)*5
                a[ctr] = 2- t1
                ctr += 1
        elif ETA == 4:
            if t0<9:
                a[ctr] = 4 - t0
                ctr += 1
            if t1<9 and ctr<l:
                a[ctr] = 4-t1
                ctr += 1
    return ctr


#################################################
# Name:        poly_uniform_eta
#
# Description: Sample polynomial with uniformly random coefficients
#              in [-ETA,ETA] by performing rejection sampling on the
#              output stream from SHAKE256(seed|nonce) or AES256CTR(seed,nonce).
#
# Arguments:   - poly a: output polynomial
#              - List[int] seed: byte array with seed of length CRHBYTES
#              - int nonce: 2-byte nonce
##################################################
if ETA == 2:
    POLY_UNIFORM_ETA_NBLOCKS = ((136 + STREAM256_BLOCKBYTES - 1)//STREAM256_BLOCKBYTES)
elif ETA == 4:
    POLY_UNIFORM_ETA_NBLOCKS = ((227 + STREAM256_BLOCKBYTES - 1)//STREAM256_BLOCKBYTES)
def poly_uniform_eta(a:poly, seed:List[int], nonce:int):
    buflen = POLY_UNIFORM_ETA_NBLOCKS*STREAM256_BLOCKBYTES
    buf = [0]*(POLY_UNIFORM_ETA_NBLOCKS*STREAM256_BLOCKBYTES)
    state = stream256_state()

    stream256_init(state, bytes(seed), nonce)
    buf = stream256_squeezeblocks(POLY_UNIFORM_ETA_NBLOCKS, state)

    ctr = rej_eta(a.coeffs, N, buf, buflen)

    while ctr<N:
        buf = stream256_squeezeblocks(1, state) 
        temp = [0]*(N-ctr)
        ctr_temp = rej_eta(temp, N-ctr, buf, STREAM256_BLOCKBYTES)
        for i in range(N-ctr):
            a.coeffs[ctr+i] = temp[i]
        ctr+=ctr_temp


#################################################
# Name:        poly_uniform_gamma1m1
#
# Description: Sample polynomial with uniformly random coefficients
#              in [-(GAMMA1 - 1), GAMMA1] by unpacking output stream
#              of SHAKE256(seed|nonce) or AES256CTR(seed,nonce).
#
# Arguments:   - poly a: output polynomial
#              - List[int] seed: byte array with seed of length CRHBYTES
#              - int nonce: 16-bit nonce
##################################################
POLY_UNIFORM_GAMMA1_NBLOCKS = ((POLYZ_PACKEDBYTES + STREAM256_BLOCKBYTES - 1)//STREAM256_BLOCKBYTES)
def poly_uniform_gamma1(a:poly, seed:List[int], nonce:int):
    buf = [0]*(POLY_UNIFORM_GAMMA1_NBLOCKS*STREAM256_BLOCKBYTES)
    state = stream256_state()

    stream256_init(state, bytes(seed), nonce)
    buf = stream256_squeezeblocks(POLY_UNIFORM_GAMMA1_NBLOCKS, state)
    polyz_unpack(a, buf)


#################################################
# Name:        challenge
#
# Description: Implementation of H. Samples polynomial with TAU nonzero
#              coefficients in {-1,1} using the output stream of
#              SHAKE256(seed).
#
# Arguments:   - poly c: output polynomial
#              - List[int] seed: byte array containing seed of length SEEDBYTES
##################################################
def poly_challenge(c:poly, seed:List[int]):
    buf = [0]*SHAKE256_RATE
    state = SHAKE256.new()

    state.update(bytes(seed))
    buf = shake256_squeezeblocks(1, state)

    signs = 0
    for i in range(8):
        signs |= buf[i] << 8*i
    pos = 8

    for i in range(N):
        c.coeffs[i] = 0
    for i in range(N-TAU, N):
        temp = 1; b = 0 # To simulate a do while loop
        while temp == 1 or b > i:
            temp = 0
            if pos >= SHAKE256_RATE:
                buf = shake256_squeezeblocks(1, state)
                pos = 0

            b = buf[pos]
            pos += 1

        c.coeffs[i] = c.coeffs[b]
        c.coeffs[b] = 1 - 2*(signs & 1)
        signs >>= 1


#################################################
# Name:        polyeta_pack
#
# Description: Bit-pack polynomial with coefficients in [-ETA,ETA].
#
# Arguments:   - List[int] r: output byte array with at least
#                            POLYETA_PACKEDBYTES bytes
#              - poly a: input polynomial
##################################################
def polyeta_pack(r:List[int], a:poly):
    t = [0]*8

    if ETA == 2:
        for i in range(N//8):
            t[0] = ETA - a.coeffs[8*i+0]
            t[1] = ETA - a.coeffs[8*i+1]
            t[2] = ETA - a.coeffs[8*i+2]
            t[3] = ETA - a.coeffs[8*i+3]
            t[4] = ETA - a.coeffs[8*i+4]
            t[5] = ETA - a.coeffs[8*i+5]
            t[6] = ETA - a.coeffs[8*i+6]
            t[7] = ETA - a.coeffs[8*i+7]

            r[3*i+0]  = (t[0] >> 0) | (t[1] << 3) | (t[2] << 6)
            r[3*i+1]  = (t[2] >> 2) | (t[3] << 1) | (t[4] << 4) | (t[5] << 7)
            r[3*i+2]  = (t[5] >> 1) | (t[6] << 2) | (t[7] << 5)
    elif ETA == 4:
        for i in range(N//2):
            t[0] = ETA - a.coeffs[2*i+0]
            t[1] = ETA - a.coeffs[2*i+1]
            r[i] = t[0] | (t[1] << 4)


#################################################
# Name:        polyeta_unpack
#
# Description: Unpack polynomial with coefficients in [-ETA,ETA].
#
# Arguments:   - poly r: output polynomial
#              - List[int] a: byte array with bit-packed polynomial
##################################################
def polyeta_unpack(r:poly, a:List[int]):
    if ETA == 2:
        for i in range(N//8):
            r.coeffs[8*i+0] =  (a[3*i+0] >> 0) & 7
            r.coeffs[8*i+1] =  (a[3*i+0] >> 3) & 7
            r.coeffs[8*i+2] = ((a[3*i+0] >> 6) | (a[3*i+1] << 2)) & 7
            r.coeffs[8*i+3] =  (a[3*i+1] >> 1) & 7
            r.coeffs[8*i+4] =  (a[3*i+1] >> 4) & 7
            r.coeffs[8*i+5] = ((a[3*i+1] >> 7) | (a[3*i+2] << 1)) & 7
            r.coeffs[8*i+6] =  (a[3*i+2] >> 2) & 7
            r.coeffs[8*i+7] =  (a[3*i+2] >> 5) & 7

            r.coeffs[8*i+0] = ETA - r.coeffs[8*i+0]
            r.coeffs[8*i+1] = ETA - r.coeffs[8*i+1]
            r.coeffs[8*i+2] = ETA - r.coeffs[8*i+2]
            r.coeffs[8*i+3] = ETA - r.coeffs[8*i+3]
            r.coeffs[8*i+4] = ETA - r.coeffs[8*i+4]
            r.coeffs[8*i+5] = ETA - r.coeffs[8*i+5]
            r.coeffs[8*i+6] = ETA - r.coeffs[8*i+6]
            r.coeffs[8*i+7] = ETA - r.coeffs[8*i+7]
    elif ETA == 4:
        for i in range(N//2):
            r.coeffs[2*i+0] = a[i] & 0x0F
            r.coeffs[2*i+1] = a[i] >> 4
            r.coeffs[2*i+0] = ETA - r.coeffs[2*i+0]
            r.coeffs[2*i+1] = ETA - r.coeffs[2*i+1]


#################################################
# Name:        polyt1_pack
#
# Description: Bit-pack polynomial t1 with coefficients fitting in 10 bits.
#              Input coefficients are assumed to be standard representatives.
#
# Arguments:   - List[int] r: output byte array with at least
#                            POLYT1_PACKEDBYTES bytes
#              - poly a: input polynomial
##################################################
def polyt1_pack(r:List[int], a:poly):
    for i in range(N//4):
        r[5*i+0] = (a.coeffs[4*i+0] >> 0)
        r[5*i+1] = (a.coeffs[4*i+0] >> 8) | (a.coeffs[4*i+1] << 2)
        r[5*i+2] = (a.coeffs[4*i+1] >> 6) | (a.coeffs[4*i+2] << 4)
        r[5*i+3] = (a.coeffs[4*i+2] >> 4) | (a.coeffs[4*i+3] << 6)
        r[5*i+4] = (a.coeffs[4*i+3] >> 2)


#################################################
# Name:        polyt1_unpack
#
# Description: Unpack polynomial t1 with 10-bit coefficients.
#              Output coefficients are standard representatives.
#
# Arguments:   - poly r: output polynomial
#              - List[int] a: byte array with bit-packed polynomial
##################################################
def polyt1_unpack(r:poly, a:List[int]):
    for i in range(N//4):
        r.coeffs[4*i+0] = ((a[5*i+0] >> 0) | (a[5*i+1] << 8)) & 0x3FF
        r.coeffs[4*i+1] = ((a[5*i+1] >> 2) | (a[5*i+2] << 6)) & 0x3FF
        r.coeffs[4*i+2] = ((a[5*i+2] >> 4) | (a[5*i+3] << 4)) & 0x3FF
        r.coeffs[4*i+3] = ((a[5*i+3] >> 6) | (a[5*i+4] << 2)) & 0x3FF


#################################################
# Name:        polyt0_pack
#
# Description: Bit-pack polynomial t0 with coefficients in [-2^{D-1}, 2^{D-1}].
#
# Arguments:   - List[int] r: output byte array with at least
#                            POLYT0_PACKEDBYTES bytes
#              - poly a: input polynomial
##################################################
def polyt0_pack(r:List[int], a:poly):
    t = [0]*8
    for i in range(N//8):
        t[0] = (1 << (D-1)) - a.coeffs[8*i+0]
        t[1] = (1 << (D-1)) - a.coeffs[8*i+1]
        t[2] = (1 << (D-1)) - a.coeffs[8*i+2]
        t[3] = (1 << (D-1)) - a.coeffs[8*i+3]
        t[4] = (1 << (D-1)) - a.coeffs[8*i+4]
        t[5] = (1 << (D-1)) - a.coeffs[8*i+5]
        t[6] = (1 << (D-1)) - a.coeffs[8*i+6]
        t[7] = (1 << (D-1)) - a.coeffs[8*i+7]

        r[13*i+ 0]  =  t[0]
        r[13*i+ 1]  =  t[0] >>  8
        r[13*i+ 1] |=  t[1] <<  5
        r[13*i+ 2]  =  t[1] >>  3
        r[13*i+ 3]  =  t[1] >> 11
        r[13*i+ 3] |=  t[2] <<  2
        r[13*i+ 4]  =  t[2] >>  6
        r[13*i+ 4] |=  t[3] <<  7
        r[13*i+ 5]  =  t[3] >>  1
        r[13*i+ 6]  =  t[3] >>  9
        r[13*i+ 6] |=  t[4] <<  4
        r[13*i+ 7]  =  t[4] >>  4
        r[13*i+ 8]  =  t[4] >> 12
        r[13*i+ 8] |=  t[5] <<  1
        r[13*i+ 9]  =  t[5] >>  7
        r[13*i+ 9] |=  t[6] <<  6
        r[13*i+10]  =  t[6] >>  2
        r[13*i+11]  =  t[6] >> 10
        r[13*i+11] |=  t[7] <<  3
        r[13*i+12]  =  t[7] >>  5


#################################################
# Name:        polyt0_unpack
#
# Description: Unpack polynomial t0 with coefficients in ]-2^{D-1}, 2^{D-1}].
#
# Arguments:   - poly r: output polynomial
#              - List[int] a: byte array with bit-packed polynomial
##################################################
def polyt0_unpack(r:poly, a:List[int]):
    for i in range(N//8):
        r.coeffs[8*i+0]  = a[13*i+0]
        r.coeffs[8*i+0] |= a[13*i+1] << 8
        r.coeffs[8*i+0] &= 0x1FFF

        r.coeffs[8*i+1]  = a[13*i+1] >> 5
        r.coeffs[8*i+1] |= a[13*i+2] << 3
        r.coeffs[8*i+1] |= a[13*i+3] << 11
        r.coeffs[8*i+1] &= 0x1FFF

        r.coeffs[8*i+2]  = a[13*i+3] >> 2
        r.coeffs[8*i+2] |= a[13*i+4] << 6
        r.coeffs[8*i+2] &= 0x1FFF

        r.coeffs[8*i+3]  = a[13*i+4] >> 7
        r.coeffs[8*i+3] |= a[13*i+5] << 1
        r.coeffs[8*i+3] |= a[13*i+6] << 9
        r.coeffs[8*i+3] &= 0x1FFF

        r.coeffs[8*i+4]  = a[13*i+6] >> 4
        r.coeffs[8*i+4] |= a[13*i+7] << 4
        r.coeffs[8*i+4] |= a[13*i+8] << 12
        r.coeffs[8*i+4] &= 0x1FFF

        r.coeffs[8*i+5]  = a[13*i+8] >> 1
        r.coeffs[8*i+5] |= a[13*i+9] << 7
        r.coeffs[8*i+5] &= 0x1FFF

        r.coeffs[8*i+6]  = a[13*i+9] >> 6
        r.coeffs[8*i+6] |= a[13*i+10] << 2
        r.coeffs[8*i+6] |= a[13*i+11] << 10
        r.coeffs[8*i+6] &= 0x1FFF

        r.coeffs[8*i+7]  = a[13*i+11] >> 3
        r.coeffs[8*i+7] |= a[13*i+12] << 5
        r.coeffs[8*i+7] &= 0x1FFF

        r.coeffs[8*i+0] = (1 << (D-1)) - r.coeffs[8*i+0]
        r.coeffs[8*i+1] = (1 << (D-1)) - r.coeffs[8*i+1]
        r.coeffs[8*i+2] = (1 << (D-1)) - r.coeffs[8*i+2]
        r.coeffs[8*i+3] = (1 << (D-1)) - r.coeffs[8*i+3]
        r.coeffs[8*i+4] = (1 << (D-1)) - r.coeffs[8*i+4]
        r.coeffs[8*i+5] = (1 << (D-1)) - r.coeffs[8*i+5]
        r.coeffs[8*i+6] = (1 << (D-1)) - r.coeffs[8*i+6]
        r.coeffs[8*i+7] = (1 << (D-1)) - r.coeffs[8*i+7]


#################################################
# Name:        polyz_pack
#
# Description: Bit-pack polynomial with coefficients
#              in [-(GAMMA1 - 1), GAMMA1].
#
# Arguments:   - List[int] r: output byte array with at least
#                            POLYZ_PACKEDBYTES bytes
#              - poly a: input polynomial
##################################################
def polyz_pack(r:List[int], a:poly):
    t = [0]*4

    if GAMMA1 == (1 << 17):
        for i in range(N//4):
            t[0] = GAMMA1 - a.coeffs[4*i+0]
            t[1] = GAMMA1 - a.coeffs[4*i+1]
            t[2] = GAMMA1 - a.coeffs[4*i+2]
            t[3] = GAMMA1 - a.coeffs[4*i+3]

            r[9*i+0]  = t[0]
            r[9*i+1]  = t[0] >> 8
            r[9*i+2]  = t[0] >> 16
            r[9*i+2] |= t[1] << 2
            r[9*i+3]  = t[1] >> 6
            r[9*i+4]  = t[1] >> 14
            r[9*i+4] |= t[2] << 4
            r[9*i+5]  = t[2] >> 4
            r[9*i+6]  = t[2] >> 12
            r[9*i+6] |= t[3] << 6
            r[9*i+7]  = t[3] >> 2
            r[9*i+8]  = t[3] >> 10
    elif GAMMA1 == (1 << 19):
        for i in range(N//2):
            t[0] = GAMMA1 - a.coeffs[2*i+0]
            t[1] = GAMMA1 - a.coeffs[2*i+1]

            r[5*i+0]  = t[0]
            r[5*i+1]  = t[0] >> 8
            r[5*i+2]  = t[0] >> 16
            r[5*i+2] |= t[1] << 4
            r[5*i+3]  = t[1] >> 4
            r[5*i+4]  = t[1] >> 12


#################################################
# Name:        polyz_unpack
#
# Description: Unpack polynomial z with coefficients
#              in [-(GAMMA1 - 1), GAMMA1].
#
# Arguments:   - poly r: output polynomial
#              - List[int] a: byte array with bit-packed polynomial
##################################################
def polyz_unpack(r:poly, a:List[int]):
    if GAMMA1 == (1 << 17):
        for i in range(N//4):
            r.coeffs[4*i+0]  = a[9*i+0]
            r.coeffs[4*i+0] |= a[9*i+1] << 8
            r.coeffs[4*i+0] |= a[9*i+2] << 16
            r.coeffs[4*i+0] &= 0x3FFFF

            r.coeffs[4*i+1]  = a[9*i+2] >> 2
            r.coeffs[4*i+1] |= a[9*i+3] << 6
            r.coeffs[4*i+1] |= a[9*i+4] << 14
            r.coeffs[4*i+1] &= 0x3FFFF

            r.coeffs[4*i+2]  = a[9*i+4] >> 4
            r.coeffs[4*i+2] |= a[9*i+5] << 4
            r.coeffs[4*i+2] |= a[9*i+6] << 12
            r.coeffs[4*i+2] &= 0x3FFFF

            r.coeffs[4*i+3]  = a[9*i+6] >> 6
            r.coeffs[4*i+3] |= a[9*i+7] << 2
            r.coeffs[4*i+3] |= a[9*i+8] << 10
            r.coeffs[4*i+3] &= 0x3FFFF

            r.coeffs[4*i+0] = GAMMA1 - r.coeffs[4*i+0]
            r.coeffs[4*i+1] = GAMMA1 - r.coeffs[4*i+1]
            r.coeffs[4*i+2] = GAMMA1 - r.coeffs[4*i+2]
            r.coeffs[4*i+3] = GAMMA1 - r.coeffs[4*i+3]
    elif GAMMA1 == (1 << 19):
        for i in range(N//2):
            r.coeffs[2*i+0]  = a[5*i+0]
            r.coeffs[2*i+0] |= a[5*i+1] << 8
            r.coeffs[2*i+0] |= a[5*i+2] << 16
            r.coeffs[2*i+0] &= 0xFFFFF

            r.coeffs[2*i+1]  = a[5*i+2] >> 4
            r.coeffs[2*i+1] |= a[5*i+3] << 4
            r.coeffs[2*i+1] |= a[5*i+4] << 12
            r.coeffs[2*i+0] &= 0xFFFFF

            r.coeffs[2*i+0] = GAMMA1 - r.coeffs[2*i+0]
            r.coeffs[2*i+1] = GAMMA1 - r.coeffs[2*i+1]


#################################################
# Name:        polyw1_pack
#
# Description: Bit-pack polynomial w1 with coefficients in [0,15] or [0,43].
#              Input coefficients are assumed to be standard representatives.
#
# Arguments:   - List[int] r: output byte array with at least
#                            POLYW1_PACKEDBYTES bytes
#              - poly a: input polynomial
##################################################
def polyw1_pack(r:List[int], a:poly):
    if GAMMA2 == (Q-1)/88:
        for i in range(N//4):
            r[3*i+0]  = a.coeffs[4*i+0];
            r[3*i+0] |= a.coeffs[4*i+1] << 6
            r[3*i+1]  = a.coeffs[4*i+1] >> 2
            r[3*i+1] |= a.coeffs[4*i+2] << 4
            r[3*i+2]  = a.coeffs[4*i+2] >> 4
            r[3*i+2] |= a.coeffs[4*i+3] << 2
    elif GAMMA2 == (Q-1)/32:
        for i in range(N//2):
            r[i] = a.coeffs[2*i+0] | (a.coeffs[2*i+1] << 4)
