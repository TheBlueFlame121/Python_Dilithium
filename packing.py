# Contains elements from packing.h and packing.c
from params import *
from polyvec import *
from poly import *


#################################################
# Name:        pack_pk
#
# Description: Bit-pack public key pk = (rho, t1).
#
# Arguments:   - List[int] pk: output byte array
#              - List[int] rho: byte array containing rho
#              - polyveck t1: vector t1
##################################################
def pack_pk(pk:List[int], rho:List[int], t1:polyveck):
    for i in range(g.SEEDBYTES):
        pk[i] = rho[i]
    
    temp = [0]*g.POLYT1_PACKEDBYTES
    for i in range(g.K):
        polyt1_pack(temp, t1.vec[i])
        for j in range(g.POLYT1_PACKEDBYTES):
            pk[g.SEEDBYTES + i*g.POLYT1_PACKEDBYTES + j] = temp[j]


#################################################
# Name:        unpack_pk
#
# Description: Unpack public key pk = (rho, t1).
#
# Arguments:   - List[int] rho: output byte array for rho
#              - polyveck t1: output vector t1
#              - List[int] pk: byte array containing bit-packed pk
##################################################
def unpack_pk(rho:List[int], t1:polyveck, pk:List[int]):
    for i in range(g.SEEDBYTES):
        rho[i] = pk[i]

    for i in range(g.K):
        polyt1_unpack(t1.vec[i], pk[g.SEEDBYTES + i*g.POLYT1_PACKEDBYTES: g.SEEDBYTES + (i+1)*g.POLYT1_PACKEDBYTES])


#################################################
# Name:        pack_sk
#
# Description: Bit-pack secret key sk = (rho, tr, key, t0, s1, s2).
#
# Arguments:   - List[int] sk: output byte array
#              - List[int] rho: byte array containing rho
#              - List[int] tr: byte array containing tr
#              - List[int] key: byte array containing key
#              - polyveck t0: vector t0
#              - polyvecl s1: vector s1
#              - polyveck s2: vector s2
##################################################
def pack_sk(sk:List[int], rho:List[int], tr:List[int], key:List[int], t0:polyveck, s1:polyvecl, s2:polyveck):
    start = 0
    for i in range(g.SEEDBYTES):
        sk[i] = rho[i]
    start += g.SEEDBYTES

    for i in range(g.SEEDBYTES):
        sk[start+i] = key[i]
    start += g.SEEDBYTES

    for i in range(g.SEEDBYTES):
        sk[start+i] = tr[i]
    start += g.SEEDBYTES

    for i in range(g.L):
        temp = [0]*g.POLYETA_PACKEDBYTES
        polyeta_pack(temp, s1.vec[i])
        for j in range(g.POLYETA_PACKEDBYTES):
            sk[start+i*g.POLYETA_PACKEDBYTES+j] = temp[j]
    start += g.L*g.POLYETA_PACKEDBYTES

    for i in range(g.K):
        temp = [0]*g.POLYETA_PACKEDBYTES
        polyeta_pack(temp, s2.vec[i])
        for j in range(g.POLYETA_PACKEDBYTES):
            sk[start+i*g.POLYETA_PACKEDBYTES+j] = temp[j]
    start += g.K*g.POLYETA_PACKEDBYTES

    for i in range(g.K):
        temp = [0]*g.POLYT0_PACKEDBYTES
        polyt0_pack(temp, t0.vec[i])
        for j in range(g.POLYT0_PACKEDBYTES):
            sk[start+i*g.POLYT0_PACKEDBYTES+j] = temp[j]


#################################################
# Name:        unpack_sk
#
# Description: Unpack secret key sk = (rho, tr, key, t0, s1, s2).
#
# Arguments:   - List[int] rho: output byte array for rho
#              - List[int] tr: output byte array for tr
#              - List[int] key: output byte array for key
#              - polyveck t0: output vector t0
#              - polyvecl s1: output vector s1
#              - polyveck s2: output vector s2
#              - List[int] sk: byte array containing bit-packed sk
##################################################
def unpack_sk(rho:List[int], tr:List[int], key:List[int], t0:polyveck, s1:polyvecl, s2:polyveck, sk:List[int]):
    start = 0
    for i in range(g.SEEDBYTES):
        rho[i] = sk[i]
    start += g.SEEDBYTES

    for i in range(g.SEEDBYTES):
        key[i] = sk[start+i]
    start += g.SEEDBYTES

    for i in range(g.SEEDBYTES):
        tr[i] = sk[start+i]
    start += g.SEEDBYTES

    for i in range(g.L):
        polyeta_unpack(s1.vec[i], sk[start+i*g.POLYETA_PACKEDBYTES:start+(i+1)*g.POLYETA_PACKEDBYTES])
    start += g.L*g.POLYETA_PACKEDBYTES

    for i in range(g.K):
        polyeta_unpack(s2.vec[i], sk[start+i*g.POLYETA_PACKEDBYTES:start+(i+1)*g.POLYETA_PACKEDBYTES])
    start += g.K*g.POLYETA_PACKEDBYTES

    for i in range(g.K):
        polyt0_unpack(t0.vec[i], sk[start+i*g.POLYT0_PACKEDBYTES:start+(i+1)*g.POLYT0_PACKEDBYTES])


#################################################
# Name:        pack_sig
#
# Description: Bit-pack signature sig = (c, z, h).
#
# Arguments:   - List[int] sig[]: output byte array
#              - List[int] c: challenge hash length SEEDBYTES
#              - polyvecl z: vector z
#              - polyveck h: hint vector h
##################################################
def pack_sig(sig:List[int], c:List[int], z:polyvecl, h:polyveck):
    start = 0
    for i in range(g.SEEDBYTES):
        sig[i] = c[i]
    start += g.SEEDBYTES

    for i in range(g.L):
        temp = [0]*g.POLYZ_PACKEDBYTES
        polyz_pack(temp, z.vec[i])
        for j in range(g.POLYZ_PACKEDBYTES):
            sig[start+i*g.POLYZ_PACKEDBYTES+j] = temp[j]
    start += g.L*g.POLYZ_PACKEDBYTES

    for i in range(g.OMEGA):
        sig[start+i] = 0

    k=0
    for i in range(g.K):
        for j in range(g.N):
            if h.vec[i].coeffs[j] != 0:
                sig[start+k] = j
                k += 1
    
        sig[start+g.OMEGA+i] = k


#################################################
# Name:        unpack_sig
#
# Description: Unpack signature sig = (c, z, h).
#
# Arguments:   - List[int] c: pointer to output challenge hash
#              - polyvecl z: pointer to output vector z
#              - polyveck h: pointer to output hint vector h
#              - List[int] sig: byte array containing
#                bit-packed signature
#
# Returns 1 in case of malformed signature; otherwise 0.
##################################################
def unpack_sig(c:List[int], z:polyvecl, h:polyveck, sig:List[int]) -> int:
    start = 0
    for i in range(g.SEEDBYTES):
        c[i] = sig[i]
    start += g.SEEDBYTES

    for i in range(g.L):
        polyz_unpack(z.vec[i], sig[start+i*g.POLYZ_PACKEDBYTES:start+(i+1)*g.POLYZ_PACKEDBYTES])
    start += g.L*g.POLYZ_PACKEDBYTES

    k = 0
    for i in range(g.K):
        for j in range(g.N):
            h.vec[i].coeffs[j] = 0

        if sig[start+g.OMEGA+i] < k or sig[start+g.OMEGA+i] > g.OMEGA:
            return 1

        for j in range(k, sig[start+g.OMEGA+i]):
            if j > k and sig[start+j] <= sig[start+j-1]:
                return 1
            h.vec[i].coeffs[sig[start+j]] = 1
        k = sig[start+g.OMEGA+i]

    for j in range(k, g.OMEGA):
        if sig[start+j]:
            return 1

    return 0
