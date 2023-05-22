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
    for i in range(SEEDBYTES):
        pk[i] = rho[i]
    
    temp = [0]*POLYT1_PACKEDBYTES
    for i in range(K):
        polyt1_pack(temp, t1.vec[i])
        for j in range(POLYT1_PACKEDBYTES):
            pk[SEEDBYTES + i*POLYT1_PACKEDBYTES + j] = temp[j]


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
    for i in range(SEEDBYTES):
        rho[i] = pk[i]

    for i in range(K):
        polyt1_unpack(t1.vec[i], pk[SEEDBYTES + i*POLYT1_PACKEDBYTES: SEEDBYTES + (i+1)*POLYT1_PACKEDBYTES])


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
    for i in range(SEEDBYTES):
        sk[i] = rho[i]
    start += SEEDBYTES

    for i in range(SEEDBYTES):
        sk[start+i] = key[i]
    start += SEEDBYTES

    for i in range(SEEDBYTES):
        sk[start+i] = tr[i]
    start += SEEDBYTES

    for i in range(L):
        temp = [0]*POLYETA_PACKEDBYTES
        polyeta_pack(temp, s1.vec[i])
        for j in range(POLYETA_PACKEDBYTES):
            sk[start+i*POLYETA_PACKEDBYTES+j] = temp[j]
    start += L*POLYETA_PACKEDBYTES

    for i in range(K):
        temp = [0]*POLYETA_PACKEDBYTES
        polyeta_pack(temp, s2.vec[i])
        for j in range(POLYETA_PACKEDBYTES):
            sk[start+i*POLYETA_PACKEDBYTES+j] = temp[j]
    start += K*POLYETA_PACKEDBYTES

    for i in range(K):
        temp = [0]*POLYT0_PACKEDBYTES
        polyt0_pack(temp, t0.vec[i])
        for j in range(POLYT0_PACKEDBYTES):
            sk[start+i*POLYT0_PACKEDBYTES+j] = temp[j]


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
    for i in range(SEEDBYTES):
        rho[i] = sk[i]
    start += SEEDBYTES

    for i in range(SEEDBYTES):
        key[i] = sk[start+i]
    start += SEEDBYTES

    for i in range(SEEDBYTES):
        tr[i] = sk[start+i]
    start += SEEDBYTES

    for i in range(L):
        polyeta_unpack(s1.vec[i], sk[start+i*POLYETA_PACKEDBYTES:start+(i+1)*POLYETA_PACKEDBYTES])
    start += L*POLYETA_PACKEDBYTES

    for i in range(K):
        polyeta_unpack(s2.vec[i], sk[start+i*POLYETA_PACKEDBYTES:start+(i+1)*POLYETA_PACKEDBYTES])
    start += K*POLYETA_PACKEDBYTES

    for i in range(K):
        polyt0_unpack(t0.vec[i], sk[start+i*POLYT0_PACKEDBYTES:start+(i+1)*POLYT0_PACKEDBYTES])


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
    for i in range(SEEDBYTES):
        sig[i] = c[i]
    start += SEEDBYTES

    for i in range(L):
        temp = [0]*POLYZ_PACKEDBYTES
        polyz_pack(temp, z.vec[i])
        for j in range(POLYZ_PACKEDBYTES):
            sig[start+i*POLYZ_PACKEDBYTES+j] = temp[j]
    start += L*POLYZ_PACKEDBYTES

    for i in range(OMEGA):
        sig[start+i] = 0

    k=0
    for i in range(K):
        for j in range(N):
            if h.vec[i].coeffs[j] != 0:
                sig[start+k] = j
                k += 1
    
        sig[start+OMEGA+i] = k


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
    for i in range(SEEDBYTES):
        c[i] = sig[i]
    start += SEEDBYTES

    for i in range(L):
        polyz_unpack(z.vec[i], sig[start+i*POLYZ_PACKEDBYTES:start+(i+1)*POLYZ_PACKEDBYTES])
    start += L*POLYZ_PACKEDBYTES

    k = 0
    for i in range(K):
        for j in range(N):
            h.vec[i].coeffs[j] = 0

        if sig[start+OMEGA+i] < k or sig[start+OMEGA+i] > OMEGA:
            return 1

        for j in range(k, sig[start+OMEGA+i]):
            if j > k and sig[start+j] <= sig[start+j-1]:
                return 1
            h.vec[i].coeffs[sig[start+j]] = 1
        k = sig[start+OMEGA+i]

    for j in range(k, OMEGA):
        if sig[start+j]:
            return 1

    return 0
