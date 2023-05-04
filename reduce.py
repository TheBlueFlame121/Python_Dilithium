# Contains elements from reduce.c and reduce.h

from params import *

MONT = -4186625 # 2^32 % Q
QINV = 58728449 # q^(-1) mod 2^32


##################################################
# Name:        montgomery_reduce
#
# Description: For finite field element a with -2^{31}Q <= a <= Q*2^31,
#              compute r \equiv a*2^{-32} (mod Q) such that -Q < r < Q.
#
# Arguments:   - int64_t: finite field element a
#
# Returns r.
##################################################
def montgomery_reduce(a:int) -> int:
    t = (a*QINV)%(2**32)
    t = (a - t*Q) >> 32
    return t


##################################################
# Name:        reduce32
#
# Description: For finite field element a with a <= 2^{31} - 2^{22} - 1,
#              compute r \equiv a (mod Q) such that -6283009 <= r <= 6283007.
#
# Arguments:   - int32_t: finite field element a
#
# Returns r.
##################################################
def reduce32(a:int) -> int:
    t = (a + (1 << 22)) >> 23
    t = a - t*Q
    return t


##################################################
# Name:        caddq
#
# Description: Add Q if input coefficient is negative.
#
# Arguments:   - int32_t: finite field element a
#
# Returns r.
##################################################
def caddq(a:int) -> int:
    a += (a >> 31) & Q
    return a


##################################################
# Name:        freeze
#
# Description: For finite field element a, compute standard
#              representative r = a mod^+ Q.
#
# Arguments:   - int32_t: finite field element a
#
# Returns r.
##################################################
def freeze(a:int) -> int:
    a = reduce32(a)
    a = caddq(a)
    return a
