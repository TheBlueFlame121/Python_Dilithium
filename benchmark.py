from sign import *
from timeit import timeit

for mode in [2, 3, 5]:
    g.set_mode(mode)
    print(f"Dilithium{mode}")

    t = timeit("crypto_sign_keypair(pk, sk)", setup="pk, sk= [0]*g.CRYPTO_PUBLICKEYBYTES, [0]*g.CRYPTO_SECRETKEYBYTES", globals=globals(), number=500)
    print(f"key generation: {round(t, 3)}s")

    stp = """pk, sk= [0]*g.CRYPTO_PUBLICKEYBYTES, [0]*g.CRYPTO_SECRETKEYBYTES;crypto_sign_keypair(pk, sk);sig = [0]*g.CRYPTO_BYTES;m = list(urandom(32))"""
    t = timeit("crypto_sign_signature(sig, None, m, None, sk)", setup = stp, globals=globals(), number=500)
    print(f"signing 32 byte msg: {round(t, 3)}s")

    stp  = """pk, sk= [0]*g.CRYPTO_PUBLICKEYBYTES, [0]*g.CRYPTO_SECRETKEYBYTES;crypto_sign_keypair(pk, sk);sig = [0]*g.CRYPTO_BYTES;m = list(urandom(32));crypto_sign_signature(sig, None, m, None, sk)
    """
    t = timeit("crypto_sign_verify(sig, None, m, None, pk)", setup=stp, globals=globals(), number=500)
    print(f"verifying 32 byte msg: {round(t, 3)}s\n")
