from sign import *


def test_dilithium2():
    f = open("KATs/KAT_Dilithium2.rsp", "r")
    g.set_mode(2)
    f.readline()
    f.readline()
    for i in range(100):
        count = int(f.readline().split()[-1])
        seed = bytes.fromhex(f.readline().split()[-1])
        mlen = int(f.readline().split()[-1])
        msg = bytes.fromhex(f.readline().split()[-1])
        pk_kat = bytes.fromhex(f.readline().split()[-1])
        sk_kat = bytes.fromhex(f.readline().split()[-1])
        siglen = int(f.readline().split()[-1])
        sig_kat = bytes.fromhex(f.readline().split()[-1])

        pk = [0]*g.CRYPTO_PUBLICKEYBYTES
        sk = [0]*g.CRYPTO_SECRETKEYBYTES
        crypto_sign_keypair(pk, sk, seed)
        assert bytes(pk) == pk_kat
        assert bytes(sk) == sk_kat


        sig = [0]*g.CRYPTO_BYTES
        crypto_sign_signature(sig, siglen, list(msg), mlen, sk)

        assert bytes(sig) == sig_kat

        assert crypto_sign_verify(sig, siglen, list(msg), mlen, pk) == 0

        f.readline()
        # print(f"Case {count} passed")
    print("Dilithium 2 passes all KATs")


def test_dilithium3():
    g.set_mode(3)
    f = open("KATs/KAT_Dilithium3.rsp", "r")
    f.readline()
    f.readline()
    for i in range(100):
        count = int(f.readline().split()[-1])
        seed = bytes.fromhex(f.readline().split()[-1])
        mlen = int(f.readline().split()[-1])
        msg = bytes.fromhex(f.readline().split()[-1])
        pk_kat = bytes.fromhex(f.readline().split()[-1])
        sk_kat = bytes.fromhex(f.readline().split()[-1])
        siglen = int(f.readline().split()[-1])
        sig_kat = bytes.fromhex(f.readline().split()[-1])

        pk = [0]*g.CRYPTO_PUBLICKEYBYTES
        sk = [0]*g.CRYPTO_SECRETKEYBYTES
        crypto_sign_keypair(pk, sk, seed)
        assert bytes(pk) == pk_kat
        assert bytes(sk) == sk_kat


        sig = [0]*g.CRYPTO_BYTES
        crypto_sign_signature(sig, siglen, list(msg), mlen, sk)

        assert bytes(sig) == sig_kat

        assert crypto_sign_verify(sig, siglen, list(msg), mlen, pk) == 0

        f.readline()
        # print(f"Case {count} passed")
    print("Dilithium 3 passes all KATs")


def test_dilithium5():
    g.set_mode(5)
    f = open("KATs/KAT_Dilithium5.rsp", "r")
    f.readline()
    f.readline()
    for i in range(100):
        count = int(f.readline().split()[-1])
        seed = bytes.fromhex(f.readline().split()[-1])
        mlen = int(f.readline().split()[-1])
        msg = bytes.fromhex(f.readline().split()[-1])
        pk_kat = bytes.fromhex(f.readline().split()[-1])
        sk_kat = bytes.fromhex(f.readline().split()[-1])
        siglen = int(f.readline().split()[-1])
        sig_kat = bytes.fromhex(f.readline().split()[-1])

        pk = [0]*g.CRYPTO_PUBLICKEYBYTES
        sk = [0]*g.CRYPTO_SECRETKEYBYTES
        crypto_sign_keypair(pk, sk, seed)
        assert bytes(pk) == pk_kat
        assert bytes(sk) == sk_kat


        sig = [0]*g.CRYPTO_BYTES
        crypto_sign_signature(sig, siglen, list(msg), mlen, sk)

        assert bytes(sig) == sig_kat

        assert crypto_sign_verify(sig, siglen, list(msg), mlen, pk) == 0

        f.readline()
        # print(f"Case {count} passed")
    print("Dilithium 5 passes all KATs")


if __name__ == "__main__":
    test_dilithium2()
    test_dilithium3()
    test_dilithium5()
