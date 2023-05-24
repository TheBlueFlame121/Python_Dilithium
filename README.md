# Python_Dilithium

A Python translation of the reference Dilithium implementation for easy experimentation. The goal of this project was to have an implementation that is as close to the latest reference (at the time of writing, v3.1) implementation of Dilithium while offering the flexibility to easily experiment and change stuff around.

Some features of this implementation:
1. Effort has been made to keep code as close to reference as possible.
2. This implementation passes all KATs generated from the reference.
3. Better handling of global parameters to switch modes on the go.
4. Key-generation accepts optional seed input for testing

## Disclaimer

This code is not suitable for real world deployment under any scenario. It is neither performant nor constant time. It should only be used for education and experimentation.

## Code Example

```python
>>> from sign import *
>>>
>>> # Set mode
>>> g.set_mode(2) # Can be 2, 3 or 5
>>>
>>> # Key generation
>>> pk = [0]*g.CRYPTO_PUBLICKEYBYTES
>>> sk = [0]*g.CRYPTO_SECRETKEYBYTES
>>> crypto_sign_keypair(pk, sk)
0
>>> # Signing
>>> msg = b'This is just a test message'
>>> sig = [0]*g.CRYPTO_BYTES
>>> crypto_sign_signature(sig, len(sig), list(msg), len(msg), sk)
0
>>> # Verification
>>> crypto_sign_verify(sig, len(sig), list(msg), len(msg), pk) == 0
True
>>> # Sanity checks
>>> crypto_sign_verify(sig, len(sig), list(urandom(len(msg))), len(msg), pk) == -1
True
>>> pk2, sk2 = [0]*g.CRYPTO_PUBLICKEYBYTES, [0]*g.CRYPTO_SECRETKEYBYTES
>>> crypto_sign_verify(sig, len(sig), list(urandom(len(msg))), len(msg), pk2) == -1
True

```

## Additional Details

### KATs

This implementation passes all the test vectors generated from v3.1 of the reference implementation. Since there were [major changes](https://github.com/pq-crystals/dilithium/commit/e989e691ae3d3f5933d012ab074bdc413ebc6fad) between versions 3 and 3.1, I had to generate the KATs myself.

The KAT files themselves can be found in the ['KATs'](KATs) folder while the test functions can be found in ['test.py'](test.py).

### Dependencies

Dilithium uses `Shake256` XOF internally in many places. It didn't seem worthwhile to reimplement it from scratch so I make use of the `pycryptodome` library for it.

You can install it by running `pip -r install requirements`.

### Global parameters

The reference C implementation puts all the parameters in one file called `params.h` and then imports it everywhere. Functions use the value of these global variables as they were at the time of import. While this is fine for a compiled language as one would only make changes to the code before compiling again, it creates problems for an interpreted language like Python.

Keeping all the parameters in one file would mean that we would not be able to change modes while in an interactive session or in a script. To change modes one would have to exit, modify the file and then run it again. To overcome this, instead of treating these parameters as global variables, I made them members of a class.

The class is still defined in params.py and an object `g` of it has been instantiated. The change we have to make now is that to access a global variable `X`, we have to use `g.X`. The upside is that I was able to write another member function called `set_mode` which allows us to change between the different parameter sets on the go.

So in short, the file params differs from the reference in some non trivial ways and all the other files use global variables with a prefix `g.`.

### Seeded key generation

The key generation algorithm `crypto_sign_keypair` accepts an additional argument `det`. It has to be 32 bytes and is used as a seed for deterministic key generation. If it is not provided, the algorithm uses `urandom(32)` as the seed. This is useful for testing and was used for verification against the KATs.

### Benchmarks

Because everyone needs numbers:

| 500 Iterations            | Mode 2 | Mode 3 | Mode 5 |
|---------------------------|--------|--------|--------|
| key generation            | 0.013s | 0.019s | 0.03 s |
| signing 32 byte message   | 0.323s | 0.224s | 0.309s |
| verifying 32 byte message | 0.015s | 0.022s | 0.034s |

The data was generated using ['benchmark.py'](benchmark.py). Note that the signing times are volatile due to rejection sampling. The test was done on an i7-10750H laptop cpu.
