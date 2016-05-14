---
bibliography: bibliography/crypto.bib
title: Plaintext slot permutations in batch ciphertexts for NTRU-based fully homomorphic encryption schemes
author: Bogdan Kulynych
email: hello@bogdankulynych.me
---

\newpage

## Abstract:
Batching using CRT allows to naturally parallelize computation for fully homomorphic encryption schemes over the integers or polynomial rings. This technique has been widely used in literature to improve performance of homomorphic algorithms \cite{GHS12b} \cite{CCK13}, but not any algorithm can be computed over batch ciphertext that support only support addition and multiplication isomorphically. Additional ability to permute plaintext slots in a batch ciphertext means that any binary circuit can be expressed as an algorithm over batch ciphertexts, leading to homomorphic computation overhead polylogarithmic in the security parameter. Plaintext slot permutations implemented as Galois group actions have only been initially proposed and explicitly described for the BGV encryption scheme. This work considers applying the Galois field-based machinery for plaintext slots permutation to NTRU and YASHE fully homomorphic encryption schemes. We also pair it with double-CRT representation, that, combined, allows for highly efficient polylog overhead homomorphic computation with popular NTRU-based FHE schemes.
\newpage

\tableofcontents
\newpage

Introduction
------------

Since the introduction of the first FHE scheme in the seminal work by Gentry \cite{Gen09}, many new Somewhat Homomorphic Encryption (SHE) and Fully Homomorphic Encryption (FHE) designs and implementations aimed at making homomorphic encryption practical for real-world applications appeared. This work focuses on schemes from NTRU family: namely, construction proposed by Doröz et al. \cite{DHS15} and variant of the original homomorphic construction by López-Alt et al. \cite{LTV12} called YASHE \cite{BLLN13}. We consider batch ciphertext permutation techniques that are originally optimizations for the Brakerski-Gentry-Vaikuntanathan (BGV) scheme \cite{BGV12}.

Batching of ciphertexts, that is, ability to simultaneously execute homomorphic operations on multiple encrypted values in SIMD ("single instruction multiple data") fashion, is one of the many optimizations that have been proposed to different FHE schemes. It was first introduced in \cite{SV11} in the context of Gentry's original scheme, and the technique was later adapted to integer-based FHE schemes \cite{CCK13}, BGV schemes \cite{GHS12a}, and recently, NTRU \cite{RC14}. Notably, it was also shown in \cite{GHS12a} that it is possible to permute plaintext slots in a batch ciphertext. Such operation allows to perform homomorphic computation with overhead that is polylogarithmic in the security parameter by expressing the computation circuit in terms of batch additions, multiplications, and permutations.

Even though the techniques from \cite{GHS12a} are trivial to apply to schemes in algebraic setting similar to BGV, the batch permutation for NTRU is missing in the literature. It is not even mentioned in "Homomorphic AES evaluation using NTRU" \cite{DHS15}, which, as the name hints, follows "Homomorphic AES evaluation" \cite{GHS12b} where batch permutations in BGV ciphertexts are extensively used to optimize AES computation. This work considers plaintext slots permutation technique applied to the NTRU and YASHE schemes, i.e. we verify that Galois group actions permute plaintext slots in batch ciphertexts in NTRU setting as well.


Preliminaries
-------------

For all schemes in this work we use rings of dimension ``m`` defined by ``m``-th cyclotomic polynomials ``\Phi_m(X)``, ``\mathbb{A} = \mathbb{Z}[X]/\Phi_m(X)``. We let ``\mathbb{A}_q`` denote the set of elements of this ring reduced modulo ``q``. The ring ``\mathbb{A}`` is the ring of integers of the ``m``-th cyclotomic number field ``\mathbb{K} = \mathbb{Q}(\zeta_m)``, where ``\zeta_m`` is primitive ``m``-th root of unity. We let ``[a]_q`` for an element ``a \in \mathbb{A}`` denote the reduction of ``a`` modulo ``q``, fixing a unique representative in ``(-q/2, q/2]`` for its equivalence class. We fix a set of ``L`` primes ``p_0, p_1, ..., p_{L-1}``, and define ``t``-th ciphertext modulus ``q_t`` to be ``q_t = \prod\limits_{i = 0}^{L - t - 1} p_i``, obtaining a decreasing moduli chain ``q_0 > q_1 > ... > q_{L-1}``. The modulus of the fresh ciphertext is ``q_0 = p_0 \cdot p_1 \cdot ... \cdot p_{L-1}``, and it decreases down to ``q_{L-1} = p_0`` as more homomorphic operations are evaluated.

Let ``\lceil \cdot \rfloor`` denote rounding to the nearest integer, ``\lfloor \cdot \rfloor`` rounding down (floor). We define scaling factor ``\Delta_q`` as ``\Delta_q = \lfloor \frac{p}{q} \rfloor``.

### Plaintext slots

We use ``p`` to denote the plaintext space modulus, thus the plaintexts are elements of ``\mathbb{A}_p``. We assume that ``p`` does not divide ``m``. As a result, polynomial ``\Phi_m(X)`` factors modulo ``p`` into ``l`` irreducible factors, ``\Phi_m(X) = F_0(X) \times F_1(X) ... \times F_{l-1}(X) \mod p``. Since ``\mathbb{K}`` is Galois, all of the factors have the same degree ``d = \phi(m)/l``. Thus, the plaintext space splits in the product of ``l`` finite fields:
```math
\mathbb{A}_p \cong \mathbb{Z}[X]_p/F_0(X) \times \mathbb{Z}_p[X]/F_1(X) \times ... \times \mathbb{Z}_p[X]/F_{l-1}(X)
```

Each factor corresponds to a so-called plaintext slot. We view a polynomial ``a \in \mathbb{A}_p`` as representing an ``l``-dimensional vector ``(a \mod F_i)_{i=0}^{l-1}``, or, equivalently, ``l``-vector of elements in ``\mathbb{F}_{p^d}``. Chinese remainder theorem is used to pack the vector into a single value in ``\mathbb{A}_p``. We call this ciphertext representation _batch_, or CRT representation. We denote procedure of packing using CRT as ``\mathsf{CRT}: \mathbb{F}_{p^d}^l \rightarrow \mathbb{A}_p``, and similarly, unpacking as ``\mathsf{CRT}^{-1}: \mathbb{A}_p \rightarrow \mathbb{F}_{p^d}^l``.

### Sampling from ``\mathbb{Z}_q[X]/\Phi_m(X)``

Let us define following probability distributions:

- Gaussian ``\mathcal{DG}_q(\sigma^2)``. Draw ``a \in \mathbb{A}_q`` with coefficients from zero-mean Gaussian distribution ``\mathcal{N}(0, \sigma^2)`` rounded to the nearest integer.
- Small polynomial ``\mathcal{HWT}(h)``. For ``h < \phi(m)``, draw ``a \in \mathbb{A}`` with coefficients uniformly drawn from ``\{0, +1, -1\}``, such that ``a`` has exactly ``h`` non-zero entries.

We will denote as ``x \leftarrow \mathcal{D}`` drawing ``x`` from the distribution ``\mathcal{D}``.

Partial scheme definitions
--------------------------

We use variants and notation close to a single-framework review by Costache and Smart \cite{CN16}. For the sake of conciseness, we only provide partial scheme definitions consisting of key generation, encryption, and decryption procedures, ommiting other important components and details that are irrelevant in the scope of this work.

The following are procedures for NTRU and YASHE encryption schemes.

``\mathsf{KeyGen}( h ).`` \hangindent=2em Sample ``f, g \leftarrow \mathcal{HWT}(h)``. If ``f`` is not invertible in ``A_{q_0}``, resample ``f``. Otherwise, set private key ``\mathfrak{sk} = p \cdot f + 1``, and public key ``\mathfrak{pk} = [p \cdot g \cdot f^{-1}]_{q_0}``. Output ``(\mathfrak{sk}, \mathfrak{pk})``.

``\mathsf{NTRU.Encrypt}( \mathfrak{pk}, m \in \mathbb{A}_p ).`` \hangindent=2em Sample ``e_1, e_2 \leftarrow \mathcal{DG}_{q_0}(\sigma^2)``. Encrypt ``m``:
```math
c = [e_1 \cdot \mathfrak{pk} + p \cdot e_2 + m]_{q_0}
```
Set the ciphertext level ``t = 0``. Output ``\mathfrak{c} = (c, 0)``.

\par

``\mathsf{NTRU.Decrypt}( \mathfrak{sk}, \mathfrak{c} \in \mathbb{A}_{q_t} ).``
\hangindent=2em Output decrypted message:
```math
m' = \big[ [\mathfrak{sk} \cdot c ]_{q_t} \big]_p
```

\par


``\mathsf{YASHE.Encrypt}( \mathfrak{pk}, m \in \mathbb{A}_p ).`` \hangindent=2em Sample ``e_1, e_2 \leftarrow \mathcal{DG}_{q_0}(\sigma^2)``. Output:
```math
c = [e_1 \cdot \mathfrak{pk} + \Delta_{q_0} \cdot e_2 + m]_{q_0}
```
Set the ciphertext level ``t = 0``. Output ``\mathfrak{c} = (c, 0)``.


``\mathsf{YASHE.Decrypt}( \mathfrak{sk}, \mathfrak{c} \in \mathbb{A}_{q_t} ).``
\hangindent=2em Output decrypted message:
```math
m' = \Big[ \big\lceil \frac{p}{q_t} [\mathfrak{sk} \cdot c ]_{q_t} \big\rfloor \Big]_p
```

\par

Permutations of plaintext slots
-------------------------------

An algorithm that allows to perform arbitrary permutation of the plaintext slots using just the homomorphic ``\mathsf{Select}`` operation and cyclic rotations of plaintext slots is shown in \cite{GHS12a}. We verify that a method of performing cyclic rotations as automorphism of ``\mathbb{K}`` can be easily used in NTRU setting.

Recall that Galois group ``\mathcal{G}\mathsf{al}(\mathbb{K}/\mathbb{Q})`` action is a result of applying transformation ``\kappa_i: f(X) \mapsto f(X^i) \mod \Phi_m(X), q_t`` for ``i \in \mathbb{Z}_m^*``. Importantly, for some values ``i`` and a vector ``\mathbf{a} = (a_0, a_1, ..., a_l) \in \mathbb{F}_{p^d}^l`` with ``f = \mathsf{CRT}(\mathbf{a})``, the transformation ``f^{(i)} = \kappa_i(f)`` produces a polynomial such that ``\mathsf{CRT}^{-1}(f^{(i)}) = (a_{l-k-1}, ..., a_l, a_0, a_1, ..., a_{k})``. Namely, ``\kappa_i`` rotates the slots when ``i`` is not in ``\{p^k ~|~ k = 0, 1, ..., d-1\}``.

### Automorphisms
We will look at the effect of the Galois group action on decryption.

##### NTRU
If NTRU ciphertext ``\mathfrak{c}`` is decryptable, we have over ``\mathbb{Z}_{q_t}[X]/\Phi_m(X)``:
```math
[\mathfrak{sk} \cdot c]_{q_t} = m + p \cdot (e'_1 \cdot g + e'_2 + f \cdot m) + p^2 \cdot e'_2 \cdot f
```

That is, for some ``r, u, v \in \mathbb{A}_{q_t}`` the following equality in ``\mathbb{Z}_{q_t}[X]`` holds:
```math
\mathfrak{sk}(X) \cdot c(X) = m(X) + p \cdot u(X) + p^2 \cdot v(X) + r(X) \cdot \Phi_m(X)
```

If we apply ``\kappa_i`` to both sides of the equation, the equality will be preserved, since ``\kappa_i`` is an automorphism in ``\mathbb{K}``:
```math
\mathfrak{sk}(X^i) \cdot c(X^i) = m(X^i) + p \cdot u(X^i) + p^2 \cdot v(X^i) + r(X^i) \cdot \Phi_m(X^i)
```
It is easy to show that for ``i \in Z_m^*``, ``\Phi_m(X)`` divides ``\Phi_m(X^i)``, thus, over the ``\mathbb{Z}_{q_t}[X]/\Phi_m(X)``, ciphertext decrypts to ``m(X^i)`` under ``\mathfrak{sk}(X^i)``:
```math
\mathfrak{sk}(X^i) \cdot c(X^i) = m(X^i) + p \cdot u(X^i) + p^2 \cdot v(X^i)
```

\par

##### YASHE
We apply the same reasoning to YASHE ciphertexts. Recall that ``\Delta_{q_t} = \lfloor \frac{q_t}{p} \rfloor = \frac{q_t}{p} - \epsilon`` for some ``\epsilon`` in ``[0, 1)``. When YASHE ciphertext ``\mathfrak{c}`` is decryptable, we have for some ``e'_1, e'_2 \in \mathbb{A}_{q_t}`` over ``\mathbb{Z}_{q_t}[X]/\Phi_m(X)``:
```math
\frac{p}{q_t} \mathfrak{sk} \cdot c = \frac{p \cdot (\frac{q_t}{p} - \epsilon)}{q_t} + \frac{p \cdot e'_1}{q_t} + p \cdot e'_2 = m + p \cdot \frac{e'_1 - \epsilon \cdot m}{q_t} + p \cdot e'_2
```

In terms of polynomials, for ``u, v, r \in \mathbb{A}_{q_t}`` we have over ``\mathbb{Z}_{q_t}[X]``:
```math
\frac{p}{q_t} \mathfrak{sk}(X) \cdot c(X) = m(X) + p \cdot u(X) + \frac{p}{q_t} \cdot v(X) + r(X) \cdot \Phi_m(X)
```
After applying, ``\kappa_i`` we have over ``\mathbb{Z}_{q_t}[X]``:
```math
\frac{p}{q_t} \mathfrak{sk}(X^i) \cdot c(X^i) = m(X^i) + p \cdot u(X^i) + \frac{p}{q_t} \cdot v(X^i) + r(X^i) \cdot \Phi_m(X^i)
```

This again reduces to a decryption of ``m(X^i)`` under ``\mathfrak{sk^i}`` over ``\mathbb{Z}_{q_t}[X]/\Phi_m(X)``:
```math
\frac{p}{q_t} \mathfrak{sk}(X^i) \cdot c(X^i) = m(X^i) + p \cdot u(X^i) + \frac{p}{q_t} \cdot v(X^i)
```

\par


Conclusion
----------

Plaintext slot rotations from Galois group actions, and therefore plaintext slot permutations, can be used in NTRU-based homomorphic encryption schemes. This follows from the similarity of algebraic settings of the BGV as described in \cite{GHS12a} and NTRU.


References
----------
