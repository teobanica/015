# Easiness over characteristic zero fields

**C. Opus and T. Banica** (November 2026)

This repository contains the paper (24 pages, PDF and LaTeX source) and the Python programs used for its computations.

## Summary

The easy quantum groups have intertwiners $T_\pi$, depending on set-theoretic partitions $\pi$, with integer entries, and so are defined over $\mathbb Q$. We study them over an arbitrary field $F$ of characteristic 0. Easiness itself does not depend on $F$, and the field comes into play via the $F$-forms, and via the analysis.

- **Easiness.** Over any $F$ the Hopf algebra is the base change of its rational version, the intertwiners are spanned by the maps $T_\pi$, and the Haar functional is $F$-valued, restricts to the rational form, and is given by the Weingarten formula when the Gram matrix is invertible. The Gram determinants of $S_N,O_N,S_N^+,O_N^+$ are given in closed form, with Beraha numbers as zeros in the free case.
- **Forms of $O_N^+$.** The $F$-forms of $O_N^+$ correspond to the central simple algebras of degree $N$ with orthogonal involution, the split ones being the Hopf algebras $B(E)$ with $E$ symmetric, up to similarity. Over $\mathbb Q$, the free bistochastic group $B_N^+$ is isomorphic to $O_{N-1}^+$ precisely when $N$ is an odd square, when $4|N$, or when $N\equiv2$ modulo 4 is a sum of 2 squares.
- **Free quadrics.** The quantum symmetry group of the free quadric given by $A\in GL_N(F)$ has Hopf algebra $B(A^{-1})$. For the free conics, the isomorphism class is determined by the similarity class of the symmetric part, and by one invariant $r\in F$, recovering the deformation parameters of the elliptic and hyperbolic free conics.
- **Forms of $S_N^+$.** The finite étale $F$-algebras $K$ produce $F$-forms of $S_N^+$, pairwise nonisomorphic for $N\geq3$, $N\neq6$. The quantum group $G_3\subset S_9^+$ coming from the color sector of the Standard Model has a $\mathbb Q$-form, with $K=\mathbb Q^3\times\mathbb Q(\omega)^3$, the area tensor making $K$ into a copy of $M_3(\mathbb Q)$.
- **Laws.** Moments, cumulants and free cumulants work over any $F$, with convolution powers indexed by $F$, and Poisson and central limit theorems over the complete valued fields.
- **p-adics.** The Poisson law $p_t$ is a bounded measure on $\mathbf Z_p$ precisely when $v_p(t)\geq1/(p-1)$, the Gaussian laws of nonzero variance never are, and neither are the laws of the main characters of $O_N^+$ with $N\geq2$ and $S_N^+$ with $N\geq4$, while those of $S_N$ are, with norms going to infinity. The paper ends with a list of open questions.

## Contents

- `015.pdf`, `015.tex`: the paper and its LaTeX source.
- Python programs, each writing a `.log` file, included here:

| File | Content | Section |
|---|---|---|
| `gram.py` | Gram determinants and their closed formulae, Weingarten matrices, Hankel determinants | 2, 6, 7 |
| `forms.py` | Hilbert symbols and Hasse-Minkowski, the rational forms $B_N^+$ of $O_{N-1}^+$ | 3 |
| `quadrics.py` | deformation parameters and congruence invariants of the free conics | 4 |
| `galois.py` | centralizers for étale algebras, descent of $G_3$, its 6 rational classical points | 5 |
| `padic.py` | valuations, Bell and involution congruences, binomial moments of the Poisson, Gaussian and free laws | 7 |

## Running the code

Python 3 is needed, with the standard library only. Each program runs from the top folder, in a few seconds, for instance:

```
python3 gram.py
```
