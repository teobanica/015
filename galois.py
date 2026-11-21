"""galois.py -- Galois forms of the quantum permutation groups.
(1) The F-points of the classical version of A_K are Aut_F(K), the centralizer
in S_N of the image of the Galois group acting on Hom(K, Fbar). We compute
these centralizers for several etale Q-algebras K.
(2) The area tensor b on F_3^2, with coefficients in Q(w), is invariant under
the semilinear map tau = (permutation s(x,y)=(x,-y)) o (complex conjugation).
We check this, compute the fixed etale algebra K = Q^3 x Q(w)^3, and check
that (K, b) is a central simple algebra over Q with a nontrivial idempotent,
hence isomorphic to M_3(Q)."""
from itertools import permutations, product
from fractions import Fraction

def compose(a, b):
    return tuple(a[b[i]] for i in range(len(b)))

def generate(gens, n):
    G = {tuple(range(n))}
    frontier = list(G)
    while frontier:
        new = []
        for g in frontier:
            for s in gens:
                h = compose(s, g)
                if h not in G:
                    G.add(h); new.append(h)
        frontier = new
    return G

def centralizer(G, n):
    return [s for s in permutations(range(n)) if all(compose(s, g) == compose(g, s) for g in G)]

def cyc(n, *cycles):
    p = list(range(n))
    for c in cycles:
        for i in range(len(c)):
            p[c[i]] = c[(i+1) % len(c)]
    return tuple(p)

# Galois images on the roots (embeddings K -> Qbar), as permutation groups
EXAMPLES = {
    "Q^4": (4, []),
    "Q(sqrt2) x Q(sqrt3)": (4, [cyc(4, (0, 1)), cyc(4, (2, 3))]),
    "Q(sqrt2) x Q(sqrt2)": (4, [cyc(4, (0, 1), (2, 3))]),
    "Q(2^(1/3))": (3, [cyc(3, (0, 1, 2)), cyc(3, (1, 2))]),
    "Q(zeta_5)": (4, [cyc(4, (0, 1, 2, 3))]),
    # roots of x^4-2: a, i a, -a, -i a, Galois group D_4 = <rotation, complex conjugation>
    "Q(2^(1/4))": (4, [cyc(4, (0, 1, 2, 3)), cyc(4, (1, 3))]),
    "Q(2^(1/3)) x Q": (4, [cyc(4, (0, 1, 2)), cyc(4, (1, 2))]),
    "Q(w) x Q(w) x Q(w) x Q^3 (the color algebra)": (9, [cyc(9, (3, 4), (5, 6), (7, 8))]),
}

# arithmetic in Q(w), w^2 = -1 - w, elements (a, b) = a + b w
def qmul(x, y):
    a, b = x; c, d = y
    return (a*c - b*d, a*d + b*c - b*d)

def qadd(x, y):
    return (x[0] + y[0], x[1] + y[1])

def qconj(x):         # w -> w^2 = -1 - w
    a, b = x
    return (a - b, -b)

ONE, ZERO = (Fraction(1), Fraction(0)), (Fraction(0), Fraction(0))
WPOW = [ONE, (Fraction(0), Fraction(1)), (Fraction(-1), Fraction(-1))]

PTS = [(x, y) for x in range(3) for y in range(3)]
IDX = {p: i for i, p in enumerate(PTS)}

def det(a, c):
    return (a[0]*c[1] - a[1]*c[0]) % 3

def sub(a, c):
    return ((a[0]-c[0]) % 3, (a[1]-c[1]) % 3)

def area_tensor():
    # b(d_a x d_c) = 1/9 sum_r w^{-det(a-r,c-r)} d_r
    B = {}
    for a in PTS:
        for c in PTS:
            for r in PTS:
                e = (-det(sub(a, r), sub(c, r))) % 3
                B[a, c, r] = (WPOW[e][0]/9, WPOW[e][1]/9)
    return B

def s(p):
    return (p[0], (-p[1]) % 3)

if __name__ == "__main__":
    log = open("galois.log", "w")
    def say(*a):
        line = " ".join(str(x) for x in a)
        print(line, flush=True); log.write(line + "\n")
    say("== Aut_Q(K) = centralizer of the Galois image in S_N")
    for name, (n, gens) in EXAMPLES.items():
        G = generate(gens, n)
        Z = centralizer(G, n)
        say("  %-46s N=%d  |Galois image|=%3d  |Aut_Q(K)|=%d" % (name, n, len(G), len(Z)))
    say("== area tensor on F_3^2")
    B = area_tensor()
    inv = all(B[s(a), s(c), s(r)] == qconj(B[a, c, r]) for a in PTS for c in PTS for r in PTS)
    say("  b(s a, s c ; s r) = conj b(a, c ; r) for all a,c,r:", inv)
    rat = all(B[a, c, r][1] == 0 for a in PTS for c in PTS for r in PTS)
    say("  b has rational coefficients in the point basis:", rat)
    fixed = [p for p in PTS if s(p) == p]
    say("  fixed points of s:", fixed, " -> K = Q^%d x Q(w)^%d" % (len(fixed), (9-len(fixed))//2))
    # Fourier basis chi_g(r) = w^{det(g,r)}, with b(chi_g x chi_h) = w^{det(g,h)} chi_{g+h}
    def chi(g):
        return {r: WPOW[det(g, r)] for r in PTS}
    def bprod(f1, f2):
        out = {r: ZERO for r in PTS}
        for a in PTS:
            for c in PTS:
                coef = qmul(f1[a], f2[c])
                if coef != ZERO:
                    for r in PTS:
                        out[r] = qadd(out[r], qmul(coef, B[a, c, r]))
        return out
    ok = True
    for g in PTS:
        for h in PTS:
            lhs = bprod(chi(g), chi(h))
            gh = ((g[0]+h[0]) % 3, (g[1]+h[1]) % 3)
            rhs = {r: qmul(WPOW[det(g, h)], chi(gh)[r]) for r in PTS}
            ok &= (lhs == rhs)
    say("  b(chi_g x chi_h) = w^det(g,h) chi_{g+h}, a twisted group algebra = M_3:", ok)
    tau = lambda f: {r: qconj(f[s(r)]) for r in PTS}
    say("  tau(chi_g) = chi_{s g} for all g:", all(tau(chi(g)) == chi(s(g)) for g in PTS))
    e = {r: ZERO for r in PTS}
    for x in range(3):
        for r in PTS:
            e[r] = qadd(e[r], (chi((x, 0))[r][0]/3, chi((x, 0))[r][1]/3))
    say("  e = (chi_0 + chi_(1,0) + chi_(2,0))/3 lies in K:", tau(e) == e)
    say("  e is a b-idempotent, e != 0, 1:", bprod(e, e) == e and e != chi((0, 0)) and any(v != ZERO for v in e.values()))
    # Q-points of the classical version: affine symplectic maps commuting with s
    def aff(M, vv):
        return lambda p: ((M[0][0]*p[0]+M[0][1]*p[1]+vv[0]) % 3, (M[1][0]*p[0]+M[1][1]*p[1]+vv[1]) % 3)
    cnt = tot = 0
    for M in product(range(3), repeat=4):
        M = ((M[0], M[1]), (M[2], M[3]))
        if (M[0][0]*M[1][1]-M[0][1]*M[1][0]) % 3 != 1:
            continue
        for vv in PTS:
            f = aff(M, vv); tot += 1
            if all(f(s(p)) == s(f(p)) for p in PTS):
                cnt += 1
    say("  |F_3^2 x SL_2(3)| =", tot, "; elements commuting with s (Q-points of the classical version of the Q-form):", cnt)
    log.close()
