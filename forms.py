"""forms.py -- forms of the free orthogonal quantum group over Q.
B_N^+ = B(C), with C the Gram matrix of the A_{N-1} root lattice, which over Q
is the diagonal form <1.2, 2.3, ..., (N-1)N>. We decide, by Hasse-Minkowski,
whether C is similar to the identity form, c*C = I over Q for some c in Q^*,
which by Theorem B is equivalent to B_N^+ = O_{N-1}^+ over Q."""
from fractions import Fraction
from itertools import product

def factor(n):
    n, f, d = abs(n), {}, 2
    while d*d <= n:
        while n % d == 0:
            f[d] = f.get(d, 0) + 1; n //= d
        d += 1
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f

def squarefree(x):
    x = Fraction(x)
    n = x.numerator * x.denominator
    s = -1 if n < 0 else 1
    for p, e in factor(n).items():
        if e % 2:
            s *= p
    return s

def legendre(a, p):
    r = pow(a % p, (p-1)//2, p)
    return -1 if r == p-1 else r

def split(a, p):
    k = 0
    while a % p == 0:
        a //= p; k += 1
    return k, a

def hilbert(a, b, p):
    a, b = squarefree(a), squarefree(b)
    if p == 0:   # the real place
        return -1 if (a < 0 and b < 0) else 1
    al, u = split(a, p)
    be, w = split(b, p)
    if p != 2:
        s = (-1)**(al*be*((p-1)//2))
        if be % 2: s *= legendre(u, p)
        if al % 2: s *= legendre(w, p)
        return s
    eps = lambda x: ((x-1)//2) % 2
    om = lambda x: ((x*x-1)//8) % 2
    return (-1)**((eps(u)*eps(w) + al*om(w) + be*om(u)) % 2)

def hasse(diag, p):
    s = 1
    for i in range(len(diag)):
        for j in range(i+1, len(diag)):
            s *= hilbert(diag[i], diag[j], p)
    return s

def primes_of(xs):
    ps = {2}
    for x in xs:
        x = Fraction(x)
        ps |= set(factor(x.numerator)) | set(factor(x.denominator))
    return sorted(ps)

def equivalent(d1, d2):
    """Hasse-Minkowski for diagonal forms over Q."""
    if len(d1) != len(d2):
        return False
    if squarefree(Fraction(eval('*'.join(str(x) for x in d1)))) != squarefree(Fraction(eval('*'.join(str(x) for x in d2)))):
        return False
    if sum(1 for x in d1 if x < 0) != sum(1 for x in d2 if x < 0):
        return False
    for p in [0] + primes_of(d1 + d2):
        if hasse(d1, p) != hasse(d2, p):
            return False
    return True

def similar_to_identity(d):
    """Is there c in Q^* with c*q = <1,...,1> ? c can be taken squarefree, built
    from -1 and the primes involved (other primes change the Hasse invariant
    at that prime only if they appear, and then the determinant/Hasse fail)."""
    n = len(d)
    ps = primes_of(d)
    for signs in product([0, 1], repeat=len(ps)):
        c = 1
        for p, e in zip(ps, signs):
            c *= p**e
        for sg in [1, -1]:
            if equivalent([sg*c*x for x in d], [1]*n):
                return sg*c
    return None


def predicted(N):
    if N % 2:
        return int(round(N**0.5))**2 == N
    if N % 4 == 0:
        return True
    return all(e % 2 == 0 for p, e in factor(N).items() if p % 4 == 3)

def trivial_direct(N):
    d = [k*(k+1) for k in range(1, N)]
    c = 1 if (N-1) % 2 == 0 else N
    return equivalent([c*x for x in d], [1]*(N-1))

if __name__ == "__main__":
    log = open("forms.log", "w")
    def say(*a):
        line = " ".join(str(x) for x in a)
        print(line, flush=True); log.write(line + "\n")
    say("== sanity: Hilbert symbol product formula, a,b in [-30,30]")
    ok = True
    for a in range(-30, 31):
        for b in range(-30, 31):
            if a and b:
                pr = hilbert(a, b, 0)
                for p in primes_of([a, b]):
                    pr *= hilbert(a, b, p)
                ok &= (pr == 1)
    say("  product formula holds:", ok)
    say("== B_N^+ over Q versus O_{N-1}^+ : C = <k(k+1)>, k=1..N-1, det C = N mod squares")
    triv = []
    for N in range(2, 41):
        d = [k*(k+1) for k in range(1, N)]
        c = similar_to_identity(d)
        status = "trivial form, c=%d" % c if c is not None else "nontrivial form"
        if c is not None:
            triv.append(N)
        say("  N=%2d  n=N-1=%2d  det=N:%3d  %s" % (N, N-1, N, status))
    say("  N<=40 with B_N^+ = O_{N-1}^+ over Q:", triv)
    say("== real forms: O^+(p,q), p+q=n, up to similarity (p,q)~(q,p); count floor(n/2)+1")
    for n in range(1, 9):
        classes = set()
        for p in range(n+1):
            classes.add(tuple(sorted((p, n-p))))
        say("  n=%d: %d classes" % (n, len(classes)))
    say("== Theorem: trivial iff (N odd and square) or (4|N) or (N=2 mod 4 and N sum of 2 squares)")
    say("  checked directly for 2<=N<=150 (c forced: c=1 for N odd, c=N for N even):",
        all(trivial_direct(N) == predicted(N) for N in range(2, 151)))
    log.close()
