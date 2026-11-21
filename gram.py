"""gram.py -- Gram matrices (with closed formulas for their determinants) G_N(pi,sigma) = N^{|pi v sigma|} for the categories
P, NC (permutation case, k points) and P_2, NC_2 (orthogonal case, 2k points).
Computes det G_N as an integer polynomial in N, its integer roots, and the
p-adic valuations of the Weingarten entries, W = G^{-1}."""
from fractions import Fraction

def set_partitions(n):
    if n == 0:
        yield []
        return
    for p in set_partitions(n-1):
        for i in range(len(p)):
            yield p[:i] + [p[i] + [n-1]] + p[i+1:]
        yield p + [[n-1]]

def noncrossing(p):
    lab = {}
    for i, b in enumerate(p):
        for x in b:
            lab[x] = i
    n = len(lab)
    for a in range(n):
        for b in range(a+1, n):
            for c in range(b+1, n):
                for d in range(c+1, n):
                    if lab[a] == lab[c] and lab[b] == lab[d] and lab[a] != lab[b]:
                        return False
    return True

def category(name, k):
    if name in ("P", "NC"):
        ps = list(set_partitions(k))
    else:
        ps = [p for p in set_partitions(2*k) if all(len(b) == 2 for b in p)]
    if name in ("NC", "NC2"):
        ps = [p for p in ps if noncrossing(p)]
    return ps

def join_blocks(p, q):
    lab = {}
    for b in p:
        for x in b:
            lab[x] = b[0]
    parent = {x: x for x in lab}
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for b in p + q:
        for x in b[1:]:
            parent[find(x)] = find(b[0])
    return len({find(x) for x in lab})

def gram_exponents(ps):
    return [[join_blocks(p, q) for q in ps] for p in ps]

def det(M):
    M = [[Fraction(x) for x in row] for row in M]
    n, d = len(M), Fraction(1)
    for c in range(n):
        piv = next((r for r in range(c, n) if M[r][c] != 0), None)
        if piv is None:
            return Fraction(0)
        if piv != c:
            M[c], M[piv] = M[piv], M[c]; d = -d
        d *= M[c][c]
        for r in range(c+1, n):
            f = M[r][c] / M[c][c]
            if f:
                for j in range(c, n):
                    M[r][j] -= f * M[c][j]
    return d

def inverse(M):
    n = len(M)
    A = [[Fraction(x) for x in row] + [Fraction(int(i == j)) for j in range(n)] for i, row in enumerate(M)]
    for c in range(n):
        piv = next(r for r in range(c, n) if A[r][c] != 0)
        A[c], A[piv] = A[piv], A[c]
        pv = A[c][c]
        A[c] = [x / pv for x in A[c]]
        for r in range(n):
            if r != c and A[r][c] != 0:
                f = A[r][c]
                A[r] = [x - f*y for x, y in zip(A[r], A[c])]
    return [row[n:] for row in A]

def interpolate(xs, ys):
    # Newton interpolation, returns coefficient list (low degree first)
    n = len(xs)
    coef = [Fraction(y) for y in ys]
    for j in range(1, n):
        for i in range(n-1, j-1, -1):
            coef[i] = (coef[i] - coef[i-1]) / (xs[i] - xs[i-j])
    poly = [Fraction(0)]
    for i in range(n-1, -1, -1):
        # poly = poly*(x - xs[i]) + coef[i]
        new = [Fraction(0)] * (len(poly) + 1)
        for d, c in enumerate(poly):
            new[d+1] += c
            new[d] -= c * xs[i]
        new[0] += coef[i]
        poly = new
    while len(poly) > 1 and poly[-1] == 0:
        poly.pop()
    return poly

def evalp(poly, x):
    r = Fraction(0)
    for c in reversed(poly):
        r = r * x + c
    return r

def divide_linear(poly, a):
    # divide by (x - a), assumes exact
    n = len(poly) - 1
    q = [Fraction(0)] * n
    r = poly[n]
    for d in range(n-1, -1, -1):
        q[d] = r
        r = poly[d] + r * a
    return q

def factor_integer_roots(poly, rng=range(-6, 20)):
    roots = {}
    for a in rng:
        while len(poly) > 1 and evalp(poly, a) == 0:
            poly = divide_linear(poly, a)
            roots[a] = roots.get(a, 0) + 1
    return roots, poly

def v(x, p):
    x = Fraction(x)
    if x == 0:
        return None
    a, b, k = x.numerator, x.denominator, 0
    while a % p == 0:
        a //= p; k += 1
    while b % p == 0:
        b //= p; k -= 1
    return k

from math import comb, factorial
def bareiss(M):
    M=[row[:] for row in M]; n=len(M); sign=1; prev=1
    for k in range(n-1):
        if M[k][k]==0:
            for r in range(k+1,n):
                if M[r][k]!=0: M[k],M[r]=M[r],M[k]; sign=-sign; break
            else: return 0
        for i in range(k+1,n):
            for j in range(k+1,n):
                M[i][j]=(M[i][j]*M[k][k]-M[i][k]*M[k][j])//prev
        prev=M[k][k]
    return sign*M[n-1][n-1]
def partitions_int(k,maxp=None):
    if maxp is None: maxp=k
    if k==0: yield []; return
    for p in range(min(k,maxp),0,-1):
        for rest in partitions_int(k-p,p): yield [p]+rest
def syt(shape):
    n=sum(shape); h=1
    conj=[sum(1 for r in shape if r>j) for j in range(shape[0])] if shape else []
    for i,r in enumerate(shape):
        for j in range(r): h*= (r-j-1)+(conj[j]-i-1)+1
    return factorial(n)//h
def U(m,x):  # Chebyshev 2nd kind, U_0=1,U_1=x
    a,b=1,x
    if m==0: return a
    for _ in range(m-1): a,b=b,x*b-a
    return b
def a(k,m): return comb(2*k,k-m)-2*(comb(2*k,k-m-1) if k-m-1>=0 else 0)+(comb(2*k,k-m-2) if k-m-2>=0 else 0)
def ff(N,b):
    r=1
    for i in range(b): r*=N-i
    return r
def f_S(N,k):
    r=1
    for p in category('P',k): r*=ff(N,len(p))
    return r
def f_O(N,k):
    r=1
    for lam in partitions_int(k):
        v=1
        for i,row in enumerate(lam,1):
            for j in range(1,row+1): v*=N+2*j-i-1
        r*=v**syt([2*x for x in lam])
    return r
def f_Op(N,k):
    r=1
    for m in range(1,k+1): r*=U(m,N)**a(k,m)
    return r
def f_Sp_sq(q,k):  # det for S_{q^2}^+, k points, q integer
    C=comb(2*k,k)//(k+1)
    r=q**C
    for m in range(1,k+1): r*=U(m,q)**a(k,m)
    return r

def pmul(a, b):          # product of integer polynomials, coefficient lists
    r = [0]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                r[i+j] += x*y
    return r
def pdivexact(a, b):      # quotient a/b, and whether the division is exact in Z[N]
    a = a[:]; q = [0]*(len(a)-len(b)+1)
    for i in range(len(a)-len(b), -1, -1):
        c, rem = divmod(a[i+len(b)-1], b[-1])
        if rem:
            return None, False
        q[i] = c
        for j in range(len(b)):
            a[i+j] -= c*b[j]
    return q, all(x == 0 for x in a)
def Upoly(m):             # coefficients of the Chebyshev polynomial U_m
    a, b = [1], [0, 1]
    if m == 0:
        return a
    for _ in range(m-1):
        nb = pmul([0, 1], b)
        nb = [(nb[i] if i < len(nb) else 0) - (a[i] if i < len(a) else 0) for i in range(max(len(nb), len(a)))]
        a, b = b, nb
    return b
if __name__ == "__main__":
    log = open("gram.log", "w")
    def say(*a):
        line = " ".join(str(x) for x in a)
        print(line, flush=True); log.write(line + "\n")
    for name, ks in [("NC", [2, 3, 4, 5]), ("P", [2, 3, 4]), ("NC2", [2, 3, 4]), ("P2", [2, 3])]:
        for k in ks:
            ps = category(name, k)
            E = gram_exponents(ps)
            deg = sum(max(row) for row in E) + 2
            xs = list(range(deg + 1))
            ys = [det([[Fraction(N)**e for e in row] for row in E]) for N in xs]
            poly = interpolate(xs, ys)
            roots, rest = factor_integer_roots(poly)
            say("%-3s k=%d  |D(k)|=%d  deg det=%d  integer roots {root: mult} = %s  remaining factor degree %d, leading coeff %s"
                % (name, k, len(ps), len(poly)-1, dict(sorted(roots.items())), len(rest)-1, rest[-1]))
    say("== Weingarten for S_N^+ at k=2, check W = 1/(N^2(N-1)) [[N^2,-N],[-N,N]], one-block partition first")
    for N in [2, 3, 4, 5, 7]:
        ps = category("NC", 2)
        E = gram_exponents(ps)
        W = inverse([[Fraction(N)**e for e in row] for row in E])
        expect = [[Fraction(N*N, N*N*(N-1)), Fraction(-N, N*N*(N-1))], [Fraction(-N, N*N*(N-1)), Fraction(N, N*N*(N-1))]]
        # ordering of NC(2): [[0,1]] (one block) and [[0],[1]]
        order = [len(p) for p in ps]
        assert order == [1, 2]
        say(" N=%d blocks order %s  W=%s  matches: %s" % (N, order, [[str(x) for x in r] for r in W], W == expect))
    say("== minimal p-adic valuation of Weingarten entries, S_N^+ (NC) and S_N (P), k=4")
    for name in ["NC", "P"]:
        ps = category(name, 4)
        E = gram_exponents(ps)
        for N in [4, 5, 6, 8, 9, 10, 12]:
            G = [[Fraction(N)**e for e in row] for row in E]
            if det(G) == 0:
                say(" %s N=%d singular" % (name, N)); continue
            W = inverse(G)
            mins = {p: min(v(x, p) for row in W for x in row if x != 0) for p in [2, 3, 5, 7, 11]}
            say(" %s k=4 N=%d  min v_p(W) = %s" % (name, N, mins))
    say("== Hankel determinants: Catalan numbers (law of chi for S_N^+, N>=4) versus fixed points of S_N")
    from math import comb, factorial
    cat = [comb(2*k, k)//(k+1) for k in range(40)]
    say("  det(C_{i+j})_{i,j<r}, r=1..12:", [int(det([[cat[i+j] for j in range(r)] for i in range(r)])) for r in range(1, 13)])
    for N in range(2, 9):
        D = [1, 0]
        for n in range(2, N+1):
            D.append((n-1)*(D[-1] + D[-2]))
        law = [Fraction(comb(N, j)*D[N-j], factorial(N)) for j in range(N+1)]
        m = [sum(c*j**k for j, c in enumerate(law)) for k in range(2*N+6)]
        rank = max(r for r in range(1, N+3) if det([[m[i+j] for j in range(r)] for i in range(r)]) != 0)
        say("  S_%d: atoms %s, Hankel rank %d" % (N, [j for j, c in enumerate(law) if c], rank))
    say('== closed formulas: Lindstrom (S_N), zonal (O_N), Di Francesco (O_N^+), fattening (S_N^+, checked at N=4,9,16)')
    for name,K,f,Ns in [('P',[2,3,4,5],f_S,range(2,8)),('P2',[2,3,4],f_O,range(2,8)),('NC2',[2,3,4,5],f_Op,range(2,8)),('NC',[2,3,4,5,6],None,range(2,5))]:
        for k in K:
            ps=category(name,k); E=gram_exponents(ps); ok=True
            for N in Ns:
                if name=='NC':
                    NN=N*N; val=f_Sp_sq(N,k)
                else:
                    NN=N; val=f(N,k)
                d=bareiss([[NN**e for e in row] for row in E])
                ok&=(d==val)
            say('  %-3s %2d points  |D|=%3d  closed formula = exact determinant at several N: %s' % (name, k if name in ('P','NC') else 2*k, len(ps), ok))
    say('== O_N^+ at 12 points (k=6), Di Francesco formula = exact determinant, N=2,3,4 (a few seconds)')
    ps = category('NC2', 6); E = gram_exponents(ps)
    say('  NC2 12 points |D|=%d: %s' % (len(ps), all(bareiss([[N**e for e in row] for row in E]) == f_Op(N, 6) for N in [2, 3, 4])))
    say('== Di Francesco exponents a_km can be negative, but prod_m U_m(N)^a_km is in Z[N]')
    for k in range(1, 9):
        num = [1]; den = [1]
        for m in range(1, k+1):
            e = a(k, m); UP = Upoly(m)
            for _ in range(abs(e)):
                if e > 0: num = pmul(num, UP)
                else: den = pmul(den, UP)
        q, ok = pdivexact(num, den)
        say('  k=%d exponents %s  polynomial in Z[N]: %s' % (k, [a(k, m) for m in range(1, k+1)], ok))
    log.close()
