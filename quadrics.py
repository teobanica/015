"""quadrics.py -- the deformation parameter of the quantum symmetry group of a
free quadric f = sum_ij H_ij x_i x_j - 1, namely q + 1/q = -tr(H^{-1} H^t),
checked against the parameters of the elliptic and hyperbolic conics, and,
for H > 0, against the matrix F = L^{-t} L^* with H = L^* L."""
import cmath, math, random

def mat(n, f):
    return [[f(i, j) for j in range(n)] for i in range(n)]

def mul(A, B):
    return [[sum(A[i][k]*B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]

def tr(A):
    return sum(A[i][i] for i in range(len(A)))

def T(A):
    return [list(r) for r in zip(*A)]

def adj(A):
    return [[x.conjugate() for x in r] for r in T(A)]

def inv(A):
    n = len(A)
    M = [list(map(complex, r)) + [complex(i == j) for j in range(n)] for i, r in enumerate(A)]
    for c in range(n):
        piv = max(range(c, n), key=lambda r: abs(M[r][c]))
        M[c], M[piv] = M[piv], M[c]
        pv = M[c][c]
        M[c] = [x/pv for x in M[c]]
        for r in range(n):
            if r != c:
                f = M[r][c]
                M[r] = [x - f*y for x, y in zip(M[r], M[c])]
    return [r[n:] for r in M]

def cholesky_upper(H):
    # H = L^* L with L upper triangular
    n = len(H)
    L = [[0j]*n for _ in range(n)]
    for i in range(n):
        s = H[i][i] - sum(abs(L[k][i])**2 for k in range(i))
        L[i][i] = complex(math.sqrt(s.real))
        for j in range(i+1, n):
            L[i][j] = (H[i][j] - sum(L[k][i].conjugate()*L[k][j] for k in range(i))) / L[i][i].conjugate()
    return L

def t_of(H):
    return tr(mul(inv(H), T(H)))

def q_of(t):
    # root of q^2 + t q + 1 = 0
    d = cmath.sqrt(t*t - 4)
    return (-t + d)/2, (-t - d)/2

if __name__ == "__main__":
    log = open("quadrics.log", "w")
    def say(*a):
        line = " ".join(str(x) for x in a)
        print(line, flush=True); log.write(line + "\n")
    random.seed(1)
    say("== elliptic conic: H = [[1, i l],[-i l, 1]], q = (l-1)/(l+1), check q + 1/q = -t")
    for l in [0.0, 0.3, 0.7, 1.5, 3.0]:
        H = [[1, 1j*l], [-1j*l, 1]]
        t = t_of(H)
        q = (l-1)/(l+1)
        say("  l=%.2f  t=%.6f  -(q+1/q)=%.6f  t real: %s" % (l, t.real, -(q+1/q), abs(t.imag) < 1e-12))
    say("== hyperbolic conic: H = [[1, i m],[-i m, -1]], m = tan(theta), check t = 2cos(2 theta), q = -e^{+-2 i theta}")
    for th in [0.0, 0.2, 0.6, 1.2]:
        m = math.tan(th)
        H = [[1, 1j*m], [-1j*m, -1]]
        t = t_of(H)
        qs = q_of(t)
        say("  theta=%.2f  t=%.6f  2cos(2theta)=%.6f  q=%s  -e^{2i theta}=%s"
            % (th, t.real, 2*math.cos(2*th), [complex(round(z.real, 6), round(z.imag, 6)) for z in qs],
               complex(round((-cmath.exp(2j*th)).real, 6), round((-cmath.exp(2j*th)).imag, 6))))
    say("== elliptic quadric in N=2n(+1) variables: t = sum_k -(q_k+1/q_k) (+1 if N odd)")
    for N in [4, 5, 6, 7]:
        n = N//2
        ls = [random.uniform(0, 3) for _ in range(n)]
        H = [[0j]*N for _ in range(N)]
        for i in range(N):
            H[i][i] = 1
        for k, l in enumerate(ls):
            H[2*k][2*k+1] = 1j*l; H[2*k+1][2*k] = -1j*l
        t = t_of(H)
        pred = sum(-((l-1)/(l+1) + (l+1)/(l-1)) for l in ls) + (1 if N % 2 else 0)
        say("  N=%d  t=%.6f  predicted=%.6f" % (N, t.real, pred))
    say("== H > 0 random: H = L^*L, F = L^{-t} L^*, check F Fbar = 1 and tr(F^*F) = t")
    for N in [2, 3, 4, 6]:
        A = mat(N, lambda i, j: complex(random.gauss(0, 1), random.gauss(0, 1)))
        H = mul(adj(A), A)
        for i in range(N):
            H[i][i] += N
        L = cholesky_upper(H)
        err = max(abs(x - y) for r1, r2 in zip(mul(adj(L), L), H) for x, y in zip(r1, r2))
        F = mul(inv(T(L)), adj(L))
        Fbar = [[x.conjugate() for x in r] for r in F]
        FFb = mul(F, Fbar)
        e1 = max(abs(FFb[i][j] - (i == j)) for i in range(N) for j in range(N))
        say("  N=%d  |L^*L-H|=%.1e  |F Fbar - 1|=%.1e  tr(F^*F)=%.6f  t=%.6f"
            % (N, err, e1, tr(mul(adj(F), F)).real, t_of(H).real))
    say("== free conics over Q: A = S + k J, d = det S, r = k^2/d, check t = tr(A^{-1}A^t) = 2(1-r)/(1+r), exactly")
    from fractions import Fraction as Fr
    ok = True
    for _ in range(2000):
        a, b, c, k = [Fr(random.randint(-9, 9), random.randint(1, 5)) for _ in range(4)]
        d = a*c - b*b
        if d == 0 or d + k*k == 0:
            continue
        A = [[a, b + k], [b - k, c]]
        dA = A[0][0]*A[1][1] - A[0][1]*A[1][0]
        Ai = [[A[1][1]/dA, -A[0][1]/dA], [-A[1][0]/dA, A[0][0]/dA]]
        At = [[A[0][0], A[1][0]], [A[0][1], A[1][1]]]
        t = sum(Ai[i][j]*At[j][i] for i in range(2) for j in range(2))
        r = k*k/d
        ok &= (t == 2*(1-r)/(1+r))
    say("  2000 random rational A:", ok)
    say("== invariance of r under A -> c P^t A P, P random rational")
    ok = True
    for _ in range(500):
        a, b, c, k, cc = [Fr(random.randint(-9, 9) or 1, random.randint(1, 5)) for _ in range(5)]
        d = a*c - b*b
        if d == 0:
            continue
        P = [[Fr(random.randint(-5, 5)) for _ in range(2)] for _ in range(2)]
        if P[0][0]*P[1][1] - P[0][1]*P[1][0] == 0:
            continue
        A = [[a, b + k], [b - k, c]]
        B = [[cc*sum(P[x][i]*A[x][y]*P[y][j] for x in range(2) for y in range(2)) for j in range(2)] for i in range(2)]
        S2 = [[(B[i][j] + B[j][i])/2 for j in range(2)] for i in range(2)]
        k2 = (B[0][1] - B[1][0])/2
        d2 = S2[0][0]*S2[1][1] - S2[0][1]**2
        ok &= (k2*k2/d2 == k*k/d)
    say("  500 random congruences:", ok)
    log.close()
