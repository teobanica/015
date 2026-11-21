"""padic.py -- p-adic valuations of the moments of the limiting laws,
congruences of Bell and involution numbers, and boundedness of the
binomial moments (Mahler test: bounded binomial moments <=> bounded
measure on Z_p; integral binomial moments, |b_n| <= 1 <=> measure of
norm <= 1, that is, Z_p-valued on the Z_p-valued continuous functions)."""
from fractions import Fraction
from math import comb, factorial

def v(x, p):
    x = Fraction(x)
    if x == 0:
        return 10**9
    a, b, k = x.numerator, x.denominator, 0
    while a % p == 0:
        a //= p; k += 1
    while b % p == 0:
        b //= p; k -= 1
    return k

def s(n, p):
    t = 0
    while n:
        t += n % p; n //= p
    return t

def dfact(k):          # (2k-1)!!
    r = 1
    for i in range(1, 2*k, 2):
        r *= i
    return r

def catalan(k):
    return comb(2*k, k) // (k+1)

def bell(n):
    row = [1]
    for _ in range(n):
        new = [row[-1]]
        for x in row:
            new.append(new[-1] + x)
        row = new
    return row[0]

def bells(n):
    return [bell(k) for k in range(n)]

def involutions(n):
    t = [1, 1]
    for k in range(2, n):
        t.append(t[-1] + (k-1)*t[-2])
    return t[:n]

def stirling1(n):      # signed s(n,k), x(x-1)...(x-n+1) = sum s(n,k) x^k
    S = [[0]*(n+1) for _ in range(n+1)]
    S[0][0] = 1
    for m in range(1, n+1):
        for k in range(1, m+1):
            S[m][k] = S[m-1][k-1] - (m-1)*S[m-1][k]
    return S

def binomial_moments(m, n):
    S = stirling1(n)
    return [sum(Fraction(S[j][k]) * m[k] for k in range(j+1)) / factorial(j)
            for j in range(n)]

def touchard(n, t):
    # moments of Poisson(t): sum_k S2(n,k) t^k
    S2 = [[0]*(n+1) for _ in range(n+1)]
    S2[0][0] = 1
    for m in range(1, n+1):
        for k in range(1, m+1):
            S2[m][k] = S2[m-1][k-1] + k*S2[m-1][k]
    return [sum(S2[j][k] * Fraction(t)**k for k in range(j+1)) for j in range(n)]

def narayana_moments(n, t):
    # moments of free Poisson(t): sum_{pi in NC(j)} t^{|pi|}
    out = [Fraction(1)]
    for j in range(1, n):
        out.append(sum(Fraction(comb(j, k) * comb(j, k-1), j) * Fraction(t)**k
                       for k in range(1, j+1)))
    return out

if __name__ == "__main__":
    log = open("padic.log", "w")
    def say(*a):
        line = " ".join(str(x) for x in a)
        print(line, flush=True); log.write(line + "\n")

    say("== (2k-1)!! : v_2 = 0, and v_p = (k+s_p(k)-s_p(2k))/(p-1) for p odd")
    ok = all(v(dfact(k), 2) == 0 for k in range(1, 300))
    for p in [3, 5, 7]:
        ok &= all(v(dfact(k), p) == (k + s(k, p) - s(2*k, p)) // (p-1) for k in range(1, 300))
    say("checked k<300:", ok)

    say("== Catalan: v_2(C_k) = s_2(k+1)-1")
    say("checked k<1000:", all(v(catalan(k), 2) == s(k+1, 2) - 1 for k in range(1000)))
    say("== Kummer: v_p(binom(2k,k)) = carries of k+k in base p = (2 s_p(k) - s_p(2k))/(p-1)")
    say("checked p=2,3,5,7, k<500:", all(v(comb(2*k, k), p) == (2*s(k, p) - s(2*k, p))//(p-1)
                                       for p in [2, 3, 5, 7] for k in range(500)))

    say("== Bell congruences (Touchard/Barsky): B_{n+k(p)p^{h-1}} = B_n mod p^h, k(p)=(p^p-1)/(p-1)")
    B = bells(400)
    for p in [3, 5]:
        kp = (p**p - 1)//(p - 1)
        for h in [1, 2]:
            T = kp * p**(h-1)
            if T + 40 < len(B):
                say(" p=%d h=%d period %d:" % (p, h, T),
                    all((B[n+T]-B[n]) % p**h == 0 for n in range(len(B)-T)))
    say(" p=2 mod 2, period 3:", all((B[n+3]-B[n]) % 2 == 0 for n in range(len(B)-3)))
    for h in [2, 3, 4, 5, 6]:
        T = 3*2**(h-1)
        say(" p=2 mod 2^%d, period 3*2^%d=%d (fails):" % (h, h-1, T),
            all((B[n+T]-B[n]) % 2**h == 0 for n in range(len(B)-T)))
        T = 3*2**h
        say(" p=2 mod 2^%d, period 3*2^%d=%d:" % (h, h, T),
            all((B[n+T]-B[n]) % 2**h == 0 for n in range(len(B)-T)))
    say(" B_n even iff n = 2 mod 3 (n<400):",
        all((B[n] % 2 == 0) == (n % 3 == 2) for n in range(400)))

    say("== involution numbers t_n = |P_12(n)| (moments of 1+N(0,1)), mod p^h periodicity, p odd")
    I = involutions(400)
    say(" p=2: I_0=%d, I_2=%d, so t_{n+2} = t_n mod 2 fails at n=0:" % (I[0], I[2]), (I[2]-I[0]) % 2 != 0)
    for p in [3, 5, 7]:
        for h in [1, 2]:
            per = None
            for T in range(1, 200):
                if all((I[n+T]-I[n]) % p**h == 0 for n in range(0, 400-T)):
                    per = T; break
            say(" p=%d h=%d period, from n=0:" % (p, h), per)

    say("== Catalan numbers are not moments of a bounded measure on Z_p:")
    say("   k=(p^m-1)/2 gives C_k a unit, k'=k+(p-1)p^(h-1), h=m-1, gives p | C_k', although")
    say("   x^(k+phi(p^h)) = x^k mod p^h on Z_p, for k>=h")
    for p in [3, 5, 7]:
        for m in [3, 4, 5]:
            k = (p**m - 1)//2; h = m - 1; k2 = k + (p-1)*p**(h-1)
            say("   p=%d m=%d: v_p(C_k)=%d v_p(C_k')=%d" % (p, m, v(catalan(k), p), v(catalan(k2), p)))
    for m in [4, 6, 8]:
        k = 2**m - 1; k2 = k + 2**(m-1)
        say("   p=2 m=%d: k=2^m-1, k'=k+2^(m-1): v_2(C_k)=%d v_2(C_k')=%d" % (m, v(catalan(k), 2), v(catalan(k2), 2)))
    say("== Poisson(t) binomial moments are t^n/n!, bounded iff v_p(t) >= 1/(p-1)")
    for p in [2, 3, 5]:
        for t in [1, p, p*p]:
            b = binomial_moments(touchard(40, t), 40)
            say("   p=%d t=%d: b_n = t^n/n! for n<40: %s, min v_p = %d" % (p, t,
                all(b[j] == Fraction(t**j, factorial(j)) for j in range(40)), min(v(x, p) for x in b)))
    say("== binomial moments b_n = L(binom(x,n)); bounded <=> bounded measure on Z_p,")
    say("   and min v_p(b_n) >= 0 <=> measure of norm <= 1 (integral)")
    n = 60
    gauss = [Fraction(dfact(k//2)) if k % 2 == 0 else Fraction(0) for k in range(n)]
    semi = [Fraction(catalan(k//2)) if k % 2 == 0 else Fraction(0) for k in range(n)]
    laws = {"Gaussian N(0,1)": gauss, "semicircle": semi,
            "Poisson(1)": touchard(n, 1), "free Poisson(1)": narayana_moments(n, 1)}
    for p in [2, 3, 5]:
        laws["Poisson(%d)" % p] = touchard(n, p)
        laws["free Poisson(%d)" % p] = narayana_moments(n, p)
    for name, m in laws.items():
        b = binomial_moments(m, n)
        for p in [2, 3, 5]:
            vals = [v(x, p) for x in b]
            say(" %-18s p=%d  min v_p(b_j), j<%d: %4d   v_p(b_j) at j=10,20,40,59: %s"
                % (name, p, n, min(vals), [vals[j] for j in (10, 20, 40, 59)]))
    say("== dilated Gaussian, sigma=25 at p=5: moments m_2k = 25^(2k) (2k-1)!!")
    n = 160
    g25 = [Fraction(25**k * dfact(k//2)) if k % 2 == 0 else Fraction(0) for k in range(n)]
    b = binomial_moments(g25, n)
    vals = [v(x, 5) for x in b]
    mn = min(vals)
    say("   max |b_n|_5 for n<%d: 5^%d, first attained at n=%d" % (n, -mn, vals.index(mn)))
    log.close()
