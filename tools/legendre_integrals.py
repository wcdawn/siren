import sympy
import sympy.functions
import sympy.functions.special
from sympy.functions.special.polynomials import legendre


def kron(m, n):
    if m == n:
        return 1
    return 0


def factorial(n):
    if n == 0:
        return 1
    return n * factorial(n - 1)


def binomial(n, k):
    return factorial(n) // (factorial(k) * factorial(n - k))


sympy.init_printing()

x = sympy.Symbol("x")

print(legendre(3, x))

NMAX = 8
for n in range(1, NMAX, 2):
    for m in range(NMAX):
        print("n=", n, "m=", m)

        if m % 2 == 1:
            integral = sympy.Integer(kron(m, n)) / sympy.Integer(2 * n + 1)
        else:
            integral = (
                (-1) ** ((m + n - 1) // 2)
                * factorial(n)
                * factorial(m)
                / sympy.Integer(
                    2 ** (m + n - 1)
                    * (m - n)
                    * (m + n + 1)
                    * factorial(m // 2) ** 2
                    * factorial((n - 1) // 2) ** 2
                )
            )

        pos = (2 * m + 1) * integral  # \mu \in [0, 1]
        neg = 2 * sympy.Integer(kron(m, n)) - pos  # \mu \in [-1, 0]
        print(pos, neg)

n = sympy.symbols("n", integer=True, positive=True)
