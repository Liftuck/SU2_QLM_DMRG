import scipy.special
import scipy.optimize
from scipy.integrate import quad

from scipy.special import mathieu_sem

import numpy as np

import math


def P(lmbd):
    # print(scipy.special.i0(4.0*lmbd))
    # print(scipy.special.i1(4.0*lmbd))
    # return (1.0/lmbd) + (2.0*scipy.special.iv(0,4.0*lmbd)/scipy.special.iv(1,4.0*lmbd))
    return 2.0 * scipy.special.iv(2, 4.0 * lmbd) / scipy.special.iv(
        1, 4.0 * lmbd)


def mathieu_b_large_q(n, q):
    h = np.sqrt(-q)
    s = 2 * (n - 1) + 1

    out = -2 * h * h
    out += 2 * s * h
    out -= 1 / 8 * (s * s + 1)
    out -= 1 / (2**7 * h) * (s**3 + 3 * s)
    out -= 1 / (2**12 * h**2) * (5 * s**4 + 34 * s**2 + 9)

    return out


def Pexc(lmbd):
    plaq = P(lmbd)
    top = scipy.integrate.quad(
        lambda t: np.sqrt(1.0 - 0.25 * t**2) * (t - plaq) *
        (t - plaq) * t * np.exp(2 * lmbd * t), -2, 2)[0]
    bot = scipy.integrate.quad(
        lambda t: np.sqrt(1.0 - 0.25 * t**2) * (t - plaq) *
        (t - plaq) * np.exp(2 * lmbd * t), -2, 2)[0]

    return top / bot


def eps0(lmbd, xi):
    return ((((3.0 * lmbd) / (8.0 * xi)) - (4.0 * xi)) * P(lmbd))


def N11(lmbd, xi):
    return (3 * P(lmbd)) / (2.0 * lmbd)


def D11(lmbd, xi):
    return (4.0 - (3.0 / (2.0 * lmbd)) * P(lmbd)) - (P(lmbd)**2)


def eps1(lmbd, xi, lmbd2):
    eps0(lmbd, xi) * P(lmbd) * lmbd2


def mass(lmbd, xi, beta):
    1


def getPlaquettes(beta):
    g2 = 4 / beta / np.sqrt(2)
    q = -8.0 / g2 / g2

    # q = -0.5 * beta * beta * 2

    def psi0(omega):
        return mathieu_sem(2, q,
                           (180 / np.pi) * omega / 4)[0] / np.sin(omega / 2)

    def psi1(omega):
        return mathieu_sem(4, q,
                           (180 / np.pi) * omega / 4)[0] / np.sin(omega / 2)

    p = [0.0, 0.0]

    for i, wv in enumerate([psi0, psi1]):
        i1 = lambda omega: (np.sin(omega / 2) * np.abs(wv(omega)))**2
        i2 = lambda omega: (np.cos(omega / 2)) * i1(omega)
        try:
            N = quad(i1, 0, 2 * np.pi)[0]
            p[i] = quad(i2, 0, 2 * np.pi)[0] / N
        except:
            p[i] = float("nan")

    return p


def getEnergies(beta, k=2):
    g2 = 4 / beta / np.sqrt(2)

    def energy(n):
        if 1 / beta > 0.05:
            out = (g2 / 8) * (scipy.special.mathieu_b(2 * n, -8 / g2 / g2) - 4)
        else:
            out = (g2 / 8) * (mathieu_b_large_q(2 * n, -8 / g2 / g2) - 4)

        return out * np.sqrt(2)

    return [energy(n + 1) for n in range(k)]


def getGroundState(beta):
    xi = beta / 8.0

    def minFunc(lmbd):
        return eps0(lmbd, xi)

    min_rslt = scipy.optimize.minimize_scalar(minFunc, bracket=(0.1, 1.0))

    lmbd_min = min_rslt.x
    gnd_state = eps0(lmbd_min, xi)

    OldPlaq = 0.5 * P(lmbd_min)
    OldP_exc = 0.5 * Pexc(lmbd_min)

    plaq, p_exc = getPlaquettes(beta)

    if math.isnan(p_exc):
        print("nan detected", lmbd_min)
        plaq, p_exc = OldPlaq, OldP_exc

    e0, e1 = getEnergies(beta)

    # mass_gap = N11(lmbd_min, xi) / D11(lmbd_min, xi) / 4.0 / xi
    mass_gap = e1 - e0

    return e0 + beta, plaq, mass_gap, p_exc


def getLambda(beta):
    xi = beta / 8.0

    def minFunc(lmbd):
        return eps0(lmbd, xi)

    min_rslt = scipy.optimize.minimize_scalar(minFunc, bracket=(0.1, 1.0))

    lmbd_min = min_rslt.x

    return lmbd_min


def getVariationalPrediction(g2, size=(3, 3)):
    """
    """

    def eps0(lmbd, xi):
        return ((((3.0 * lmbd) / (8.0 * xi)) - (4.0 * xi)) * P(lmbd))

    def P(lmbd):
        return 2.0 * scipy.special.iv(2, 4.0 * lmbd) / scipy.special.iv(
            1, 4.0 * lmbd)

    def N11(lmbd):
        return (3 * P(lmbd)) / (2 * lmbd)

    def N12(lmbd):
        return (9 / (4 * lmbd)) * (P(lmbd)**2)

    def N22(lmbd):
        return (9 / lmbd) * P(lmbd) - (27 / (
            (4 * lmbd)**2)) * (P(lmbd)**2) + (9 / (4 * lmbd)) * (P(lmbd)**3)

    def D11(lmbd):
        return (4 - (3 / (2 * lmbd)) * P(lmbd)) - (P(lmbd)**2)

    def D12(lmbd):
        return 8 * P(lmbd) - (3 / lmbd) * (P(lmbd)**2) - 2 * (P(lmbd)**3)

    def D22(lmbd):
        return 8 - (6 / lmbd) * P(lmbd) + (12 + 9 / (8 * (lmbd**2))) * (
            P(lmbd)**2) - (9 /
                           (2 * lmbd)) * (P(lmbd)**3) - (7 / 2) * (P(lmbd)**4)

    xi = 0.5 / g2

    def minFunc(lmbd):
        return eps0(lmbd, xi)

    min_rslt = scipy.optimize.minimize_scalar(minFunc, bracket=(0.1, 1.0))
    lmbd_min = min_rslt.x
    gnd_state = eps0(lmbd_min, xi) * (size[0] - 1) * (size[1] - 1)
    plaq = 0.5 * P(lmbd_min)

    def rootFunc(aM):
        D = np.array([[D11(lmbd_min), D12(lmbd_min)],
                      [D12(lmbd_min), D22(lmbd_min)]])
        N = np.array([[N11(lmbd_min), N12(lmbd_min)],
                      [N12(lmbd_min), N22(lmbd_min)]])
        return np.linalg.det(0.25 * N - xi * aM * D)

    mass_gap = float("nan")
    root_rslt = scipy.optimize.root_scalar(rootFunc,
                                           x0=1,
                                           x1=12,
                                           method="secant")

    mass_gap = root_rslt.root

    return gnd_state + (size[0] - 1) * (size[1] - 1) * 4 / g2, plaq, mass_gap


# print(getVariationalPrediction(0.05, (4, 6)))
# print(np.array(getEnergies(1.0779123358892524)) + 4/1.0779123358892524)
# print(
#     getGroundState(4.0 / 0.1)[0],
#     getGroundState(4.0 / 0.1)[0] +
#     getGroundState(4.0 / 0.1)[2],
#     getGroundState(4.0 / 0.1)[1],
#     getGroundState(4.0 / 0.1)[3])

# print(
#     getGroundState(4.0 / 0.2)[0],
#     getGroundState(4.0 / 0.2)[0] +
#     getGroundState(4.0 / 0.2)[2])


# print(getVariationalPrediction(0.7, (2,3)))

# print(
#     getGroundState(4.0 / 0.09102821015130397)[0],
#     getGroundState(4.0 / 0.09102821015130397)[0] +
#     getGroundState(4.0 / 0.09102821015130397)[2])