# coding:        utf-8
# filename:      vdw.py
# author:        Henry Mitchell
# creation date: 26 Mar 2017

# description:   Find the van der Waals forces between
#                a sheet of a material and an atom.

import argparse
import numpy as np
import scipy.integrate as integrate

from pympler import tracker

tr = tracker.SummaryTracker()

π = np.pi


def c3(α0, ω0, g, z, Δ, v):
    '''Calculates c3'''
    Ω = 2 * ω0 * z / v
    δ = 2 * z * Δ / v
    return integrate.nquad(kernel, [[0, np.inf], [0, 2 * π], [0, np.inf]],
                           args=(g, Ω, δ))[0] * α0 * ω0 / (8 * π)


def kernel(q, φ, ω, g, Ω, δ):
    '''The kernel of the integral used to find c3'''
    fin = np.abs(vpi(g, q, ω, Ω, δ)) / (1 - vpi(g, q, ω, Ω, δ))
    return q**2 * np.exp(-q) * (1 / (1 + ω**2)) * (1 / (2 * π)) * fin


def vpi(g, q, ω, Ω, δ):
    '''Calculates the potential times the polarizability'''
    w = q**2 + (ω**2) * (Ω**2)
    return -2 * g * q * ((δ / w) + 1 / (2 * np.sqrt(w)) * (1 - 4 * δ**2 / w)
                         * np.arctan(np.sqrt(w) / (2 * δ)))


if __name__ == "__main__":
    elements = {"H": {"α0": 4.5 * 1.4818e-4, "ω0": 11.65},
                "Na": {"α0": 162.6 * 1.4818e-4, "ω0": 2.15},
                "He": {"α0": 1.38 * 1.4818e-4, "ω0": 27}}
    materials = {"graphene": {"v": 0.658, "Δ": 1e-100, "g": 2.2},
                 "MoS2": {"v": 0.351, "Δ": 0.835, "g": 4.1}}
    parser = argparse.ArgumentParser(description="Create a list of C3s")
    parser.add_argument("atom", metavar="X", type=str, nargs=1,
                        help="the atom which we're analyzing")
    parser.add_argument("material", metavar="S", type=str, nargs=1,
                        help="MoS2 or graphene")
    parser.add_argument("-g", action="store_true",
                        help="run through the values of g",
                        required=False)
    args = parser.parse_args()

    zs = np.arange(2, 140, 2)
    element = elements[args.atom[0]]
    material = materials[args.material[0]]
    α0 = element["α0"]          # nm^3
    ω0 = element["ω0"]          # eV
    γ = [material["g"]]         # dimensionless
    v = material["v"]           # eV nm
    Δ = material["Δ"]           # eV

    if args.g:
        gs = [0.5]
        gs = np.append(gs,
                       [g for g in np.arange(1, np.floor(material["g"]), 1)])
        gs = np.append(gs, [γ])
    else:
        gs = γ

    for g in gs:
        for z in zs:
            c = c3(α0, ω0, g, z, Δ, v) / 4.032e-3  # a.u. of C3
            tr.print_diff()
            print(z)
            with open("{}-{}-{}.txt".format(args.atom[0],
                                            args.material[0],
                                            g), "ab") as data_file:
                np.savetxt(data_file, np.array([[z, c]]))
