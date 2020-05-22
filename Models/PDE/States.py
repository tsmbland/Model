import numpy as np

"""
Functions for qualitatively defining the state of a model
For use with ParamSpaceQual2D function

"""

"""
PAR model states

"""


def par_state_pol(a, p):
    if sum(a > p) == len(a):
        # A dominant
        return 1
    elif sum(a > p) == 0:
        # P dominant
        return 1
    else:
        # Polarised
        return 2


def par_state_pol2(a, p):
    if sum(a > p) == len(a):
        # A dominant
        return 3
    elif sum(a > p) == 0:
        # P dominant
        return 1
    else:
        # Polarised
        return 2


def par_state_bp(a, p):
    if sum(a > p) == len(a):
        # Unpolarised, A dominant
        return 1
    elif sum(a > p) == 0:
        # Unpolarised, P dominant
        return 1
    else:
        if sum(a > p) < len(a) // 2:
            # Polarised, boundary anterior of centre
            return 2
        elif sum(a > p) > len(a) // 2:
            # Polarised, boundary posterior of centre
            return 3
        elif abs(a[len(a) // 2] - p[len(a) // 2]) > abs(a[(len(a) // 2) - 1] - p[(len(a) // 2) - 1]):
            # Polarised, boundary just anterior of centre
            return 2
        else:
            # Polarised, boundary just posterior of centre
            return 3


def par_state_asi_p(a, p):
    if sum(a > p) == len(a):
        # A dominant
        return 1
    elif sum(a > p) == 0:
        # P dominant
        return 1
    else:
        # Polarised
        ant = np.mean(p[:50])
        post = np.mean(p[50:])
        asi = abs((ant - post) / (2 * (ant + post)))
        if asi < 0.2:
            return 2
        elif asi < 0.35:
            return 3
        elif asi < 0.45:
            return 4
        elif asi < 0.49:
            return 5
        else:
            return 6


def par_state_asi_a(a, p):
    if sum(a > p) == len(a):
        # A dominant
        return 1
    elif sum(a > p) == 0:
        # P dominant
        return 1
    else:
        # Polarised
        ant = np.mean(a[:50])
        post = np.mean(a[50:])
        asi = abs((ant - post) / (2 * (ant + post)))
        if asi < 0.2:
            return 2
        elif asi < 0.35:
            return 3
        elif asi < 0.45:
            return 4
        elif asi < 0.49:
            return 5
        else:
            return 6


def par_state_asi_a2(a, p):
    # Polarised
    ant = np.mean(a[:50])
    post = np.mean(a[50:])
    asi = abs((ant - post) / (2 * (ant + post)))
    if asi < 0.2:
        return 2
    elif asi < 0.35:
        return 3
    elif asi < 0.45:
        return 4
    elif asi < 0.49:
        return 5
    else:
        return 6

# def par_state_dc(a, p, thresh):
#     if sum(a > p) == len(a):
#         # A dominant
#         return 4
#     elif sum(a > p) == 0:
#         # P dominant
#         return 1
#     else:
#         # Polarised
#         if p[-1] > thresh:
#             # Concentration too high
#             return 2
#         else:
#             # Concentrated too low
#             return 3
#
#
# def par_state_dc2(a, p):
#     if sum(a > p) == len(a):
#         # A dominant
#         return 1
#     elif sum(a > p) == 0:
#         # P dominant
#         return 1
#     else:
#         # Polarised
#         if p[-1] < 2:
#             return 2
#         elif p[-1] < 5:
#             return 3
#         elif p[-1] < 10:
#             return 4
#         elif p[-1] < 15:
#             return 5
#         else:
#             return 6
