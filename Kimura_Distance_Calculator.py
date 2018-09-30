# -*- coding: utf-8 -*-
"""
@author: Ty Medina
"""

# Import necessary libraries.
import numpy as np
from scipy.optimize import minimize
from scipy.linalg import expm

# Set aligned sequences.
seq1 = "agtccatgat"
seq2 = "acgtcgtgct"


# Function to calculate probability of distance.
#   Note: "sign" is used to find negative "maximum" by using
#          scipy.optimize.minimize.
def Kimura_Distance(tkappa, sign=1):
    t = tkappa[0]
    k = tkappa[1]

# Sets the unscaled generator matrix.
    generator = np.matrix([[0, 1, 1, k],
                           [1, 0, k, 1],
                           [1, k, 0, 1],
                           [k, 1, 1, 0]])

# Calculates scale factor; note that stationary distribution is 1/4 and is thus
# ignored, since it will cancel out anyway.
    scale_factor = np.sum(generator[0])
    print(generator[0])
    print(scale_factor)
# Sets diagonal indices.
    diagonal = np.diag_indices(4)
# Scales the generator matrix.
    generator = (generator/scale_factor)
    #print(generator)
# Sets the diagonal of the scaled matrix so that sum of each row is zero; note
# that it does not matter when the diagonals are set up till now, as they will
# always end up as -1 after scaling.
    generator[diagonal] = -1

# Produces the transitional probability matrix, depending on "t"
    tpm = expm(generator * t)
# Sets up a dictionary to convert between the nucleotides in the sequence and
# the indices of the matrix.
    tpm_dex = dict(zip(["a", "t", "c", "g"], [0, 1, 2, 3]))

# Dummy variable to start the ongoing product in the loop below.
    product = 1
# Loop to reference the TPM locations for the nucleotides at each position,
# then multiply their values together.
    for i in range(len(seq1)):
        product = product * tpm[tpm_dex[seq1[i]], tpm_dex[seq2[i]]]

# As mentioned, "sign" allows the minimization function to find the maximum by
# by reporting the greatest negative value.
    answer = product * sign
    return(answer)


# Optimizer to find the maximum value (minimum value * -1).
# !!!!VERY VERY IMPORTANT!!!! Default method is not suitable for this function!
# The function is not smooth. As a result, the optimizer, which depends on
# derivative values, will get stuck on the first iteration and report the
# initial value only. Nelder-Mead method is not immune to this, but search
# window happens to be wide enough to bridge to next non-smooth step, allowing
# iteration.
optimization = minimize(Kimura_Distance, x0=[1, 0],
                        args=-1, method="Nelder-Mead")

# ==Display results==================================================
print("\n===Estimating Optimal Kimura-Model Evolutionary Distance===\n")
print("Aligned Sequences:\nSequence 1: " +
      "{}\nSequence 2: {}\n".format(seq1.upper(), seq2.upper()))
print("---Results of Nelder-Mead Optimization---\n")
print("Optimal t: {:.5}".format(optimization["x"][0]))
print("Optimal k: {:.6}".format(optimization["x"][1]))
print("Max Probability: {:.5}".format(-1 * optimization["fun"]))

# ==Distance from formula=====
p = 2/10
q = 3/10
d = (-0.5)*(np.log((1 - 2*p - q)*np.sqrt(1 - 2*q)))
print("\n\nTest Distance (t) from formula: {:.5}".format(d))
