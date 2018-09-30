# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 10:18:30 2018

@author: Tyler Medina
Forwards-Backwards Algorithm
"""

import numpy as np

# ===MANUALLY-ENTERED DATA=====================================================

# TPM of the hidden states
tprob = np.matrix([[0.5, 0.5], [0.2, 0.8]])
# List of the names of the hidden states
hid = ["Active", "Dormant"]
# List of the observed symbols
symbols = ["Low", "Med", "Norm"]
# Emission probabilities
eprob1 = [0.7, 0.2, 0.1]
eprob2 = [0.2, 0.3, 0.5]

# Structured array holding tuples with symbol, emission probability for state 1
# and emission probability for state 2
obs = np.array(list(zip(symbols, eprob1, eprob2)),
               dtype=[('symbol', 'U10'), ('eprobs1', 'f4'), ('eprobs2', 'f4')])

# Observed symbol sequence
seq = ["Low", "Norm", "Med", "Low", "Norm", "Norm", "Norm", "Med"]
# =============================================================================


# ===STATIONARY DISTRIBUTION===================================================

# Array, where [0] is the eigenvalues and [1] is the eigenvectors matrix
vector = np.linalg.eig(np.transpose(tprob))
# The index of eigenvalue == 1
eigen_loc = np.argmin((abs(vector[0]-1)))
# The row vector of eigenvalue == 1
raw_pi = np.asarray(np.transpose(vector[1][:, eigen_loc]))
# The stationary distribution vector
pi_vector = np.array([x/np.sum(raw_pi) for x in np.nditer(raw_pi)])
# =============================================================================


# ===FORWARDS ALGORITHM========================================================

# Initializes the forward algorithm by storing the start value for each
# hidden state, which is emission prob times stationary dist.
initial = obs[int(np.where(obs['symbol'] == seq[0])[0])]
forward = [[pi_vector[0]*initial[1]], [pi_vector[1]*initial[2]]]

i = 1
while i < len(seq):
    observed_x = obs[int(np.where(obs['symbol'] == seq[i])[0])]
    for state in range(len(hid)):
        prev_sum = 0.0
        for prev in range(len(hid)):
            prev_sum += forward[prev][i-1]*tprob[prev, state]
        forward[state].append(prev_sum*observed_x[state+1])
    i += 1
# =============================================================================


# ===BACKWARDS ALGORITHM=======================================================

# Initializes the backwards algorithm by storing a dummy 0 for each position
# and a 1 for the last entry for each hidden state.
back = [([0]*(len(seq)-1)) + [1], ([0]*(len(seq)-1)) + [1]]

j = len(seq) - 2
while j >= 0:
    obs_after = obs[int(np.where(obs['symbol'] == seq[j+1])[0])]
    for state in range(len(hid)):
        aftersum = 0.0
        for after in range(len(hid)):
            aftersum += back[after][j+1]*tprob[state, after]*obs_after[after+1]
        back[state][j] = aftersum
    j -= 1
# =============================================================================


# ===P(x)======================================================================
ends = []
for hidden in forward:
    ends.append(hidden[-1])

px = sum(ends)
# =============================================================================


# ===ASSIGNMENT ANSWER=========================================================
print("MA461 Assignment 2: Forwards-Backwards Algorithm\n")
print("Hidden states:")
for name in hid:
    print(name)
print("\nTransition Probability Matrix:\n{}\n".format(tprob))
print("Observed symbols:\n{}\n".format(symbols))
print("Emission probabilities:")
print("Active: {}\nDormant: {}\n".format(eprob1, eprob2))
print('Observed sequence:\n{}\n'.format(seq))
print("Probability of being in the \'active\' state at month 6:")
print((forward[0][5]*back[0][5])/(px))

print(forward)
print(back)
print(px)
# =============================================================================
