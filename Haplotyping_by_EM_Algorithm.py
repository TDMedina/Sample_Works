# -*- coding: utf-8 -*-
"""
Created on Wed Mar 28 17:07:58 2018

@author: Tyler
"""

import itertools as itt

# Genotypes for haplotyping. "0" and "1" represent homozygous loci, "2"
# represents heterozygous loci.
genos = ["0112", "1022", "0122", "1122", "1022"]
# genos = ["112", "111", "201", "021"]

# Empty list to add haplotypes to.
possible_haplotypes = []

# Forms every possible haplotype based on the genotypes.
for geno in genos:
    geno = list(geno)
    for locus in range(len(geno)):
        if geno[locus] == "1":
            geno[locus] = ["0", "1"]
        elif geno[locus] == "2":
            geno[locus] = "1"
    permutation = list(itt.product(*geno))
    permutation = ["".join(haplo) for haplo in permutation]
    possible_haplotypes.append(permutation)

# Flattens the list, removes duplicates, and sorts.
possible_haplotypes = [x for l in possible_haplotypes for x in l]
unique_haplos = list(set(possible_haplotypes))
unique_haplos.sort()

# Produces every combination of 2 haplotypes.
combos = list(itt.combinations_with_replacement(unique_haplos, 2))

# Creates genotypes from the 'sum' of each haplotype combination.
hap_sums = []
for pair in combos:
    hap_sum = int(pair[0]) + int(pair[1])
    hap_sums.append(hap_sum)

# List of haplotype names e.g. "h1", "h2", and a list of names corresponding to
# all combinations of the haplotypes.
names = ["h" + str(x + 1) for x in range(len(unique_haplos))]
name_combos = list(itt.combinations_with_replacement(names, 2))

# Named list of haplotypes.
haplo_dic = dict(zip(names, unique_haplos))

# Zipped list of combo names and combo genotypes.
pair_ids = list(zip(name_combos, hap_sums))

# Creates a list of the haplotype-pairs that are possible for each genotype.
possible_pairs = ["x"] * len(genos)
k = 0
for geno in genos:
    geno_matches = []
    for pair in pair_ids:
        if pair[1] == int(geno):
            geno_matches.append(pair[0])
    possible_pairs[k] = list(geno_matches)
    k += 1

# Creates a dictionary of haplotype-name: population-frequency. Initialized
# with equal frequencies for each haplotype.
pop_freqs = dict(zip(names, [1/(len(unique_haplos))] * len(unique_haplos)))


# Function to perform one EM algorithm iteration for haplotype frequencies.
def Haplotype_EM():

    # List to hold the haplotype-pair frequency per genotype.
    freqs_by_geno = []

    # Calculates the frequency of each haplotype-pair for each genotype.
    for geno in possible_pairs:
        pair_prods = []
        for pair in geno:
            pair_prod = pop_freqs[pair[0]] * pop_freqs[pair[1]]
            pair_prods.append(pair_prod)
        freqs = [x/sum(pair_prods) for x in pair_prods]
        freqs = list(zip(geno, freqs))
        freqs_by_geno.append(freqs)

    # Calculates the total weighted population frequency of each haplotype, and
    # overwrites the current population frequency.
    for haplotype in names:
        hap_freq = []
        for geno_freq in freqs_by_geno:
            for pair in geno_freq:
                weighted_occurence = pair[0].count(haplotype) * pair[1]
                hap_freq.append(weighted_occurence)
        pop_freqs[haplotype] = (sum(hap_freq)) / (2 * len(genos))
    return(pop_freqs)


# Prints out naive initial haplotype frequencies.
print("Initial Haplotype Frequencies:")
for value in pop_freqs:
    print("  {}: {:.3}".format(value, pop_freqs[value]))
print("\n")

# Loop to perform EM iterations. Stops once the population frequency converges
# (rounded to 6 decimal places). Does this by storing the current population
# frequencies in a temporary list, then performs an EM iteration, then
# compares the two. If there is no change, the algorithm terminates and prints.
t = 0
end_it = 0
while end_it < 1:
    temp = list(pop_freqs.values())
    temp = [round(x, 6) for x in temp]
    first_it = list(Haplotype_EM().values())
    first_it = [round(x, 6) for x in first_it]
    t += 1
    if first_it == temp:
        end_it += 1
        print("EM Algorithm Haplotype Frequency Results:")
        for value in pop_freqs:
            print("  {}: {:.3f}".format(value, pop_freqs[value]))
        print("\n")
        print("Iterations until convergence:\n  {}\n".format(t))

# Recalculates the final haplotype-pair frequencies per genotype using the
# final population haplotype frequencies.
freqs_by_geno = []
for geno in possible_pairs:
    pair_prods = []
    for pair in geno:
        pair_prod = pop_freqs[pair[0]] * pop_freqs[pair[1]]
        pair_prods.append(pair_prod)
    freqs = [x/sum(pair_prods) for x in pair_prods]
    freqs = list(zip(geno, freqs))
    freqs_by_geno.append(freqs)

# Prints the haplotype pair with the highest probability for each genotype, as
# well as its probability.
for geno, pair in zip(genos, freqs_by_geno):
    best = max(pair)
    print("Genotype: {}".format(geno))
    print("  Most likely haplotype pair:")
    print("    {}: {}\n    {}: {}".format(best[0][0], haplo_dic[best[0][0]],
                                          best[0][1], haplo_dic[best[0][1]]))
    print("  Probability: {:.3}\n".format(best[1]))
