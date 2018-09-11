import numpy


def greedymotifsearch(dna, k, t):
    bestmotifs = []

    for i in range(t):
        bestmotifs.append(dna[i][0:k])

    for i in range(len(dna[0]) - k + 1):
        motifs = []
        motifs.append(dna[0][i:i+k])

        for j in range(1, t):
            profile = buildprofile(motifs, k)

            motifs.append(profilemostprobable(profile, dna[j], k, i))

        if score(motifs) < score(bestmotifs):
            bestmotifs = motifs.copy()

    return bestmotifs

def buildprofile(motifs, k):
    profile = numpy.zeros((4, k))

    for r in profile:
        for c in range(len(r)):
            r[c] = 1

    for col in range(k):
        for r in range(len(motifs)):
            if motifs[r][col] == 'A':
                profile[0][col] += 1
            if motifs[r][col] == 'C':
                profile[1][col] += 1
            if motifs[r][col] == 'G':
                profile[2][col] += 1
            if motifs[r][col] == 'T':
                profile[3][col] += 1

    return profile / (len(motifs) + 4)


# return the kmer with the highest probability given the profile
def profilemostprobable(profile, dna, k, index):
    maxprob = -1
    probable = []

    # loop through each kmer in dna
    for i in range(len(dna) - k + 1):
        kmer = dna[i:i+k]
        prob = 1
        j = 0

        for n in kmer:
            if n == 'A':
                prob *= profile[0][j]
            if n == 'C':
                prob *= profile[1][j]
            if n == 'G':
                prob *= profile[2][j]
            if n == 'T':
                prob *= profile[3][j]
            j += 1
        if prob > maxprob:
            maxprob = prob
            probable.append([kmer, maxprob])
            mostprobable = kmer

    return mostprobable


def score(motif):
    score = 0

    consensus = max(set(motif), key=motif.count)

    for i in motif:
        score += hamming_distance(i, consensus)

    return score


def hamming_distance(str1, str2):
    hd = 0

    for i in range(len(str1)):
        if str1[i] != str2[i]:
            hd += 1

    return hd


def readdna(dna, k, t):
    dnalist = []
    n = len(dna)/t
    for i in range(0, len(dna), int(n)):
        dnalist.append(dna[i: i + int(n)][0])
    return dnalist


if __name__ == '__main__':
    # Set values for k and t
    # The file name should be a text file containing a DNA sequence
    filename = ''
    k = 0
    t = 0
    with open(filename) as f:
        dna = f.read().splitlines()

    x = readdna(dna, k, t)

    for i in greedymotifsearch(x, k, t):
        print(i)

