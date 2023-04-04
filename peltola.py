#program to calculate frame alignment using Petola, Soderlund, and Ukkonen's algorithm
#source: Algorithms for the search of amino acid patterns in nucleic acid sequences. Nucleaic Acids Research, Volume 14, Number 18, 1986.

from helpers import codonToAmino
from helpers import aminoToCodon

INF = 1e9

del1 = 2
del_1 = 2
del2 = 4
del_2 = 4

# a is the amino acid sequence
a = input()

# b is the nucleotide sequence
b = input()

# n is the length of the amino acid sequence
n = len(a)

# m is the length of the nucleotide sequence
m = len(b)

#calculates the penalty of matching amino1 to amino2
def calculatePenalty(a, amino1, amino2):
    return 0 	if (amino1 == amino2) or (amino1 == "X" and amino2 not in a) 	else 3
	

#main peltola algorithm
def peltolaAlg(amino_seq, nucleotide_seq, n, m):
	#initialize memo
	memo = [[INF for _ in range(m+1)] for _ in range(n+1)]
 
	for i in range(n):
		for j in range(-2,m,1):
    
			#base cases
			if i == 0 and j in range():
				memo[i][j] = INF
		
			if i == 0 and j > 0:
				memo[i][j] = 0
    
			if i > 0 and j in range(-2,3):
				memo[i][j] = INF
    
			else:
				#calculate penalty
				penalty = calculatePenalty(amino_seq, amino_seq[i], codonToAmino(nucleotide_seq[j-2:]))
				#calculate memo
				memo[i][j] = penalty + min(
								memo[i-1][j-1] + del_2,
								memo[i-1][j-2] + del_1,
								memo[i-1][j-3],
								memo[i-1][j-4] + del1,
								memo[i-1][j-5] + del2,
							)
    
	memo = peltolaCont(memo, n, m, amino_seq, nucleotide_seq)
	return memo[n-1][m-1]




#finds the minimum edit distance between a single amino acid and a codon within a nucleotide sequence seq
def minCodonPenalty(amino, seq):
	minPenalty = 0
 
 	#todo: create amino to codon function, returns a list of codons that encode to the amino
	amino_codons = aminoToCodon(amino)

	for i in range(len(seq)):
		#get codon from seq
		seq_codon = seq[i-2:]

		#todo: implement/call an edit distance implementation 
  		# to find minimum edit distance between amino and codon, for all shifts in seq
		for codon in amino_codons:
			penalty = editDistance(codon, seq_codon)
		
			if penalty < minPenalty:
				minPenalty = penalty
   
	return minPenalty




# continuation of main algorithm, used to find more homologies in the sequence
def peltolaCont(memo, n, m, amino_seq, nucleotide_seq):
    for i in range(n):
        for j in range(m):
            #update matrix E
            memo[i][j] = min(
				memo[i-1][j-1] + minCodonPenalty(amino_seq[i], nucleotide_seq[j]),
				memo[i-1][j-2] + minCodonPenalty(amino_seq[i], nucleotide_seq[j-1:]),
				memo[i-1][j-3] + minCodonPenalty(amino_seq[i], nucleotide_seq[j-2:]),
				memo[i-1][j-4] + minCodonPenalty(amino_seq[i], nucleotide_seq[j-3:]),
				memo[i-1][j-5] + minCodonPenalty(amino_seq[i], nucleotide_seq[j-4:])
			)
            
    return memo
    