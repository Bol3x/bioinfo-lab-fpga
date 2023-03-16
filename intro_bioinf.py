#program to calculate frame alignment using Petola, Soderlund, and Ukkonen's algorithm
#source: Algorithms for the search of amino acid patterns in nucleic acid sequences. Nucleaic Acids Research, Volume 14, Number 18, 1986.

INF = 1e9

del1 = 2
del_1 = 2
del2 = 4
del_2 = 4


#convert codon to amino acid
#time complexity: O(n) - https://wiki.python.org/moin/TimeComplexity#dict
def codonToAmino(codon):
    codon_table = {
		'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
		'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
		'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
		'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
		'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
		'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
		'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
		'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
		'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
		'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
		'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
		'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
		'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
		'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
		'TAC':'Y', 'TAT':'Y', 'TAA':'-', 'TAG':'-',
		'TGC':'C', 'TGT':'C', 'TGA':'-', 'TGG':'W',
	}
    return codon_table[codon]

# a is the amino acid sequence
a = input()

# b is the nucleotide sequence
b = input()

# n is the length of the amino acid sequence
n = len(a)

# m is the length of the nucleotide sequence
m = len(b)

def calculatePenalty(a, amino1, amino2):
    return 0 if amino1 == amino2 or (amino1 == "X" and amino2 not in a) else 3

def originalFrameAlgo(a, b, n, m):
	#initialize memo
	memo = [[INF for i in range(m+1)] for j in range(n+1)]
 
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
				penalty = calculatePenalty(a, a[i], codonToAmino(b[j-3:j]))
				#calculate memo
				memo[i][j] = penalty + min(
								memo[i-1][j-1]+del_2,
								memo[i-1][j-2]+del_1,
								memo[i-1][j-3],
								memo[i-1][j-4]+del1,
								memo[i-1][j-5]+del2,
							)
    
	return memo[n-1][m-1]