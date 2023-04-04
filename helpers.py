# converts a codon (3 nucleotide) sequence into its respective amino acid
# stop codons are encoded with '-'
def codonToAmino(codon):
    codon_table = {
		'UUU': 'F', 'UUC': 'L', 'UUA': 'L', 'UUG': 'L',
		'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
		'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
		'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
		'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
		'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
		'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
		'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
		'UAU': 'Y', 'UAC': 'Y', 'UAA': '-', 'UAG': '-',
		'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
		'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
		'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
		'UGU': 'C', 'UGC': 'C', 'UGA': '-', 'UGG': 'W',
		'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
		'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
		'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
	}
	#return X if codon is not in the table
    return codon_table.get(codon, 'X')

# PAM250 Dictionary Matrix, taken from https://gist.github.com/skyfox/5123558
PAM250 = {
'A': {'A': 2,  'C': -2, 'D':  0, 'E': 0, 'F': -3, 'G':  1, 'H': -1, 'I': -1, 'K': -1, 'L': -2, 'M': -1, 'N':  0, 'P':  1, 'Q':  0, 'R': -2, 'S':  1, 'T':  1, 'V':  0, 'W': -6, 'Y': -3},
'C': {'A': -2, 'C': 12, 'D': -5, 'E':-5, 'F': -4, 'G': -3, 'H': -3, 'I': -2, 'K': -5, 'L': -6, 'M': -5, 'N': -4, 'P': -3, 'Q': -5, 'R': -4, 'S':  0, 'T': -2, 'V': -2, 'W': -8, 'Y':  0},
'D': {'A': 0,  'C': -5, 'D':  4, 'E': 3, 'F': -6, 'G':  1, 'H':  1, 'I': -2, 'K':  0, 'L': -4, 'M': -3, 'N':  2, 'P': -1, 'Q':  2, 'R': -1, 'S':  0, 'T':  0, 'V': -2, 'W': -7, 'Y': -4},
'E': {'A': 0,  'C': -5, 'D':  3, 'E': 4, 'F': -5, 'G':  0, 'H':  1, 'I': -2, 'K':  0, 'L': -3, 'M': -2, 'N':  1, 'P': -1, 'Q':  2, 'R': -1, 'S':  0, 'T':  0, 'V': -2, 'W': -7, 'Y': -4},
'F': {'A': -3, 'C': -4, 'D': -6, 'E':-5, 'F':  9, 'G': -5, 'H': -2, 'I':  1, 'K': -5, 'L':  2, 'M':  0, 'N': -3, 'P': -5, 'Q': -5, 'R': -4, 'S': -3, 'T': -3, 'V': -1, 'W':  0, 'Y':  7},
'G': {'A': 1,  'C': -3, 'D':  1, 'E': 0, 'F': -5, 'G':  5, 'H': -2, 'I': -3, 'K': -2, 'L': -4, 'M': -3, 'N':  0, 'P':  0, 'Q': -1, 'R': -3, 'S':  1, 'T':  0, 'V': -1, 'W': -7, 'Y': -5},
'H': {'A': -1, 'C': -3, 'D':  1, 'E': 1, 'F': -2, 'G': -2, 'H':  6, 'I': -2, 'K':  0, 'L': -2, 'M': -2, 'N':  2, 'P':  0, 'Q':  3, 'R':  2, 'S': -1, 'T': -1, 'V': -2, 'W': -3, 'Y':  0},
'I': {'A': -1, 'C': -2, 'D': -2, 'E':-2, 'F':  1, 'G': -3, 'H': -2, 'I':  5, 'K': -2, 'L':  2, 'M':  2, 'N': -2, 'P': -2, 'Q': -2, 'R': -2, 'S': -1, 'T':  0, 'V':  4, 'W': -5, 'Y': -1},
'K': {'A': -1, 'C': -5, 'D':  0, 'E': 0, 'F': -5, 'G': -2, 'H':  0, 'I': -2, 'K':  5, 'L': -3, 'M':  0, 'N':  1, 'P': -1, 'Q':  1, 'R':  3, 'S':  0, 'T':  0, 'V': -2, 'W': -3, 'Y': -4},
'L': {'A': -2, 'C': -6, 'D': -4, 'E':-3, 'F':  2, 'G': -4, 'H': -2, 'I':  2, 'K': -3, 'L':  6, 'M':  4, 'N': -3, 'P': -3, 'Q': -2, 'R': -3, 'S': -3, 'T': -2, 'V':  2, 'W': -2, 'Y': -1},
'M': {'A': -1, 'C': -5, 'D': -3, 'E':-2, 'F':  0, 'G': -3, 'H': -2, 'I':  2, 'K':  0, 'L':  4, 'M':  6, 'N': -2, 'P': -2, 'Q': -1, 'R':  0, 'S': -2, 'T': -1, 'V':  2, 'W': -4, 'Y': -2},
'N': {'A': 0,  'C': -4, 'D':  2, 'E': 1, 'F': -3, 'G':  0, 'H':  2, 'I': -2, 'K':  1, 'L': -3, 'M': -2, 'N':  2, 'P':  0, 'Q':  1, 'R':  0, 'S':  1, 'T':  0, 'V': -2, 'W': -4, 'Y': -2},
'P': {'A': 1,  'C': -3, 'D': -1, 'E':-1, 'F': -5, 'G':  0, 'H':  0, 'I': -2, 'K': -1, 'L': -3, 'M': -2, 'N':  0, 'P':  6, 'Q':  0, 'R':  0, 'S':  1, 'T':  0, 'V': -1, 'W': -6, 'Y': -5},
'Q': {'A': 0,  'C': -5, 'D':  2, 'E': 2, 'F': -5, 'G': -1, 'H':  3, 'I': -2, 'K':  1, 'L': -2, 'M': -1, 'N':  1, 'P':  0, 'Q':  4, 'R':  1, 'S': -1, 'T': -1, 'V': -2, 'W': -5, 'Y': -4},
'R': {'A': -2, 'C': -4, 'D': -1, 'E':-1, 'F': -4, 'G': -3, 'H':  2, 'I': -2, 'K':  3, 'L': -3, 'M':  0, 'N':  0, 'P':  0, 'Q':  1, 'R':  6, 'S':  0, 'T': -1, 'V': -2, 'W':  2, 'Y': -4},
'S': {'A': 1,  'C':  0, 'D':  0, 'E': 0, 'F': -3, 'G':  1, 'H': -1, 'I': -1, 'K':  0, 'L': -3, 'M': -2, 'N':  1, 'P':  1, 'Q': -1, 'R':  0, 'S':  2, 'T':  1, 'V': -1, 'W': -2, 'Y': -3},
'T': {'A': 1,  'C': -2, 'D':  0, 'E': 0, 'F': -3, 'G':  0, 'H': -1, 'I':  0, 'K':  0, 'L': -2, 'M': -1, 'N':  0, 'P':  0, 'Q': -1, 'R': -1, 'S':  1, 'T':  3, 'V':  0, 'W': -5, 'Y': -3},
'V': {'A': 0,  'C': -2, 'D': -2, 'E':-2, 'F': -1, 'G': -1, 'H': -2, 'I':  4, 'K': -2, 'L':  2, 'M':  2, 'N': -2, 'P': -1, 'Q': -2, 'R': -2, 'S': -1, 'T':  0, 'V':  4, 'W': -6, 'Y': -2 }}

#Make sure to `pip install blosum` before running this code
import blosum as bl

BLOSUM62 = bl.BLOSUM(62)

# calculates the score for a given amino acid pair
# returns their score based on the matrix used
# can use either PAM250 or BLOSUM62 scoring matrices, selected via scoringType parameter
def getScore(a, b, scoringType):
    
    if scoringType == "PAM250":
        #todo: return score from PAM250 "matrix"
        return PAM250[a][b]
    
    if scoringType == "BLOSUM62":
        #todo: fill in BLOSUM62 scoring matrix
        return BLOSUM62[a][b]
    

    