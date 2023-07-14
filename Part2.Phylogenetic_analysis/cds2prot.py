from sys import argv

script, infile, outfile1, outfile2 = argv

code = {"UUU":"F","UUC":"F","UUA":"L","UUG":"L","CUU":"L","CUC":"L","CUA":"L","CUG":"L","AUU":"I","AUC":"I","AUA":"I","AUG":"M","GUU":"V","GUC":"V","GUA":"V","GUG":"V","UCU":"S","UCC":"S","UCA":"S","UCG":"S","CCU":"P","CCC":"P","CCA":"P","CCG":"P","ACU":"T","ACC":"T","ACA":"T","ACG":"T","GCU":"A","GCC":"A","GCA":"A","GCG":"A","UAU":"Y","UAC":"Y","UAA":"Stop","UAG":"Stop","CAU":"H","CAC":"H","CAA":"Q","CAG":"Q","AAU":"N","AAC":"N","AAA":"K","AAG":"K","GAU":"D","GAC":"D","GAA":"E","GAG":"E","UGU":"C","UGC":"C","UGA":"Stop","UGG":"W","CGU":"R","CGC":"R","CGA":"R","CGG":"R","AGU":"S","AGC":"S","AGA":"R","AGG":"R","GGU":"G","GGC":"G","GGA":"G","GGG":"G"}

seq = ""
nucl_seq = ""
prot_seq = ""

with open(infile) as fh:
    while True:
        line = fh.readline().strip()
        if len(line) == 0:
            break
        if line.startswith(">"):
            name = line
        else:
            seq += line

count = len(seq)/3
n = 0
for i in range(int(count)):
    codon = seq[n:n+3]
    codon_rna = codon.replace("T","U")
    if "N" in codon_rna:
        prot = "-"
    else:
        prot = code[codon_rna]
    if prot != "Stop":
        prot_seq += prot
        nucl_seq += codon
    else:
        pass
    n += 3

output = open(outfile1, "w")
output.write("%s\n" % name)
output.write("%s\n" % prot_seq)
output.close()

output = open(outfile2, "w")
output.write("%s\n" % name)
output.write("%s\n" % nucl_seq)
output.close()
