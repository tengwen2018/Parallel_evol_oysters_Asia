from sys import argv

script, infile1, infile2, outfile = argv

# 
seq = ''
codon1 = ''
codon2 = ''
codon3 = ''

with open(infile1) as fh:
    while True:
        line = fh.readline().strip()
        if len(line) == 0:
            break
        if line.startswith(">"):
            name = line.strip(">")
            pass
        else:
            seq += line

count = len(seq)/3
n = 0
for i in range(int(count)):
    codon = seq[n:n+3]
    codon1 += codon[0]
    codon2 += codon[1]
    codon3 += codon[2]
    n += 3

# 
blt_codon1 = ''
blt_codon2 = ''
blt_codon3 = ''
with open(infile2) as fh:
    line = fh.readline()
    while True:
        line = fh.readline().strip().split()
        if len(line) == 0:
            break
        if line[0] == name:
            n = 0
            for i in line[1]:
                if i != "-":
                    blt_codon1 += codon1[n]
                    blt_codon2 += codon2[n]
                    blt_codon3 += codon3[n]
                    n += 1
                else:
                    blt_codon1 += "-"
                    blt_codon2 += "-"
                    blt_codon3 += "-"

output1 = open(outfile+".blt_codon1.fa","w")
output2 = open(outfile+".blt_codon2.fa","w")
output3 = open(outfile+".blt_codon3.fa","w")

output1.write(">%s\n" % outfile)
output1.write("%s\n" % blt_codon1)
output2.write(">%s\n" % outfile)
output2.write("%s\n" % blt_codon2)
output3.write(">%s\n" % outfile)
output3.write("%s\n" % blt_codon3)

output1.close()
output2.close()
output3.close()
