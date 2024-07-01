import dnaio
import sys
import pandas as pd

FQfile = sys.argv[1]
samfile = sys.argv[2]

seqdict = dict()

with dnaio.open(FQfile, mode='r') as readers:
    for record in readers:
        rn = (record.name.split(" ")[0])
        seqdict[rn] = str(record.sequence)

oname = str(FQfile.split(".")[0])+ "_umi_bwtMap.csv"
out = open(oname,"w+")
out.write("UMI_Template,Name,SAM_flag,Subject,Pos,MapQ,CIGAR,mate,QueryStart,Query\n")
with open(samfile, "r") as infil:
    inf = infil.readlines()
    for i in inf:
        i = i.strip()
        i = i.split("\t")
        j = seqdict[i[0]]
        z = [j] + i[:10]
        n = ",".join(z)
        out.write(n+ "\n")
        #print(n)

df = pd.read_csv(oname,index_col=False)
print(df.columns)
print(df.head())
new = df.groupby(['UMI_Template','SAM_flag','Subject'])['Query'].agg('count').reset_index()
onamecounts = str(FQfile.split(".")[0])+ "_umi_bwtMapCounts.csv"
#print(new)
new.to_csv(onamecounts, index=False)

#SRR5234336.17	0	hsa-miR-490-3p	1	255	22M	*	0	0	CAACCTGGAGGACTCCATGCTG	BFFFFFFFFFFFFFFFFFFFFF	XA:i:0	MD:Z:22	NM:i:0	XM:i:2
#SRR5234336.18	0	hsa-miR-181a-2-5p	2	255	22M	*	0	0	ACATTCAACGCTGTCGGTGAGT	BFFFFFFFFFFFFFFFFFFFFF	XA:i:0	MD:Z:22	NM:i:0	XM:i:2
#SRR5234336.19	0	hsa-miR-500a-5p	1	255	23M	*	0	0	TAATCCTTGCTACCTGGGTGAGA	BFBFFFFFFFFFFFFFFFFFFF<	XA:i:0	MD:Z:23	NM:i:0	XM:i:2
#SRR5234336.20	0	hsa-miR-196b-5p	1	255	22M	*	0	0	TAGGTAGTTTCCTGTTGTTGGG	BFFFFFFFFFFFFFFFFFFFFF	XA:i:0	MD:Z:22	NM:i:0	XM:i:2
#AACACAACCTGGAGGACTCCATGCTGGTCC  SRR5234336.17   0       hsa-miR-490-3p  1       255     22M     *       0       0       CAACCTGGAGGACTCCATGCTG
#NTCAACATTCAACGCTGTCGGTGAGTTCTT  SRR5234336.18   0       hsa-miR-181a-2-5p       2       255     22M     *       0       0       ACATTCAACGCTGTCGGTGAGT
#NAGGTAATCCTTGCTACCTGGGTGAGATGTA SRR5234336.19   0       hsa-miR-500a-5p 1       255     23M     *       0       0       TAATCCTTGCTACCTGGGTGAGA
