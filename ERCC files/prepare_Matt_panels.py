import sys
inf = sys.argv[1]
infOut  = str(inf.split("_")[0]) + "_UMI_panel_A.txt"
infOut2 = str(inf.split("_")[0]) + "_UMI_panel_B.txt"

out = open(infOut, "w+") 
out.write("Template\tUMI1\tUMI2\tmiRNA\tCount\n")
panelb = dict()
with open(inf,"r") as infile:
    inf = infile.readlines()
    for i in inf:
        i = i.strip()
        i = i.split(",")
        if i[0] != "UMI_Template":
            #print(i)
            #print(i[0]+"\t"+i[0][4:-4]+"\t"+str(i[0][:4])+"\t"+str(i[0][-4:])+"\t"+i[2]+"\t"+i[3]+"\n")
            req_seq = i[0][4:-4]+"#"+ i[2]
            try:
                panelb[req_seq].append(str(i[3]))
            except KeyError:
                panelb[req_seq] = [str(i[3])]
            out.write(i[0][4:-4]+"\t"+str(i[0][:4])+"\t"+str(i[0][-4:])+"\t"+i[2]+"\t"+i[3]+"\n")


#out2 = 
out2 = open(infOut2, "w+") 
out2.write("Template\tmiRNA\tTotalCount\tUMI_counts\n")
for x, y in panelb.items():
    #print(x, y)
    #print(str(sum(int(i) for i in y)) +"\t" + str(len(y))+"\n")
    if x.split("#")[1] == "*":
        out2.write(x.split("#")[0] + "\t-\t" + str(sum(int(i) for i in y)) +"\t" + str(len(y))+"\n")
    else:
        out2.write(x.split("#")[0] + "\t" + x.split("#")[1] + "\t" + str(sum(int(i) for i in y)) +"\t" + str(len(y))+"\n")
    #print(str(x.split("#")[0])+"\t"+str(x.split("#")[1])+"\t"+str(sum(y)) +"\t" + str(len(y))+"\n")


#reads_templates_umi_bwtMapCounts.csv
#['UMI_Template', 'SAM_flag', 'Subject', 'Query\n']
#['AAAAAAAAAAAAAAGTGC', '0', 'HBA2-15', '1\n']
#['AAAAAAAAAAAAAATTGC', '0', 'HBA2-15', '2\n']
#['AAAAAAAAAAAAAGCTG', '0', 'HBA2-15', '1\n']
#['AAAAAAAAAAAAAGTGA', '0', 'HBA2-15', '1\n']
#['AAAAAAAAAACGTGAGGCGCTGCTATCGTT', '4', '*', '1\n']
#['AAAAAAAAAATCCAACAATAAAATAACCATGGAAATAAAAAAAAAAAATCAA', '4', '*', '1\n']
#['AAAAAAAAAATCCCCTCACAAACTATCACCTGAAAGCCGGGTAAACAGTCTC', '4', '*', '1\n']
#['AAAAAAAAACAGGTGCAATAGTATTGTCAAAGCTCTC', '4', '*', '1\n']
#['AAAAAAAAAGTTAAAA', '4', '*', '1\n']
#['AAAAAAAAATTAAAATAGACAAATAAGATGGAATAGGAGTAAAAAATGTTTA', '4', '*', '1\n']
#['AAAAAAAACATCCTCGACTGGAAGGCAC', '4', '*', '1\n']
#['AAAAAAAACATCCTCGACTGGAAGGTGT', '4', '*', '2\n']
#['AAAAAAAACATTCGCGGTGCACTTCTTCGAA', '4', '*', '1\n']
#['AAAAAAAACATTCGCGGTGCACTTCTTCGAT', '4', '*', '1\n']
#['AAAAAAAACATTCGCGGTGCACTTCTTGCGC', '4', '*', '1\n']
#['AAAAAAAACATTCGCGGTGCACTTCTTGGAT', '4', '*', '1\n']
#['AAAAAAAACATTCGCGGTGCACTTCTTGGCA', '4', '*', '1\n']
#['AAAAAAAACATTCGCGGTGCACTTCTTGGTC', '4', '*', '1\n']
#['AAAAAAAACATTCGCGGTGCACTTCTTGTAC', '4', '*', '1\n']
#['AAAAAAAACATTCGCGGTGCACTTCTTGTTG', '4', '*', '1\n']
