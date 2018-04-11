import pysam
import argparse
import csv
import sys


ap = argparse.ArgumentParser()
ap.add_argument('inBam', help='Mapped reads in bam format')
ap.add_argument('rl', help='read length')
ap.add_argument('out', help='out')
ap.add_argument("-ed", "--ed", default=1,help="Min edit distance", type=int)



args = ap.parse_args()


cigar_full=str(args.rl)+"M"


#SRR062635.19983732	163	2	60058	60	100M	2	60129	100	TGGTAGGGGACAAAACAGATACAGCACCAAACATGAAAAGCATTGATTTATCTCTCTCTGCTAATGCATGTGAAGTGTTGATCCTGGGGTTGAAAGTCGG	array('B', [15, 35, 36, 32, 36, 36, 36, 34, 29, 37, 32, 37, 37, 37, 37, 34, 37, 36, 37, 36, 37, 34, 37, 36, 36, 37, 34, 36, 35, 37, 37, 32, 37, 36, 34, 37, 37, 37, 33, 34, 36, 37, 34, 37, 36, 33, 36, 37, 37, 36, 36, 36, 37, 36, 35, 36, 37, 32, 37, 33, 36, 37, 37, 38, 33, 34, 36, 34, 37, 36, 34, 35, 32, 35, 35, 32, 35, 32, 34, 34, 33, 32, 33, 33, 34, 33, 34, 34, 34, 33, 32, 34, 35, 30, 32, 33, 31, 29, 23, 30])	[('X0', 1), ('X1', 1), ('MD', '0G99'), ('RG', 'SRR062635'), ('AM', 23), ('NM', 1), ('SM', 23), ('MQ', 60), ('XT', 'U'), ('BQ', 'ARP@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')]	0	0	14636624	0	1	0	0



#MD tag will be always in this form : 0G99, Since we allow for only one mistmath. number{ACTG}number


snpSet=set()
dict={}


n=0

out=open(args.out,"w")

samfile = pysam.AlignmentFile(args.inBam, "rb" )

#XA are alternative alignments (~secondary hits).

for read in samfile.fetch("chr3"):
        if read.has_tag("NM"):
            ed=read.get_tag("NM")
            cigar=read.cigarstring
            phreadQ = pileupread.alignment.query_qualities[pileupread.query_position]


            if ed<=args.ed and cigar==cigar_full and not read.has_tag("XA") and not read.mate_is_unmapped and phreadQ>17:
                
                for i in read.get_aligned_pairs(with_seq=True):
                    
                    
                    
                    
                    
                    if i[2].islower():
                        if read.query_qualities[i[0]]>17:
                            ALT=i[2]
                            posRef=int(i[1])
                            #print i, posRef, ALT

                            if posRef not in snpSet:
                                dict[posRef]=[0,0,0,0]
                                if ALT=="a":
                                    dict[posRef][0]+=1
                                elif ALT=="c":
                                    dict[posRef][1]+=1
                                elif ALT=="t":
                                    dict[posRef][2]+=1
                                elif ALT=="g":
                                    dict[posRef][3]+=1
                            else:
                                if ALT=="a":
                                    dict[posRef][0]+=1
                                elif ALT=="c":
                                    dict[posRef][1]+=1
                                elif ALT=="t":
                                    dict[posRef][2]+=1
                                elif ALT=="g":
                                    dict[posRef][3]+=1
                            snpSet.add(posRef)








print "Number of SNPs", len(snpSet)
print "---------"


for k,v in dict.items():
    
        out.write("chr3,"+str(k)+","+str(v[0])+","+str(v[1])+","+str(v[2])+","+str(v[3]) )
        out.write("\n")




samfile.close()





