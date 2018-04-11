import pysam
import argparse
import csv
import sys


ap = argparse.ArgumentParser()
ap.add_argument('inBam', help='Mapped reads in bam format')
ap.add_argument('out', help='out')
ap.add_argument('chr', help='chr')

args = ap.parse_args()










out=open(args.out,"w")

reads_all=set()
reads_ed_0=set()
reads_multi_mapped=set()
reads_discordant=set()
reads_partially_mapped=set()
reads_ed_gt_1=set()
reads_PHREAD=set()
reads_good=set()



samfile = pysam.AlignmentFile(args.inBam, "rb" )




out.write("readName,refPos,readPos,refAllele,readAllele,phread,mapq,ed,n.alternative.mapping, portion.mapped,flag.mate_pair\n")


for read in samfile.fetch(args.chr):
    alignedRefPositions = read.get_reference_positions()

    if alignedRefPositions:
        #print ("----------------------------")


        refStart = alignedRefPositions[0]

        refSequence = read.get_reference_sequence()
        readSequence = read.query_alignment_sequence
        readQ=read.query_qualities

        ed = read.get_tag("NM")
        l_original=read.infer_query_length()
        l_mapped=len(readSequence)
        portion_mapped=l_mapped/float(l_original)
        flag_mate_pair=1-read.mate_is_unmapped
        mapq=read.mapping_quality

        if read.is_read1:
            readName=read.query_name + "/1"
        else:
            readName=read.query_name + "/2"

        if not read.has_tag("XA"):
            n_alternative_mapping = 1
        else:
            n_alternative_mapping = len(read.get_tag("XA").split(";"))


        if ed!=0:



            k=0
            for f, b,q in zip(refSequence, readSequence,readQ):
                k+=1
                if f!=b:
                    #print(readName,refStart+k,k, f,b,q,mapq,ed, n_alternative_mapping, portion_mapped,flag_mate_pair)
                    t_list=[]
                    t_list.clear()
                    t_list.append(readName)
                    t_list.append(str(refStart+k))
                    t_list.append(str(k))
                    t_list.append(str(f))
                    t_list.append(str(b))
                    t_list.append(str(q))
                    t_list.append(str(mapq))
                    t_list.append(str(ed))
                    t_list.append(str(n_alternative_mapping))

                    t_list.append("{:.2f}".format(portion_mapped))
                    t_list.append(str(flag_mate_pair))



                    t_string=','.join(t_list)
                    out.write(t_string)
                    out.write("\n")



















