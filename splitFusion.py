import pysam
import sys
import argparse
import re

def parse_args():
    AP = argparse.ArgumentParser("get likely fusion by split reads")
    AP.add_argument('-bam',help='tumor bam',dest='bam')
    AP.add_argument('-mpq',help='map quality cutoff',dest='mpq',default=10)
    AP.add_argument('-od',help='output dir',dest='outdir')

    return AP.parse_args()

def count_sa_tag_num(SA_tag):
    num = 0
    for i in SA_tag:
        if i == ';':
            num += 1

    return num

def count_S_num(samflag):
    num = 0
    for i in samflag:
        if i == 'S':
            num += 1

    return num

def cigarTotuples(CIGAR):
    '''
    M   BAM_CMATCH      0
    I   BAM_CINS        1
    D   BAM_CDEL        2
    N   BAM_CREF_SKIP   3
    S   BAM_CSOFT_CLIP  4
    H   BAM_CHARD_CLIP  5
    P   BAM_CPAD        6
    =   BAM_CEQUAL      7
    X   BAM_CDIFF       8
    B   BAM_CBACK       9
    '''

    cigar_tuple = []

    cigar_s = []
    CIGAR_VALUE = 'MIDNSHP=XB'
    
    single_int = {}
    val = 0
    for i in CIGAR_VALUE:
        single_int[i] = val
        val += 1


    cigar_exists = []
    for s in CIGAR:
        if s in CIGAR_VALUE:
            cigar_exists.append(s)

    CIGAR_COPY = CIGAR

    for s in cigar_exists:
        v = single_int[s]
        idx = CIGAR_COPY.index(s) # 该字符第一次出现时的位置
        base_n = int(CIGAR_COPY[:idx]) # 碱基数
        CIGAR_COPY = CIGAR_COPY[(idx+1):] # 更新CIGAR
        tuple_val = (v,base_n)
        cigar_s.append(tuple_val)

        if len(CIGAR_COPY) == 0:
            break
    
    return(cigar_s)

        
def main():
    args = parse_args()
    samfile = pysam.AlignmentFile(args.bam,'rb')

    outfile = args.outdir + '/likely.fusion.xls'
    of = open(outfile,'w')
    header = ['']

    for read in samfile.fetch():
        #print(read)
        # get SECONDARY or SUPPLEMENTARY reads
        # Note: for BWA MEM's -M arg, split read will be marked as secondary read
        if read.mapping_quality < 10:
            continue

        cigar = read.cigarstring
        cigar_tuple = read.cigartuples
        #print(cigar_tuple)
        #print(cigar)

        if read.is_secondary or read.is_supplementary:
            # check SA tag
            if read.has_tag('SA'):
                sa_tag = read.get_tag('SA')
                sa_tag = sa_tag[:-1] # remove last ';'

                #print(sa_tag)
                '''
                sa_tag's value maybe:
                    13,25049714,-,64S87M,27,0;
                    11,94219169,-,33S84M34S,41,0;10,95507679,-,33M118S,41,0;
                '''
                if read.is_reverse:
                    read_p1_strand = '-'
                else:
                    read_p1_strand = '+'

                #print(read_p1_strand)

                # count ';' num
                read_p2_strand = []
                sa_num = count_sa_tag_num(sa_tag)
                #print(sa_num)
                if sa_num == 0:
                    # check mapq
                    q = int(sa_tag.split(',')[-2])
                    if q <= 10:
                        continue

                    # check 'H'
                    sa_cigar = sa_tag.split(',')[3]
                    #print(sa_cigar)
                    #if 'H' not in sa_cigar:
                    #    # only use 'H'
                    #    continue

                    S_num = count_S_num(sa_cigar)
                    if S_num > 1:
                        # skip flag if it contains >= 2 'S'
                        continue

                    # SA is 'S'
                    sa_strand = sa_tag.split(',')[2]
                    #print(sa_strand)

                    if sa_strand == '+':
                        if read_p1_strand == '+':
                            #print(read)
                            # ++
                            # check who is the hgene and tgene
                            if sa_cigar[-1] == 'S':
                                # sa is the hgene
                                # check if read_p1 is the tgene
                                if cigar[0] == 'H':
                                    # read_p1 is the tgene
                                    pass
                                    # output fusion info
                                    #print(read)
                                else:
                                    continue
                            elif sa_cigar[0] == 'S':
                                # sa is the tgene
                                # check if read_p1 is the hgene
                                if cigar[-1] == 'H':
                                    # read_p1 us the hgene

                                    # output fusion info
                                    #print(read)
                                    pass

                                else:
                                    continue
                            else:
                                continue
                        else:
                            # +-
                            # check who is the hgene and tgene
                            #print(read)
                            #print(cigar_tuple)
                            if sa_cigar[-1] == 'S':
                                #print(read)
                                #print(cigar_tuple)
                                # sa is the hgene
                                # check if read_p1 is the tgene
                                #if cigar[0] == 'H':

                                if cigar_tuple[-1][0] == 5:
                                    #print(read)
                                    #print(cigar_tuple)
                                    # H is the first CHR in CIGAR
                                    # read_p1 is the tgene
                                    hgene_chr = sa_tag.split(',')[0].split(':')[-1]
                                    hgene_left_align_pos = int(sa_tag.split(',')[1])

                                    # count how many ref base(s) consumed by CIGAR
                                    ref_consume_len = 0`
                                    cigar2tupele = cigarTotuples(sa_tag.split(',')[3])
                                    #print(cigar2tupele)

                                    for i in cigar2tupele:
                                        if i[0] == 0 or i[0] == 2 or i[0] == 3 or i[0] == 7 or i[0] == 8:
                                            ref_consume_len += i[1]

                                    hgene_breakpoint_pos = hgene_left_align_pos + ref_consume_len + 1 # add 1

                                    tgene_chr = read.reference_name
                                    tgene_left_align_pos = read.reference_start + 1 # 0-based

                                    ref_consume_len = 0
                                    for i in cigar_tuple:
                                        if i[0] == 0 or i[0] == 2 or i[0] == 3 or i[0] == 7 or i[0] == 8:
                                            ref_consume_len += i[1]

                                    tgene_breakpoint_pos = tgene_left_align_pos + ref_consume_len + 1 # add 1

                                    fs_info = read.query_name + '\t' + str(hgene_chr) + '\t' + str(hgene_breakpoint_pos) + '\t' + str(tgene_chr) + '\t' + str(tgene_breakpoint_pos)

                                    print(fs_info)

                                else:
                                    continue
                            elif sa_cigar[0] == 'S':
                                # sa is the tgene
                                # check if read_p1 is the hgene
                                if cigar[0] == 'H':
                                    # read_p1 is the hgene

                                    # output fusion info
                                    #print(read)
                                    pass
                                else:
                                    continue
                            else:
                                continue
                    else:
                        if read_p1_strand == '+':
                            # +-
                            # check who is the hgene and tgene
                            if sa_cigar[0] == 'S':
                                # sa is the hgene
                                # check if read_p1 is the tgene
                                if cigar[0] == 'H':
                                    # read_p1 is the tgene

                                    # output fusion info
                                    #print(read)
                                    pass


                                else:
                                    continue
                            elif sa_cigar[-1] == 'S':
                                # sa is the tgene
                                # check if read_p1 is the hgene
                                if cigar[-1] == 'H':
                                    # read_p1 is the hgene

                                    # output fusion info
                                    #print(read)
                                    pass
                            else:
                                continue
                        else:
                            # --
                            # check who is the hgene and tgene
                            if sa_cigar[0] == 'S':
                                # sa is the hgene
                                # check read_p1 is the tgene
                                if cigar[-1] == 'H':
                                    # read_p1 is the tgene

                                    # output fusion info
                                    #print(read)
                                    pass

                                else:
                                    continue
                            elif sa_cigar[-1] == 'S':
                                # sa is the tgene
                                # check if read_p1 is the hgene
                                if cigar[0] == 'H':
                                    # read_p1 is the hgene

                                    # output fusion info
                                    #print(read)
                                    pass
                                else:
                                    continue
                            else:
                                continue
                else:
                    sa_arr = sa_tag.split(';')
                    for sa_tag in sa_arr:
                        # ckeck mapq
                        q = int(sa_tag.split(',')[-2])
                        if q <= 10:
                            continue

            else:
                # do not have SA tag
                continue
        else:
            # not sec or sup
            continue



    samfile.close()
    of.close()

if __name__ == "__main__":
    main()