import pysam
import sys
import argparse
import re
import os

def parse_args():
    AP = argparse.ArgumentParser("check the uniqness of alignment of fusion support split reads")
    AP.add_argument('-bam',help='tumor bam',dest='bam')
    AP.add_argument('-f',help='*.fusion.xls file',dest='fsFile')
    AP.add_argument('-mq',help='map quality cutoff',dest='mapq',default=20)
    AP.add_argument('-refFlag',help='refFlat to annot gene strand info',dest='refFlat',default='/home/fulongfei/workdir/git_repo/GeExCNV/public_db/refFlat.txt')
    AP.add_argument('-od',help='output dir',dest='outdir')

    return AP.parse_args()

def check_read_cigar(gene_type,strand,cigar_tuple):
    '''
    根据下面两个条件判断CIGAR是否异常:
    1) hgene or tgene?
    2) forward or reverse strand?
    

    30M -> (0,30)
    
    M 0
    S 4

    '''

    cigar_check = ''

    if gene_type == 'hgene':
        if strand == '+': # hgene is on forward strand
            # 右端为S
            # EML4->ALK +-
            if cigar_tuple[-1][0] != 4:
                cigar_check = 'false'
            else:
                cigar_check = 'true'
        else:
            # 左端为S
            # KIF5B->RET -+
            if cigar_tuple[0][0] != 4:
                cigar_check = 'false'
            else:
                cigar_check = 'true'
    else:
        # for tgene
        if strand == '+':
            # 左端为S
            if cigar_tuple[0][0] != 4:
                cigar_check = 'false'
            else:
                cigar_check = 'true'
        else:
            # 右端为S
            if cigar_tuple[-1][0] != 4:
                cigar_check = 'false'
            else:
                cigar_check = 'true'

    return(cigar_check)
            


def check_sa_ori(gene_type,read_strand,sa_tag,hgene_strand,tgene_strand):
    '''
    根据下面条件判断SA +/-链是否异常:
    1) hgene strand
    2) tgene strand
    3) this read strand

    '''

    sa_ori_check = ''

    sa_strand = sa_tag.split(',')[2]
    sa_cigar = sa_tag.split(',')[3]
    sa_cigar_tuple = cigarTotuples(sa_cigar)


    if gene_type == 'hgene':
        if hgene_strand == '+' and tgene_strand == '-':
            # EML4->ALK
            # https://cancer.sanger.ac.uk/cosmic/fusion/overview?fid=180308&gid=150139

            if read_strand == 'FORWARD':
                # split read need on '-'
                if sa_strand == '-':
                    sa_ori_check = 'true'
                else:
                    sa_ori_check = 'false'
            else:
                # split read need on '+'
                if sa_strand == '+':
                    sa_ori_check = 'true'
                else:
                    sa_ori_check = 'false'


        if hgene_strand == '+' and tgene_strand == '+':
            # FGFR3->TACC3
            # https://cancer.sanger.ac.uk/cosmic/fusion/overview?fid=90231&gid=9677

            if read_strand == 'FORWARD':
                if sa_strand == '+':
                    sa_ori_check = 'true'
                else:
                    sa_ori_check = 'false'
            else:
                if sa_strand == '-':
                    sa_ori_check = 'true'
                else:
                    sa_ori_check = 'false'


        if hgene_strand == '-' and tgene_strand == '-':
            # CD74->ROS1
            # https://cancer.sanger.ac.uk/cosmic/fusion/overview?fid=141941&gid=72432

            if read_strand == 'FORWARD':
                if sa_strand == '+':
                    sa_ori_check = 'true'
                else:
                    sa_ori_check = 'false'
            else:
                if sa_strand == '-':
                    sa_ori_check = 'true'
                else:
                    sa_ori_check = 'false'


        if hgene_strand == '-' and tgene_strand == '+':
            # KIF5B->RET
            # https://cancer.sanger.ac.uk/cosmic/fusion/overview?fid=9508&gid=190031

            if read_strand == 'FORWARD':
                if sa_strand == '-':
                    sa_ori_check = 'true'
                else:
                    sa_ori_check = 'false'
            else:
                if sa_strand == '+':
                    sa_ori_check = 'true'
                else:
                    sa_ori_check = 'false'
    else:
        # for tgene
        if hgene_strand == '+' and tgene_strand == '-':
            # EML4->ALK

            if read_strand == 'FORWARD':
                if sa_strand == '-':
                    sa_ori_check = 'true'
                else:
                    sa_ori_check = 'false'
            else:
                if sa_strand == '+':
                    sa_ori_check = 'true'
                else:
                    sa_ori_check = 'false'


        if hgene_strand == '+' and tgene_strand == '+':
            # FGFR3->TACC3
            # https://cancer.sanger.ac.uk/cosmic/fusion/overview?fid=90231&gid=9677

            if read_strand == 'FORWARD':
                if sa_strand == '+':
                    sa_ori_check = 'true'
                else:
                    sa_ori_check = 'false'
            else:
                if sa_strand == '-':
                    sa_ori_check = 'true'
                else:
                    sa_ori_check = 'false'

        if hgene_strand == '-' and tgene_strand == '-':
            # CD74->ROS1
            # https://cancer.sanger.ac.uk/cosmic/fusion/overview?fid=141941&gid=72432

            if read_strand == 'FORWARD':
                if sa_strand == '+':
                    sa_ori_check = 'true'
                else:
                    sa_ori_check = 'false'
            else:
                if sa_strand == '-':
                    sa_ori_check = 'true'
                else:
                    sa_ori_check = 'false'

        if hgene_strand == '-' and tgene_strand == '+':
            # KIF5B->RET
            # https://cancer.sanger.ac.uk/cosmic/fusion/overview?fid=9508&gid=190031

            if read_strand == 'FORWARD':
                if sa_strand == '-':
                    sa_ori_check = 'true'
                else:
                    sa_ori_check = 'false'
            else:
                if sa_strand == '+':
                    sa_ori_check = 'true'
                else:
                    sa_ori_check = 'false'

    return(sa_ori_check)
               

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

def cal_ref_consume_len(cigar_tuple):
    ref_consume_len = 0
	
    '''
    CIGAR is a little complex to parse
        58S30M3I60M
        53S45M4S
        50S100M
        100M50S
    '''
    
    if cigar_tuple[0][0] == 4 and cigar_tuple[-1][0] != 4:
        # S开头，且末尾不为S
        # 50S100M or 58S30M3I60M
        ref_consume_len = 0
    elif cigar_tuple[0][0] == 4 and cigar_tuple[-1][0] == 4:
        # S开头S结尾
        # 53S45M4S or 4S45M53S
        # 比较开头、结尾的长度
        if int(cigar_tuple[0][1]) > int(cigar_tuple[-1][1]): # 53S45M4S
            ref_consume_len = 0
        else: # 4S45M53S
            for i in cigar_tuple:
                if i[0] == 0 or i[0] == 2 or i[0] == 3 or i[0] == 7 or i[0] == 8:
                    ref_consume_len += i[1]
    elif cigar_tuple[0][0] == 0 and cigar_tuple[-1][0] == 4:
        # M开头S结尾
		# 100M50S
        for i in cigar_tuple:
            if i[0] == 0 or i[0] == 2 or i[0] == 3 or i[0] == 7 or i[0] == 8:
                ref_consume_len += i[1]
    else:
        # 其余特殊情况，暂时还未发现
        ref_consume_len = 0
    
    return(ref_consume_len)
            

def main():
    args = parse_args()
    fsFH = open(args.fsFile,'r')
    name = os.path.basename(args.fsFile).split('.')[0]
    outfile = '%s/%s.sr_mapq_uniqness_check.txt' % (args.outdir,name)

    ofFH = open(outfile,'w')

    outfile_header = ['fusionGene','breakpoint_info','Gene','Hgene_OR_Tgene','read_name','read_info','cigar_check','mate_chr','mate_chr_check','SA_info','SR_mapq','SR_uniq_check','sa_ori_check','hgene_sr_pos','tgene_sr_pos','fs_pos_check(split.pos.VS.sv.pos','final_check']
    ofFH.write('\t'.join(outfile_header)+'\n') # write outfile header

    # get gene strand info
    gene_strand = {}
    refFlatFH = open(args.refFlat,'r')
    for line in refFlatFH:
        arr = line.strip().split('\t')
        gene_strand[arr[0]] = arr[3]

    refFlatFH.close()

    gene_strand['SEPTIN14'] = '-' # SEPTIN14 is alias to SEPT14
    gene_strand['SEPT14'] = '-'

    gene_strand['WISP1'] = '+' # # WISP1 is alias to CCN4
    gene_strand['CCN4'] = '+'


    

    for line in fsFH:
        line_count = 0 # 如果fsFH中的基因没有找到一条split read，也要在结果文件中输出该基因的信息

        # each line is a fusion candidate
        arr = line.strip().split('\t')
        if arr[0] == 'SampleID':
            # skip header
            continue

        if arr[9] != 'YES':
            # skip not report fusion
            continue

        print("[INFO]: check %s fusion info..." % (arr[1]))

        # 确定hgene/tgene
        if arr[1].count('-') == 1:
            hgene = arr[1].split('-')[0] # not use ->
            tgene = arr[1].split('-')[1]
        elif arr[1].count('-') == 2:
            # 处理特殊基因名
            # FMO4-PLCG1-AS1,hgene is FM04, tgene is PLCG1-AS1 (PLCG1-AS1是一个RNA基因，可以过滤掉）
            gene_pair = []
            gene_pair.append(arr[2])
            gene_pair.append(arr[3])

            hg1 = arr[1].split('-')[0]
            tg1 = arr[1].split('-')[1] + '-' + arr[1].split('-')[2]

            hg2 = arr[1].split('-')[0] + '-' + arr[1].split('-')[1]
            tg2 = arr[1].split('-')[2]

            if hg1 in gene_pair and tg1 in gene_pair:
                hgene = hg1
                tgene = tg1
            elif hg2 in gene_pair and tg2 in gene_pair:
                hgene = hg2
                tgene = tg2
            else:
                print("arr[1] skipped")
                continue
        else:
            print("arr[1] skipped")
            continue

        # 染色体信息
        fs_chr = {}
        fs_chr[arr[11]] = arr[12] # gene => chr
        fs_chr[arr[17]] = arr[18]


        # 断点信息
        fs_pos = {}
        fs_pos[arr[11]] = int(arr[13]) # gene -> pos
        fs_pos[arr[17]] = int(arr[19])

        # +/-链信息
        if hgene in gene_strand:
            hgene_strand = gene_strand[hgene]
        else:
            print("%s skipped for hgene do not have strand info" % (arr[1]))
            continue

        if tgene in gene_strand:
            tgene_strand = gene_strand[tgene]
        else:
            print("%s skipped for tgene do not have strand info" % (arr[1]))
            continue
        

        fs_info = "%s,%s;%s;%s" % (fs_chr[hgene],fs_pos[hgene],fs_chr[tgene],fs_pos[tgene]) # chr,pos;chr,pos


        samfile = pysam.AlignmentFile(args.bam,'rb')
        
        ########################################### 处理hgene ###########################################
        
        # hgene 断点
        hgene_bp = fs_pos[hgene]
        tgene_bp = fs_pos[tgene]
        hgene_left = hgene_bp - 10
        hgene_right = hgene_bp + 10

        print("[INFO]: check Hgene %s (%s) info..." % (hgene,hgene_strand))

        for read in samfile.fetch(arr[12],hgene_left,hgene_right):
            #print(read)
            
            if read.is_unmapped:
                continue

            if read.is_secondary or read.is_supplementary or read.is_duplicate:
                continue

            aln_chr = read.reference_name
            #mate_chr = read.next_reference_name
            
            if read.is_read1:
                r1r2 = 'READ1'
            else:
                r1r2 = 'READ2'

            cigarS = read.cigarstring

            if read.is_reverse:
                read_rev = 'REVERSE'
            else:
                read_rev = 'FORWARD'

            MAPQ = read.mapping_quality
            
            if MAPQ  < args.mapq: # default: 20
                continue
            
            read_info = ';'.join([str(aln_chr),str(read.reference_start),r1r2,read_rev,cigarS,str(MAPQ)]) # chr;left_pos,read1/2;reverse/forward;CIGAR;MAPQ;
            

            if read.has_tag('SA'):
                #print(read.query_name)
                sa_tag = read.get_tag('SA')
                # count how many SA
                sa_count = sa_tag.count(';')
                sa_strand = sa_tag.split(',')[2] # +/-

                # check read cigar
                read_cigar_check = ''
                cigar_check = check_read_cigar('hgene',hgene_strand,read.cigartuples)
                if cigar_check == 'false':
                    read_cigar_check = 'NOPASS'
                else:
                    read_cigar_check = 'PASS'

                sa_arr = sa_tag.split(';')
                for SA in sa_arr:
                    if len(SA) == 0:
                        continue # skip empty

                    sr_mapq = int(SA.split(',')[4])
                    sa_chr = SA.split(',')[0]

                    # 检查比对唯一性, default: 20
                    if sr_mapq >= args.mapq:
                        sr_uniq_check = 'PASS'
                    else:
                        sr_uniq_check = 'NOPASS'

                    # 检查SA ori
                    sa_ori = check_sa_ori('hgene',read_rev,SA,hgene_strand,tgene_strand)
                    
                    if sa_ori == 'false':
                        sa_ori_check = 'NOPASS'
                    else:
                        sa_ori_check = 'PASS'

                    # 计算hgene断点位置
                    left_pos = read.reference_start + 1 # 0-based leftmost coordinate
                    cigar_tuple = read.cigartuples
                    ref_consume_len = cal_ref_consume_len(cigar_tuple)
                    
                    hgene_fs_pos = left_pos + ref_consume_len - 1 + 1
                    hgene_fs_info = read.reference_name + ',' + str(hgene_fs_pos)
                    
                    # 计算tgene断点位置
                    tgene_left_pos = int(SA.split(',')[1])
                    sa_cigar = SA.split(',')[3]
                    sa_chr = SA.split(',')[0]
                    sa_cigar_tuple = cigarTotuples(sa_cigar)
                    
                    ref_consume_len = cal_ref_consume_len(sa_cigar_tuple)
                    tgene_fs_pos = tgene_left_pos + ref_consume_len - 1 + 1
                    tgene_fs_info = sa_chr + ',' + str(tgene_fs_pos)

                    # for hgene, skip if the the split pos is far away from the sv pos
                    # if the split pos is close the sv pos, then check the tgene's split pos

                    hgene_fs_chr = read.reference_name
                    hgene_sv_chr = fs_chr[hgene]
                    tgene_sv_chr = fs_chr[tgene]

                    if hgene_fs_chr == hgene_sv_chr and abs(hgene_fs_pos - hgene_bp) <= 10: # 保证根据split read计算的hgene的chr、pos与SV结果一致
                        # check sa chr and pos
                        if sa_chr == tgene_sv_chr and abs(tgene_fs_pos - tgene_bp) <= 10: # 保证split read计算的tgene的chr、pos与SV结果一致
                            fs_pos_check = 'PASS'
                        else:
                            fs_pos_check = 'NOPASS'
                    else:
                        # skip because these reads are not the fusion-support
                        print("%s skipped for chr conflict (SV chr & SR chr) OR pos conflict (SV pos VS SR pos)" % (read.query_name))
                        continue
                    
                    # check chr set
                    if read.next_reference_name not in [sa_chr,read.reference_name]:
                        print("%s skipped: mate chr is not included by this read and its split read! chrom conflict" % (read.query_name))
                        continue # 如果read1比对到chr6,sup read比对到chr5，read2比对到chr7，那么这个是有问题的
                    
                    chr_check = 'PASS' #
                    
                    # check sr uniq and fusion pos
                    if sr_uniq_check == 'PASS' and fs_pos_check == 'PASS' and read_cigar_check == 'PASS' and sa_ori_check == 'PASS':
                        final_check = 'PASS'
                    else:
                        final_check = 'NOPASS'

                    
                    ofFH.write('\t'.join([arr[1],fs_info,hgene,'hgene',str(read.query_name),read_info,read_cigar_check,read.next_reference_name,chr_check,sa_tag,str(sr_mapq),sr_uniq_check,sa_ori_check,str(hgene_fs_info),str(tgene_fs_info),fs_pos_check,final_check])+'\n')

                    line_count += 1

        # finished hgene screen
        
        ########################################### 处理tgene ###########################################

        # tgene 断点
        hgene_bp = fs_pos[hgene]
        tgene_bp = fs_pos[tgene]
        tgene_left = tgene_bp - 10
        tgene_right = tgene_bp + 10

        print("[INFO]: check Tgene %s (%s) info..." % (tgene,tgene_strand))
        
        for read in samfile.fetch(arr[18],tgene_left,tgene_right):
            if read.is_secondary or read.is_supplementary or read.is_duplicate:
                continue

            if read.is_unmapped:
                continue
            
            aln_chr = read.reference_name

            if read.is_read1:
                r1r2 = 'READ1'
            else:
                r1r2 = 'READ2'

            cigarS = read.cigarstring

            if read.is_reverse:
                read_rev = 'REVERSE'
            else:
                read_rev = 'FORWARD'

            MAPQ = read.mapping_quality

            if MAPQ < args.mapq: # default: 20
                continue

            read_info = ';'.join([str(aln_chr),str(read.reference_start),r1r2,read_rev,cigarS,str(MAPQ)]) # chr;left_pos,read1/2;reverse/forward;CIGAR;MAPQ;
            
            if read.has_tag('SA'):
                #print(read.query_name)
                sa_tag = read.get_tag('SA')
                sa_strand = sa_tag.split(',')[2] # +/-
                # count how many SA
                sa_count = sa_tag.count(';')

                # check read cigar

                read_cigar_check = ''
                cigar_check = check_read_cigar('tgene',tgene_strand,read.cigartuples)
                if cigar_check == 'false':
                    read_cigar_check = 'NOPASS'
                else:
                    read_cigar_check = 'PASS'


                sa_arr = sa_tag.split(';')
                for SA in sa_arr:
                    if len(SA) == 0:
                        continue # skip empty
                        '''

                        >>> a='2,42553276,+,45M106S,60,0;2,29447572,-,37M114S,60,0;'
                        >>> a.split(';')
                        ['2,42553276,+,45M106S,60,0', '2,29447572,-,37M114S,60,0', '']
                        >>> a.split(';')[-1]
                        ''
                        >>> len(a.split(';')[-1])
                        0
                        >>> len(a.split(';')[-2])
                        25

                        '''

                    sr_mapq = int(SA.split(',')[4])
                    sa_chr = SA.split(',')[0]

                    # 检查比对唯一性,default: 20
                    if sr_mapq >= args.mapq:
                        sr_uniq_check = 'PASS'
                    else:
                        sr_uniq_check = 'NOPASS'

                    # 检查SA ori
                    sa_ori_check = ''
                    sa_ori = check_sa_ori('tgene',read_rev,SA,hgene_strand,tgene_strand)
                    if sa_ori == 'false':
                        sa_ori_check = 'NOPASS'
                    else:
                        sa_ori_check = 'PASS'


                    # 计算tgene断点位置
                    left_pos = read.reference_start + 1 # 0-based leftmost coordinate
                    cigar_tuple = read.cigartuples
                    ref_consume_len = cal_ref_consume_len(cigar_tuple)
                    tgene_fs_pos = left_pos + ref_consume_len - 1 + 1
                    tgene_fs_info = read.reference_name + ',' + str(tgene_fs_pos)

                    
                    # 计算hgene断点位置
                    hgene_left_pos = int(SA.split(',')[1])
                    sa_cigar = SA.split(',')[3]
                    sa_chr = SA.split(',')[0]
                    sa_cigar_tuple = cigarTotuples(sa_cigar)
                    ref_consume_len = cal_ref_consume_len(sa_cigar_tuple)

                    hgene_fs_pos = hgene_left_pos + ref_consume_len - 1 + 1
                    hgene_fs_info = sa_chr + ',' + str(hgene_fs_pos)


                    # for tgene, skip if the the split pos is far away from the sv pos
                    # if the split pos is close the sv pos, then check the hgene's split pos

                    tgene_fs_chr = read.reference_name
                    hgene_sv_chr = fs_chr[hgene]
                    tgene_sv_chr = fs_chr[tgene]

                    if tgene_fs_chr == tgene_sv_chr and abs(tgene_fs_pos - tgene_bp) <= 10:
                        # check sa chr and pos for hgene
                        if sa_chr == hgene_sv_chr and abs(hgene_fs_pos - hgene_bp) <= 10:
                            fs_pos_check = 'PASS'
                        else:
                            fs_pos_check = 'NOPASS'
                    else:
                        #print("%s was skipped for tgene fusion chr and pos are not exact" % (read.query_name))
                        continue

                    # check chr set
                    if read.next_reference_name not in [sa_chr,read.reference_name]:
                        continue
                    
                    chr_check = 'PASS'
 
                    # check sr uniq and fusion pos
                    if sr_uniq_check == 'PASS' and fs_pos_check == 'PASS' and read_cigar_check == 'PASS' and sa_ori_check == 'PASS':
                        final_check = 'PASS'
                    else:
                        final_check = 'NOPASS'


                    # output
                    ofFH.write('\t'.join([arr[1],fs_info,tgene,'tgene',str(read.query_name),read_info,read_cigar_check,read.next_reference_name,chr_check,SA,str(sr_mapq),sr_uniq_check,sa_ori_check,str(hgene_fs_info),str(tgene_fs_info),fs_pos_check,final_check])+'\n')
                    line_count += 1
        samfile.close()

        if line_count == 0:
            ofFH.write('\t'.join([arr[1],fs_info,hgene,'hgene','query_name','read_info','read_cigar_check','mate_chr','chr_check','SA','sr_mapq','sr_uniq_check','sa_ori_check','hgene_fs_info','tgene_fs_info','fs_pos_check','NOPASS'])+'\n')

    fsFH.close()
    ofFH.close()




if __name__ == "__main__":
    main()



