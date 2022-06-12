#!/python
import re
import sys
import os
########################################inclusion or exclusion for a given read############################################
#intron is larger than 20bp
#j_l = 4338022#left side of intron
#j_r = 4338126#right side of intron
#read_start = 4337091
#cigar = "83M788N61M103N6M"
#IR_inclusion(4338022, 4338126,4337091,"83M788N61M103N6M")

def IR_inclusion(j_l, j_r,read_start,cigar):# 5bp overhangs#
    in_l_count,in_r_count,junc_count = 0,0,0#IR_left, IR_right and intron_junction
    matchs = re.findall(r'(\d+)(\w)', cigar)
    matchs = [i for i in matchs if "S" not in i]
    matchs = list(map(lambda x: (int(x[0]),x[1]),matchs))
    ops = [i[1] for i in matchs]
    if len(ops) > ops.count("M") + ops.count("N"):
        return 0,0,0
    m_blocks, n_blocks = [], []
    #left_overhangs, right_overhangs = [], []
    i = 0
    for match in matchs:
        i += 1
        if i == 1:
            start,end = read_start, read_start + match[0] - 1
            m_blocks.append([start,end])
        else:
            start = end + 1
            end = start + match[0] -1
            if match[1] == 'M':
                m_blocks.append([start,end])
            else:
                n_blocks.append([start,end])
    ##############non-splicing read########################
    ### 1111111-------1111111    1111111-------1111111      1111111-------1111111
    ### RRRRRRRRRRR                       RRRRRRRRRRRR         RRRRRRRRRRRRRRRR
    if len(m_blocks) == 1:
        m_block = m_blocks[0]
        if m_block[0] < j_l and m_block[1] > j_r:
            return 1,1,0
        elif m_block[0] < j_l and m_block[1] >= j_l + 1 + 5:
            return 1,0,0
        elif m_block[0] <= j_r - 1 - 5 and m_block[1] > j_r:
            return 0,1,0
        else:
            return 0,0,0
    ###########################splicing read############################
    else:
        ####splicing junction equal to intron junction####
        if [j_l + 1, j_r -1] in n_blocks:
            intron_index = n_blocks.index([j_l + 1, j_r -1])
            left_overhang = m_blocks[intron_index][1] - m_blocks[intron_index][0]
            right_overhang = m_blocks[intron_index + 1][1] - m_blocks[intron_index + 1][0]
            if left_overhang >= 5 and right_overhang >=5:
                #1111111-------1111111
                #RRRRRRR-------RRRRRRR
                return 0,0,1
            else:
                return 0,0,0#useless read, so treated as 0,0,0
        else:
            if_over_intron = 0
            for n_block in n_blocks:
                ####splicing junction overlap to intron junction####
                if n_block[0] <= j_r - 1 and n_block[1] >= j_l + 1:
                    if_over_intron = 1
                    break
            if if_over_intron == 1:
                #1111111-------1111111    1111111-------1111111      1111111-------1111111
                #..RRR----RRR.........    .........RRR----RRR..      ..RRR------------RR..
                return 0,0,0#useless read, so treated as 0,0,0
            else:
                #like non-splicing read
                type = 0,0,0
                for m_block in m_blocks:
                    if m_block[0] < j_l and m_block[1] > j_r:
                        type =  1,1,0
                    elif m_block[0] < j_l and m_block[1] >= j_l + 1 + 5:
                        type = 1,0,0
                    elif m_block[0] <= j_r - 1 - 5 and m_block[1] > j_r:
                        type = 0,1,0
                    else:
                        continue
                return type
#############################for a IR event, return the number of inclusion and exclusion and PSI value##################################
def inc_exc(j_l,j_r,sam):
    in_l_count_tot, in_r_count_tot, j1_count_tot = 0,0,0
    with open(sam) as f:
        for l in f:
            inf = l.rstrip().split("\t")
            read_start,cigar = int(inf[3]),inf[5]
            in_l_count,in_r_count,j1_count = IR_inclusion(j_l,j_r,read_start,cigar)
            in_l_count_tot = in_l_count_tot + in_l_count
            in_r_count_tot = in_r_count_tot + in_r_count
            j1_count_tot = j1_count_tot + j1_count
            #if in_r_count == 1:
                #print l
    n = in_l_count_tot + in_r_count_tot
    d = in_l_count_tot + in_r_count_tot + j1_count_tot * 2
    ##calculate PSI
    if (in_l_count_tot >= 4 and in_r_count_tot >= 4) or j1_count_tot >=4:
        psi = float(n)/d
    else:
        psi = "NA"
    ## if PSI is meaningful
    #if in_l_count_tot < 4 and in_r_count_tot < 4 and j1_count_tot < 4:
        #if_conf = "amb"
    #elif in_l_count_tot >= 4 or in_r_count_tot >=4:
        #if_conf = "yes"
    #else:
        #if_conf = "no"
    return  str(in_l_count_tot), str(in_r_count_tot), str(j1_count_tot), str(psi)
###############################################main###############################################################
#######################################main########################################################
bam_in = sys.argv[1]#"RNA-seq.bam"
ioe_IR = sys.argv[2]#sca.AS.events_RI_strict.ioe
sam_out = "tmp.sam"
samtools = "/exdata/Users/zhibin/program/anaconda2/bin/samtools"
print "\t".join(["info",'in_l','in_r','J1','psi'])
with open(ioe_IR) as fh:
    for line in fh:
        if "seqname" in line:
            continue
        info = line.split("\t")
        lis = re.search('RI:(.*?):(\d+):(\d+)-(\d+):(\d+)',info[2]).groups()
        chr, j_l, j_r = lis[0], int(lis[2]), int(lis[3])
        os.system("{0} view {1} {2}:{3}-{4} > {5}".format(samtools,bam_in,chr,j_l,j_r,sam_out))
        print(info[2] + "\t" + "\t".join(list(inc_exc(j_l,j_r,sam_out))))
        #break
