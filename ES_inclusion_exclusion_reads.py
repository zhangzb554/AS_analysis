#!/python
import re
import sys
import os
########################################inclusion or exclusion for a given read############################################
def ES_inclusion(j_ref1_l, j_ref1_r, j_ref2_l, j_ref2_r, read_start, cigar):# 5bp overhangs#
    junc_l_count,junc_r_count,ex_junc_count = 0,0,0
    intron_ref1 = [j_ref1_l + 1, j_ref1_r - 1]
    intron_ref2 = [j_ref2_l + 1, j_ref2_r - 1]
    intron_alt = [j_ref1_l + 1, j_ref2_r - 1]
    exon_mid_len = int(j_ref2_l) - int(j_ref1_r) + 1
    matches = re.findall(r'(\d+)(\w)', cigar)
    matches = [i for i in matches if "S" not in i]
    matches = list(map(lambda x: (int(x[0]),x[1]),matches))
    ops = [i[1] for i in matches]
    if len(ops) > ops.count("M") + ops.count("N"):
        return 0,0,0
    m_blocks, n_blocks = [], []
    #left_overhangs, right_overhangs = [], []
    i = 0
    for match in matches:
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
    if intron_ref1 in n_blocks and intron_ref2 in n_blocks:
        intron_ref1_index = n_blocks.index(intron_ref1)
        intron_ref2_index = n_blocks.index(intron_ref2)
        if intron_ref2_index - intron_ref1_index != 1:
            return 0,0,0
        left_overhang = m_blocks[intron_ref1_index][1] - m_blocks[intron_ref1_index][0]
        right_overhang = m_blocks[intron_ref2_index + 1][1] - m_blocks[intron_ref2_index + 1][0]
        if left_overhang >= 5 and right_overhang >= 5:
            #11111111111------------11111111111------------11111111111
            #RRRRRRRRRRR------------RRRRRRRRRRR------------RRRRRRRRRRR
            return 1,1,0
        else:
            return 0,0,0#useless read, so treated as 0,0,0
    elif intron_ref1 in n_blocks:
        intron_ref1_index = n_blocks.index(intron_ref1)
        left_overhang = m_blocks[intron_ref1_index][1] - m_blocks[intron_ref1_index][0]
        right_overhang = m_blocks[intron_ref1_index + 1][1] - m_blocks[intron_ref1_index + 1][0]
        if left_overhang >= 5 and right_overhang >= 5 and right_overhang <= exon_mid_len:
            #11111111111------------11111111111------------11111111111
            #RRRRRRRRRRR------------RRRRRRRR
            return 1,0,0
        else:
            return 0,0,0#useless read, so treated as 0,0,0
    elif intron_ref2 in n_blocks:
        intron_ref2_index = n_blocks.index(intron_ref2)
        left_overhang = m_blocks[intron_ref2_index][1] - m_blocks[intron_ref2_index][0]
        right_overhang = m_blocks[intron_ref2_index + 1][1] - m_blocks[intron_ref2_index + 1][0]
        if left_overhang >= 5 and right_overhang >= 5 and left_overhang <= exon_mid_len:
            #11111111111------------11111111111------------11111111111
            #                          RRRRRRRR------------RRRRRRRRRRR
            return 0,1,0
        else:
            return 0,0,0#useless read, so treated as 0,0,0
    elif intron_alt in n_blocks:
        intron_alt_index = n_blocks.index(intron_alt)
        left_overhang = m_blocks[intron_alt_index][1] - m_blocks[intron_alt_index][0]
        right_overhang = m_blocks[intron_alt_index + 1][1] - m_blocks[intron_alt_index + 1][0]
        if left_overhang >= 5 and right_overhang >=5:
            #11111111111------------11111111111------------11111111111
            #RRRRRRRRRRR-----------------------------------RRRRRRRRRRR
            return 0,0,1
        else:
            return 0,0,0#useless read, so treated as 0,0
    else:
        return 0,0,0#useless read, so treated as 0,0

#############################for a ES event, return the number of inclusion and exclusion and PSI value##################################
def inc_exc(j_ref1_l, j_ref1_r,j_ref2_l, j_ref2_r,sam):
    junc_l_count_tot,junc_r_count_tot,ex_junc_count_tot = 0,0,0
    with open(sam) as f:
        for l in f:
            inf = l.rstrip().split("\t")
            read_start,cigar = int(inf[3]),inf[5]
            junc_l_count,junc_r_count,ex_junc_count = ES_inclusion(j_ref1_l, j_ref1_r,j_ref2_l, j_ref2_r,read_start,cigar)
            junc_l_count_tot = junc_l_count_tot + junc_l_count
            junc_r_count_tot = junc_r_count_tot + junc_r_count
            ex_junc_count_tot = ex_junc_count_tot + ex_junc_count
    n = ex_junc_count_tot * 2
    d = junc_l_count_tot + junc_r_count_tot + ex_junc_count_tot * 2

    ##calculate PSI
    if (junc_l_count_tot >= 4 and junc_r_count_tot >= 4) or ex_junc_count_tot >=4:
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
    return  str(junc_l_count_tot), str(junc_r_count_tot), str(ex_junc_count_tot), str(psi)
#######################################main########################################################
bam_in = sys.argv[1]#"RNA-seq.bam"
ioe_SE = sys.argv[2]#SE.ioe
sam_out = "tmp.sam"
samtools = "/exdata/Users/zhibin/program/anaconda2/bin/samtools"
print "\t".join(["info",'j1','j2','ex','psi'])
with open(ioe_SE) as fh:
    for line in fh:
        if "seqname" in line:
            continue
        info = line.split("\t")
        lis = re.search('SE:(.*?):(\d+)-(\d+):(\d+)-(\d+)',info[2]).groups()
        chr, j_ref1_l, j_ref1_r,j_ref2_l, j_ref2_r = lis[0], int(lis[1]), int(lis[2]), int(lis[3]), int(lis[4])
        os.system("{0} view {1} {2}:{3}-{4} > {5}".format(samtools,bam_in,chr,j_ref1_l,j_ref2_r,sam_out))
        print(info[2] + "\t" + "\t".join(list(inc_exc(j_ref1_l, j_ref1_r,j_ref2_l, j_ref2_r,sam_out))))
        #break
