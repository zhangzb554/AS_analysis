#!/python
import re
import sys
import os
########################################inclusion or exclusion for a given read############################################
def A3_inclusion(j_ref_l, j_ref_r, j_alt_l, j_alt_r, read_start, cigar):# 5bp overhangs#
    in_count, ex_count = 0, 0
    intron_ref = [j_ref_l + 1, j_ref_r - 1]
    intron_alt = [j_alt_l + 1, j_alt_r - 1]
    matches = re.findall(r'(\d+)(\w)', cigar)
    matches = [i for i in matches if "S" not in i]
    matches = list(map(lambda x: (int(x[0]),x[1]),matches))
    ops = [i[1] for i in matches]
    if len(ops) > ops.count("M") + ops.count("N"):
        return 0,0
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
    if intron_ref in n_blocks:
        intron_ref_index = n_blocks.index(intron_ref)
        left_overhang = m_blocks[intron_ref_index][1] - m_blocks[intron_ref_index][0]
        right_overhang = m_blocks[intron_ref_index + 1][1] - m_blocks[intron_ref_index + 1][0]
        if left_overhang >= 5 and right_overhang >= 5:
            return 0,1
        else:
            return 0,0#useless read, so treated as 0,0,0
    elif intron_alt in n_blocks:
        intron_alt_index = n_blocks.index(intron_alt)
        left_overhang = m_blocks[intron_alt_index][1] - m_blocks[intron_alt_index][0]
        right_overhang = m_blocks[intron_alt_index + 1][1] - m_blocks[intron_alt_index + 1][0]
        if left_overhang >= 5 and right_overhang >= 5:
            return 1,0
        else:
            return 0,0#useless read, so treated as 0,0,0
    else:
        return 0,0
#############################for a A5 event, return the number of inclusion and exclusion and PSI value##################################
def inc_exc(j_ref_l, j_ref_r, j_alt_l, j_alt_r,sam):
    in_count_tot, ex_count_tot = 0,0
    with open(sam) as f:
        for l in f:
            inf = l.rstrip().split("\t")
            read_start,cigar = int(inf[3]),inf[5]
            in_count, ex_count = A3_inclusion(j_ref_l, j_ref_r, j_alt_l, j_alt_r,read_start,cigar)
            in_count_tot = in_count_tot + in_count
            ex_count_tot = ex_count_tot + ex_count
    ##calculate PSI
    if in_count_tot >= 4 or ex_count_tot >=4:
        psi = float(in_count_tot)/(in_count_tot + ex_count_tot)
    else:
        psi = "NA"
    ## if PSI is meaningful
    #if in_l_count_tot < 4 and in_r_count_tot < 4 and j1_count_tot < 4:
        #if_conf = "amb"
    #elif in_l_count_tot >= 4 or in_r_count_tot >=4:
        #if_conf = "yes"
    #else:
        #if_conf = "no"
    return  str(in_count_tot), str(ex_count_tot), str(psi)

###############################################main###############################################################
bam_in = sys.argv[1]#"test.bam"
ioe_A3 = sys.argv[2]#A3
sam_out = "tmp.sam"
samtools = "/exdata/Users/zhibin/program/anaconda2/bin/samtools"
print "\t".join(["info",'in','ex','psi'])
with open(ioe_A3) as fh:
    for line in fh:
        if "seqname" in line:
            continue
        info = line.split("\t")
        lis = re.search('A3:(.*?):(\d+)-(\d+):(\d+)-(\d+)',info[2]).groups()
        chr, j_ref_l, j_ref_r, j_alt_l, j_alt_r = lis[0], int(lis[1]), int(lis[2]), int(lis[3]), int(lis[4])
        os.system("{0} view {1} {2}:{3}-{4} > {5}".format(samtools,bam_in,chr,min(j_ref_l, j_ref_r, j_alt_l, j_alt_r), max(j_ref_l, j_ref_r, j_alt_l, j_alt_r), sam_out))
        print(info[2] + "\t" + "\t".join(list(inc_exc(j_ref_l, j_ref_r, j_alt_l, j_alt_r,sam_out))))
        #break
