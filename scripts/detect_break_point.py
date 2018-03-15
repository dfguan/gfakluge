#!/usr/bin/python

#check break point of dnaNexus gfa according to primary contigs in dirk's assembly 
# Input: <PAF> (alignment result using minimap2)
# Output: rule: <CONTIG_ID>\t<CONTIG_LEN>\t<# of Gaps>\t<GAPPED_BASES>\t<COVERED_INTERVALS>\n
# Example: python detech_break_point.py ALN.paf

import sys
def merge(cors):
    a = sorted([sorted(t) for t in cors]) 
    s = list(a[0]) 
    for st, en in a:
        if st <= s[1]:
            s[1] = max(s[1], en)
        else:
            yield tuple(s)
            s[0] = st
            s[1] = en
    yield tuple(s)

def process_paf(flname):
    fl = open(flname, "r")
    lnlist = fl.readline().strip().split('\t')
    pre_id = lnlist[5]
    cors = [[int(lnlist[7]), int(lnlist[8])]]
    len_s = lnlist[6]
    for ln in fl:
        lnlist = ln.strip().split('\t')
        r_id = lnlist[5]
        if r_id == pre_id:
            cors.append([int(lnlist[7]),int(lnlist[8])])
        else:
            pre_cor = 0
            tup_list = list(merge(cors)) 
            gap_n = 0;
            gap_bases = 0
            for e in tup_list:
                # print e[0], pre_cor
                z = e[0] - pre_cor
                if z > 100:
                    gap_n += 1
                    gap_bases += z
                pre_cor = e[1]
            z = int(len_s) - pre_cor
            if z > 100:
                gap_n += 1
                gap_bases += z
            print pre_id, len_s, str(gap_n), str(gap_bases), tup_list 
            cors = [[int(lnlist[7]), int(lnlist[8])]]
            pre_id = lnlist[5]
            len_s = lnlist[6]
    pre_cor = 0
    tup_list = list(merge(cors)) 
    gap_n = 0;
    gap_bases = 0
    for e in tup_list:
        z = e[0] - pre_cor
        if z > 100:
            gap_n += 1
            gap_bases += z
        pre_cor = e[1]
    z = int(len_s) - pre_cor
    if z > 100:
        gap_n += 1
        gap_bases += z
    print pre_id, len_s, str(gap_n), str(gap_bases), tup_list 


if __name__ == "__main__":
    flname = sys.argv[1]
    process_paf(flname)    
    

