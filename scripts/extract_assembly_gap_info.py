#!/usr/bin/python



import sys

class aln_unit:
    q_id = ""
    p_infor = {}
    def add_cords(self, p_id, cord):
        if pid not in p_infor:
            self.p_infor[pid] = []
        self.p_infor[pid].append(cord)
    def cal_cl(self):
        for e in self.p_infor:
            self.p_infor[e] = merge(self.p_infor[e]) 
    def merge(self, cors):
        a = cors.sort(lambda x:x[0]) 
        s = a[0] 
        z = []
        length = 0
        for st, en in a:
            if st <= s[1]:
                s[1] = max(s[1], en)
            else:
                z.append(s)
                length += s[1] -s[0]
                s[0] = st
                s[1] = en
        z.append(s)
        length += s[1] -s[0]
        z.append(length)
        return z
class cord:
    q_x = q_y = 0 #p_x = p_y = 0
    orient = True
    def _init_(self, q_x, q_y, orient):
        self.q_x = q_x
        self.q_y = q_y
        # self.p_x = p_x
        # self.p_y = p_y
        self.orient = orient
    
def process_paf(flname):
    fl = open(flname, "r")
    lnlist = fl.readline().strip().split('\t')
    pre_id = lnlist[5]
    aln_units = []
    cors = [[int(lnlist[7]), int(lnlist[8])]]
    len_s = lnlist[6]
    for ln in fl:
        lnlist = ln.strip().split('\t')
        r_id = lnlist[5]
        if r_id == pre_id:
            cors.append([int(lnlist[7]),int(lnlist[8])])
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
    

