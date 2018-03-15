#use bed file to represent breaks in GFA file, mainly used to show breaks in gEval
import sys
if __name__ == "__main__":
    for ln in sys.stdin:
        lnlist = ln.strip('\n').split('\t')
        if lnlist[0] == 'S' and len(lnlist) > 5:
            i = 4
            while i < len(lnlist):
                if (lnlist[i][0:2] == "BK"):
                    bk_list = lnlist[i].split(':')
                    point_list = bk_list[2].split(',')
                    print lnlist[1]+"\t"+point_list[0][1:] + "\t" + point_list[1]
                i += 1


