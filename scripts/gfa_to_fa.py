#! env python2

import sys

if __name__ == "__main__":
    if len(sys.argv) < 1:
        print "No GFA file specified, exit 1"
    else:
        fl = open(sys.argv[1])
        #first line get version
        frst_ln_list = fl.readline().strip('\n').split('\t')
        vn = -1
        for i in frst_ln_list:
            if i[0:2] == "VN":
                vn = float(i[5:])
        if vn >= 1.0 and vn < 2.0:
            for ln in fl:
                ln_list = ln.strip('\n').split('\t')
                if ln_list[0] == 'S':
                    print ">"+ln_list[1]
                    print ln_list[2]
        elif vn >= 2.0:
            for ln in fl:
                ln_list = ln.strip('\n').split('\t')
                if ln_list[0] == 'S':
                    print ">"+ln_list[1]
                    print ln_list[3]
        else:
            print "No version info found in GFA file, quit converting..."
