#! env python2

import sys

if __name__ == "__main__":
    if len(sys.argv) < 1:
        print "No GFA file specified, exit 1"
    else:
        fl = open(sys.argv[1])
        for ln in fl:
            ln_list = ln.strip('\n').split('\t')
            if ln_list[0] == 'S':
                print ">"+ln_list[1]
                print ln_list[2]

