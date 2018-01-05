import sys

if __name__ == "__main__":

    id_ctr = ""
    print "\t".join(["H", "VZ:i:2.0"])
    with open(sys.argv[1],"r") as ifi:
        for line in ifi:
            if not line.startswith(">"):
                print "\t".join(["S", id_ctr, str(len(line.strip())), line.strip(), ""])
            else:
                id_ctr = line[1:].strip().split(' ')[0]

