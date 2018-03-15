#!env python
# Extract scaffold boundary and write into a bed file
import sys
def readfq(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break
def split_contigs_by_Ns(seq):
    p = []
    seq_len = len(seq)
    i = 0
    while i < seq_len:
        if seq[i] != 'N' and seq[i] != 'n':
            s = i
            i += 1
            while i < seq_len and seq[i] != 'N' and seq[i] != 'n':
                i += 1
        
            e = i 
            p.append([s, e])
        i += 1
    return p


def extract_Ns_cords(seq):
    p = []
    seq_len = len(seq)
    i = 0
    while i < seq_len:
        if seq[i] == 'N' or seq[i] == 'n':
            s = i
            i += 1
            while i < seq_len and (seq[i] == 'N' or seq[i] == 'n'):
                i += 1
        
            e = i - 1
            p.append([s, e])
        i += 1
    return p
import argparse
if __name__ == "__main__":
    # parser = argparse.ArgumentParser(description="Process Ns in reference or query file")
    # parser.add_argument('--'
    if sys.argv[1] == "r":
        print "ref\tref_start\tref_end\tname\tstrand"
        for name, seq, qual in readfq(sys.stdin):
            p = extract_Ns_cords(seq)
            for e in p:
                print name+"\t"+str(e[0])+"\t"+str(e[1])+"\tNone\t+"
    elif sys.argv[1] == 'q':
        print "query\tquery_start\tquery_end\tname\tstrand"
        for name, seq, qual in readfq(sys.stdin):
            p = extract_Ns_cords(seq)
            for e in p:
                print name+"\t"+str(e[0])+"\t"+str(e[1])+"\tNone\t+"
    else:
        for name, seq, qual in readfq(sys.stdin):
            p = split_contigs_by_Ns(seq)
            i = 0
            for e in p:
                print ">"+name+"_dg_s"+str(i)
                print seq[e[0]:e[1]]
                i += 1
