#cat file and remove same seq
import re, os, sys
def readFasta(path):
    count = 0
    seqnames, seqs = [], []
    Loading = ["\\","\\","|", "|","/", "/", "-", "-"]
    with open(path) as FILE:
        seq = ""
        for line in FILE:
            line = line.strip()
            if line.find(">") != -1:
                if seq != "":
                    seqnames.append(seqname)
                    seqs.append(seq)
                    seq = ""
                seqname = line
                count += 1
                print(f"[{Loading[count % 8]}]", end = "\r")
                if count % 100 == 0 or count == 1:
                    print("   Loading",count, f"sequences from {path}", end = "\r")
            else:
                seq += line
    print("   Loading",count, f"sequences from {path}", end = "\r", flush = True)
    print(f"[*]", end = "\r", flush = True)
    seqnames.append(seqname)
    seqs.append(seq)
    print()
    return seqnames, seqs
def writeNewFasta(path, seqnames, seqs):
    with open(path, 'w') as FILE:
        for i in range(len(seqnames)):
            FILE.write(seqnames[i] + "\n")
            FILE.write(seqs[i] + "\n")

            print(f"   Write to file '{path}' ... ({i+1}/{len(seqnames)}, {(i+1)/len(seqnames)*100:.2f}%)",end = "\r", flush=True)
    print()

def Main():
    folder = "./data/"
    files = os.listdir(folder)
    files.sort()
    seqnames, seqs = [], []
    seqnames_dic = {}
    same_c = 0
    for file in files:
        tempseqnames, tempseqs = readFasta(folder + file)
        for i in range(len(tempseqnames)):
            seqname, seq = tempseqnames[i], tempseqs[i]
            if seqname not in seqnames_dic:
                seqnames.append(seqname)
                seqs.append(seq)
                seqnames_dic.setdefault(seqname, file)
            else:
                same_c += 1
                print(f"[{same_c:>2}]\t'{file}':\t{seqname}\twas found in\t'{seqnames_dic[seqname]}'")
    if len(seqnames) == len(seqs):
        print(f"Sequences number = {len(seqnames)}")
        writeNewFasta(f"./allData_n{len(seqnames)}.fas", seqnames, seqs)
    else:
        sys.exit("Something get wrong.")

Main()
