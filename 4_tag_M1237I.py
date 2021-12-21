#tag M1237I
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

    return seqnames, seqs
def writeNewFasta(path, seqnames, seqs):
    with open(path, 'w') as FILE:
        for i in range(len(seqnames)):
            FILE.write(seqnames[i] + "\n")
            FILE.write(seqs[i] + "\n")

            print(f"   Write to file '{path}' ... ({i+1}/{len(seqnames)}, {(i+1)/len(seqnames)*100:.2f}%)",end = "\r", flush=True)
    print()

M1237I = 25273

seqnames, seqs =  readFasta("./final.fas")

count = 0
count_C, count_T = 0, 0
for i in range(len(seqs)):
    seq = seqs[i]
    seqnames[i] = seqnames[i].replace("M1237I.", "")
    if seq[M1237I-1] == 'C':
        seqnames[i] =  ">M1237I.C." + seqnames[i][1:].replace("hCoV-19/", "")
        count += 1
        count_C += 1
    elif seq[M1237I-1] == 'T':
        seqnames[i] =  ">M1237I.T." + seqnames[i][1:].replace("hCoV-19/", "")
        count += 1
        count_T += 1
    
        # print(f"[{count}] {seqnames[i]}")
print(f"\n   {count} sequences with M1237I.")
print(f"G->C: {count_C}")
print(f"G->T: {count_T}")
writeNewFasta("./final.fas", seqnames, seqs)