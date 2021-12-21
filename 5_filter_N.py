#...
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

seqnames, seqs = readFasta("./final.fas")
count = 0
with open("./final_filtered_N_0.01.fas", 'w') as FILE:
    for i in range(len(seqs)):
        seqname = seqnames[i]
        seq = seqs[i]
        if seq.count("N") / len(seq) <= 0.01 or seqname.find("B.1.1.529") != -1:
            FILE.write(f"{seqname}\n")
            FILE.write(f"{seq}\n")
            count+=1
        else:
            print(f"Ignore {seqname}, {seq.count('N') / len(seq)} of N.")
            

print(f"{count} / {len(seqs)}, {count/len(seqs)*100}%")