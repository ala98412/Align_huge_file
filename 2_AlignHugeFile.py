import numpy as np
import sys, math, os, datetime, re, time


filename = "./allData_n2.fas"
seq_per_file = 10

def minDate(mindate, date2):
    f1 = time.strptime(mindate, "%Y-%m-%d")
    f2 = time.strptime(date2, "%Y-%m-%d")
    if f2 < f1:
        return date2
    else:
        return mindate
def maxDate(maxdate, date2):
    f1 = time.strptime(maxdate, "%Y-%m-%d")
    f2 = time.strptime(date2, "%Y-%m-%d")
    if f2 > f1:
        return date2
    else:
        return maxdate
def readLineageReport(path):
    n_count, alpha_count = 0, 0
    min_date, max_date = "2222-12-31", "1911-12-31"
    with open(path) as FILE:
        for line in FILE:
            line = line.split(",")
            name = line[0]
            if name != "taxon":
                n_count += 1
                strain = line[4].split(" ")[0]
                if strain == "Alpha":
                    alpha_count += 1
                    date = re.findall(r"\|(\d+\-\d+\-\d+)", name)
                    if date:
                        date = date[0]
                        min_date = minDate(min_date, date)
                        max_date = maxDate(max_date, date)
  
    return min_date, max_date, n_count, alpha_count
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
def findDel(refseq):
    del_sites = []
    start = -1
    for i in range(len(refseq)):
        nt = refseq[i]
        if nt == '-':
            if start == -1:
                start = i
        else:
            if start != -1:
                del_sites.append([start, i])
                start = -1   
    if start != -1:
        del_sites.append([start, len(refseq)])

    cutting_sites = []
    start = 0
    for del_site in del_sites:
        end = del_site[0]
        cutting_sites.append([start, end])
        start = del_site[1]
  
    if start != len(refseq):
        cutting_sites.append([start, len(refseq)])
  
    return cutting_sites
def cutDelRegion(cutting_sites, seqs, seqnames):
    new_seqs = []
    for i in range(len(seqs)):
        seq = seqs[i]
        newseq = ""
        for cut in cutting_sites:
            newseq += seq[cut[0]:cut[1]]
        #print(seqnames[i], len(newseq), newseq.count("-"))
        new_seqs.append(newseq)
  
    return new_seqs
def writeNewFasta(path, seqnames, seqs):
    with open(path, 'w') as FILE:
        for i in range(len(seqnames)):
            FILE.write(seqnames[i] + "\n")
            FILE.write(seqs[i] + "\n")

            print(f"   Write to file '{path}' ... ({i+1}/{len(seqnames)}, {(i+1)/len(seqnames)*100:.2f}%)",end = "\r", flush=True)
    print()

def createFolder(path):
    if os.path.isdir(path):
        os.system(f"rm -r {path}")
        os.mkdir(path)
        print(f"Clear {path}")
    else:
        os.mkdir(path)
        print(f"mkdir {path}")
        
seqnames, seqs = readFasta(filename)

if len(seqnames) == len(seqs):
    print(f"\nFound {len(seqnames)} sequences. ({len(seqnames)} = {len(seqs)})")
else:
    sys.exit("Fasta format Error.")

len_total = 0
for i in range(len(seqs)):
    len_total += len(seqs[i])

print("Avg len =",len_total/len(seqs))


# >>> divide
divide = round(len(seqs)/seq_per_file)
if divide == 0:
    divide = 1
print(f"Divide into {divide} files.")

f = np.linspace(0, len(seqs), divide+1)
print(f)
new_f = []
for i in f:
    new_f.append(math.floor(i))

createFolder("./split/")
print(new_f)

for i in range(len(new_f)-1):
    with open(f"./split/file{i+1}.fas", 'w') as FILE:
        if i != 0:
            FILE.write(seqnames[0]+"\n")
            FILE.write(seqs[0]+"\n")
        for k in range(new_f[i], new_f[i+1]):
            FILE.write(seqnames[k]+"\n")
            FILE.write(seqs[k] +"\n")
            print(f"=> file{i+1} / {divide}: {k+1} / {len(seqs)} {(k+1)/len(seqs)*100:.2f}%", end = "\r")
print()
# <<< divide

#Align
createFolder("./align/")
for i in range(divide):
    print(f"=> file {i+1} / {divide} ({(i+1)/divide*100:.2f}%):")

    print(f"   Aligning ... ", end = "", flush=True)
    os.system(f"mafft --thread -1 --reorder --anysymbol --auto --quiet ./split/file{i+1}.fas > ./align/file{i+1}_mafft.fas")
    print(f"Done", flush = True)

    seqnames, seqs = readFasta(f"./align/file{i+1}_mafft.fas")
    print(f"\n   Checking alignment ... ", end= "", flush = True)
    if len(seqs[0]) != 29903:
        print(" Fail", flush = True)
        # sys.exit(f"file{i+1}.fas didn't pass QC.")
        print(f"   Cut off DEL region in Ref seq ... ", flush=True)
        cutting_sites = findDel(seqs[0])
        seqs = cutDelRegion(cutting_sites, seqs, seqnames)
        writeNewFasta(f"./align/file{i+1}_mafft.fas", seqnames, seqs)
    else:
        print(f" Pass [len(seq[0]) == {len(seqs[0])}]", flush = True)
        writeNewFasta(f"./align/file{i+1}_mafft.fas", seqnames, seqs)
    if i == 0:
        writeNewFasta(f"./align/ref.fas", [seqnames[0]], [seqs[0]])

# merge all file
os.system("awk 'FNR>2' ./align/file* > f.fas")
os.system("cat ./align/ref.fas f.fas > final.fas")
os.system("rm f.fas")
