# tag the lineage of each sequence
import os, re

def createFolder(path):
    if os.path.isdir(path):
        os.system(f"rm -r {path}")
        os.mkdir(path)
        print(f"Clear {path}")
    else:
        os.mkdir(path)
        print(f"mkdir {path}")

def readReport(path):
    dic = {}
    with open(path) as FILE:
        for line in FILE:
            if line.find("taxon") == -1:
                line = line.strip().split(",")
                seqname = line[0]
                lineage = line[1]
                strain = line[4]
                if strain != "":
                    strain = "_" + strain.split(" ")[0]
                dic.setdefault(seqname, f"|{lineage}{strain}")
    return dic

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

folder = "./align/"
files = os.listdir(folder)
for i in range(len(files)):
    if files[i].find("file") == -1:
        files.pop(i)
files.sort()
report_folder = "./report/"
createFolder(report_folder)
strain_dic = {}
count = 1
for file in files:
    if file.find("file") != -1:
        print(f"=> {file} ({count} / {len(files)}, {count/len(files)*100:.2f}%)")
        os.system(f"bash -i ini_pangolin.sh {folder + file} {report_folder}")
        readReport(report_folder + "lineage_report.csv")
        num = int(re.findall(r"(\d+)", file)[0])
        strain_dic.update(readReport(f"{report_folder}lineage_report.csv"))
        os.system(f"mv {report_folder + 'lineage_report.csv'} {report_folder + f'lineage_report_{num}.csv'}")

        count += 1
        print(f"=> {len(strain_dic)}")


seqnames, seqs = readFasta("./final.fas")
print("\nTagging Strain on seqname ...")
with open("./final_taglineage.fas", 'w') as FILE:
    for i in range(len(seqnames)):
        seqname = seqnames[i]
        seq = seqs[i]

        newseqname = seqname + strain_dic[seqname[1:]]
        FILE.write(newseqname + "\n")
        FILE.write(seq + "\n")
        print(f"{i+1} / {len(seqnames)}, {(i+1)/len(seqnames)*100:.2f}%", end = "\r", flush=True)
print()
os.system("rm ./final.fas")
os.system("mv ./final_taglineage.fas ./final.fas")