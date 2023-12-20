## this script is to find read that have CCCATGTACT  and remain CCC  at 3'
## for example: "GAGTGAAGCCATTCAAAAGACAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGGATATCCGACTCTGCGTTGATACCACCTGTCCCTTATACACAACTCCGAGCCCACCAGACTCCTTAGCATCTCGTATGCCGTCTTCTGCTTGAAA" will change into "GAGTGAAGCCATTCAAAAGACAAAAA"
## this file is to find the last poly A. If want to find the first poly A, then [m.start() for m in re.finditer("A{5,150}", seq)][0].  If want to find the last poly A, then [m.start() for m in re.finditer("A{5,150}", seq)][-1]
## If want to find A20, change A{10,30} to A{20,30}, A{10,150} to A{20,150}
## If want to remain 5 index =  [m.start() for m in re.finditer("A{10,150}", seq)][0] + 5;  remain 10   index =  [m.start() for m in re.finditer("A{10,150}", seq)][0] + 10
import os
import gzip, argparse, re

parser = argparse.ArgumentParser(usage="it's usage tip.", description="help info.")
parser.add_argument("-i", "--input", help="the port number.")
parser.add_argument("-o", "--out", help="the file type")

args = parser.parse_args()
i = args.input
o = args.out
fq_i = open(i)
fq_o = open(o, "w")

for line in fq_i:
        seq_name = line.strip("\n")
        seq = next(fq_i).strip("\n")
        seq_inf = next(fq_i).strip("\n")
        seq_qua = next(fq_i).strip("\n")
        if re.findall("AACGCAGAGT[A|G|C|T][A|G|C|T][A|G|C|T][A|G|C|T][A|G|C|T][A|G|C|T][A|G|C|T][A|G|C|T][A|G|C|T][A|G|C|T][A|G|C|T][A|G|C|T]AGTTAA",seq):
                index =  [m.start() for m in re.finditer("AACGCAGAGT[A|G|C|T][A|G|C|T][A|G|C|T][A|G|C|T][A|G|C|T][A|G|C|T][A|G|C|T][A|G|C|T][A|G|C|T][A|G|C|T][A|G|C|T][A|G|C|T]AGTTAA", seq)][0] 
                seq = seq[index:]
                seq_qua = seq_qua[index:]
                if index != 0:
                        fq_o.write(seq_name + "\n" + seq + "\n" + seq_inf + "\n" + seq_qua + "\n")


