#!/usr/bin/env python3

import os, sys
import argparse
from Bio import SearchIO
from Bio import SeqIO
import re
import csv

parser = argparse.ArgumentParser(description='finding domain containing proteins')
parser.add_argument('-i', '--infile', help='.fasta file with the organism total pridicted proteins', required=True)
parser.add_argument('-o', '--outfolder', help='outfolder', required= True)
parser.add_argument('-t', '--threads', help='no_threads', default= 10)
parser.add_argument('-hmm', '--hmm', help='hmm model file', required= True)
parser.add_argument('-dom', '--dom', help='hmm model file', required= True)

args = parser.parse_args()
infile = args.infile
dom_hmm = args.hmm
dom_n = args.dom 
output = args.outfolder
threads = args.threads

#read folder files function
def files(path):
    for file in os.listdir(path):
        if os.path.isfile(os.path.join(path, file)):
            yield file

if not os.path.exists(output):
    os.mkdir(os.path.join('./', output))

out_path = os.path.join('./', (output +'/'))
db_path = './DBs/'


#VEFS_domain search and parse results

DB_l =[]
dom_in = os.path.join('./', infile)       
dom_out = os.path.join('./', infile.replace('.fasta','.domhits'))

#you need to download, build and press the aligment first using:
# hmmbuild -n SET PF00856.hmm PF00856_full.txt
# hmmpress PF00856.hmm


os.system('hmmsearch --tblout %s --noali --cpu %s %s %s'%(dom_out,threads,dom_hmm,dom_in))
for qresult in SearchIO.parse(dom_out, 'hmmer3-tab'):
    for hit in qresult.hits:
        DB_l.append(hit.id)
os.remove(dom_out)
with open(os.path.join('./', infile.replace('.fasta','.fa')), "w") as f:
    seqs_DBs = SeqIO.parse(open(infile),'fasta')
    for seq in seqs_DBs:
        if seq.id in DB_l:
            SeqIO.write([seq], f, "fasta")
f.close()

#domain screening

dom_D={}
for file_n in files('./'):
    if file_n.endswith(".fa"):
        if os.stat(file_n).st_size != 0:
            dom_in = os.path.join('./', file_n)       
            dom_out = os.path.join('./', file_n.replace('.fa','.domhits')) 
            os.system('hmmscan --tblout %s --noali --cpu %s %s %s'%(dom_out,threads,os.path.join(db_path, 'Pf_Sm'),dom_in))
            os.remove(dom_in)
            for qresult in SearchIO.parse(dom_out, 'hmmer3-tab'):
                query_id = qresult.id
                hits = qresult.hits
                num_hits = len(hits)
                if num_hits > 0:
                    for i in range(0,num_hits):
                        if query_id in dom_D:
                            dom_D[query_id].append(hits[i].id)
                        else:
                            dom_D[query_id]=[hits[i].id]
                else:
                    os.remove(os.path.join('./', file_n))
                    print('VEFS_domain homologous NOT FOUND')
                    sys.exit(1)          
            os.remove(dom_out)
        else:
            os.remove(os.path.join('./', file_n))
            print('VEFS_domain homologous NOT FOUND')
            sys.exit(1)

#remove false-positives
def result(dom):
    pattern = re.compile('.*'+dom+'.*')
    dic = {key:val for key, val in dom_D.items() if any(pattern.match(line) for line in val)}
    out_repo = csv.writer(open(os.path.join('./', infile.replace('.fasta','.csv')), "a"))    
    out_repo.writerow(['Proteins','Domains'])
    for key, value in dic.items():
        out_repo.writerow([key, '|'.join(value)])
    os.rename(os.path.join('./', infile.replace('.fasta','.csv')), os.path.join(out_path, infile.replace('.fasta','.csv')))
    with open(os.path.join('./', infile.replace('.fasta','.faa')), "w") as f:
        seqs_DBs = SeqIO.parse(open(infile),'fasta')
        for seq in seqs_DBs:
            if seq.id in dic.keys():
                SeqIO.write([seq], f, "fasta")
    f.close()
    os.rename(os.path.join('./', infile.replace('.fasta','.faa')), os.path.join(out_path, infile.replace('.fasta','.faa')))


result(dom_n)