#!/usr/bin/env python3

import os, sys
import argparse
from Bio import SearchIO
from Bio import SeqIO
from Bio import AlignIO
import re
import csv

parser = argparse.ArgumentParser(description='finding histone H3 variants(H3_1&H3_3) homologous')
parser.add_argument('-i', '--infile', help='.fasta file with the organism total pridicted proteins', required=True)
parser.add_argument('-hv', '--Histone', help='histone H3 variants(H3_1&H3_3)', required= True)
parser.add_argument('-o', '--outfolder', help='outfolder', required= True)
parser.add_argument('-t', '--threads', help='no_threads', default= 10)


args = parser.parse_args()
infile = args.infile

hv_name = args.Histone
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

#rename duplcate sequences
records = set()
of = open(os.path.join('./', infile.replace('.fasta','ed.fasta')), "w")
for record in SeqIO.parse(infile, "fasta"):
    ID = record.id
    num = 1
    while ID in records:
        ID = "{}_{}".format(record.id, num)
        num += 1
    records.add(ID)
    record.id = ID
    record.name = ID
    record.description = ID
    SeqIO.write(record, of, "fasta")
of.close()

new_input = os.path.join('./', infile.replace('.fasta','ed.fasta'))

#hmmrsearch and parse results
def hmmr(hv):
    H3_1_l=[]
    H3_3_l=[]
    hmmr_q = os.path.join('./DBs/queries/histone_H3.fasta')
    hmmr_out = os.path.join('./', infile.replace('.fasta','_'+hv+'.sto'))
    os.system('jackhmmer -A %s --cpu %s -N 10 %s %s'%(hmmr_out, threads, hmmr_q, new_input))
    for alignment in AlignIO.parse(hmmr_out, 'stockholm'):
        for alg in alignment:
           id = alg.id
           id = id.split('/')[0]
           if alg.id != 'AT1G09200.1':
               if 'APA' in alg.seq and 'QSSAVA' in alg.seq and 'ARKST' in alg.seq:
                   H3_1_l.append(id)
               elif 'APT' in alg.seq and 'QSHAVL' in alg.seq and 'STG' in alg.seq and 'KQL' in alg.seq:
                   H3_3_l.append(id)
    os.remove(hmmr_out)
    
    with open(os.path.join('./', infile.replace('.fasta','.fa')), "w") as f:
        seqs_DBs = SeqIO.parse(open(new_input),'fasta')
        for seq in seqs_DBs:
            if hv == 'H3_1' and seq.id in H3_1_l:
                SeqIO.write([seq], f, "fasta")
            elif hv == 'H3_3' and seq.id in H3_3_l:
                SeqIO.write([seq], f, "fasta")
    f.close()

hmmr(hv_name)

print ('hmmsearch done')



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
                    os.remove(new_input)
                    print(hv_name +' homologous NOT FOUND')
                    sys.exit(1)          
            os.remove(dom_out)
        else:
            os.remove(os.path.join('./', file_n))
            os.remove(new_input)
            print(hv_name +' homologous NOT FOUND')
            sys.exit(1)


#remove false-positives
def result(hv,dom):
    pattern = re.compile('.*'+dom+'.*')
    dic = {key:val for key, val in dom_D.items() if any(pattern.match(line) for line in val)}
    out_repo = csv.writer(open(os.path.join('./', infile.replace('.fasta','_'+hv+'.csv')), "a"))    
    out_repo.writerow(['Proteins','Domains'])
    for key, value in dic.items():
        out_repo.writerow([key, '|'.join(value)])
    os.rename(os.path.join('./', infile.replace('.fasta','_'+hv+'.csv')), os.path.join(out_path, infile.replace('.fasta','_'+hv+'.csv')))
    with open(os.path.join('./', infile.replace('.fasta','_' + hv +'.faa')), "w") as f:
        seqs_DBs = SeqIO.parse(open(new_input),'fasta')
        for seq in seqs_DBs:
            if seq.id in dic.keys():
                SeqIO.write([seq], f, "fasta")
    f.close()
    os.rename(os.path.join('./', infile.replace('.fasta','_' + hv +'.faa')), os.path.join(out_path, infile.replace('.fasta','_' + hv +'.faa')))
    os.remove(new_input)


result(hv_name,'Pfam_Histone')
