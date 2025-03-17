#!/usr/bin/env python3

import os, sys
import argparse
from Bio import SearchIO
from Bio import SeqIO

import re
import csv

parser = argparse.ArgumentParser(description='finding histone H3 variants(H3.1&H3.3) homologous')
parser.add_argument('-i', '--infile', help='.fasta file with the organism total pridicted proteins', required=True)
parser.add_argument('-o', '--outfolder', help='outfolder', required= True)
parser.add_argument('-t', '--threads', help='no_threads', default= 10)


args = parser.parse_args()
infile = args.infile

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
def jack(quary,hv):
    DB_l =[]
    hv_q = os.path.join('./DBs/queries/',quary)
    hv_out = os.path.join('./', infile.replace('.fasta','_'+hv+'.jack'))
    os.system('jackhmmer --tblout %s --cpu %s -N 10 --noali %s %s'%(hv_out, threads, hv_q, new_input))
    for qresult in SearchIO.parse(hv_out, 'hmmer3-tab'):
        for hit in qresult.hits:
            DB_l.append(hit.id)
    os.remove(hv_out)
    with open(os.path.join('./', infile.replace('.fasta','.fa')), "a") as f:
        seqs_DBs = SeqIO.parse(open(new_input),'fasta')
        for seq in seqs_DBs:
            if seq.id in DB_l:
                seq.description = seq.id + "_" + hv  
                seq.id = seq.description
                SeqIO.write([seq], f, "fasta")
    f.close()    

jack('histone_H3_1.fasta','H3_1')
jack('histone_H3_3.fasta','H3_3')
print ('hmmsearch done')

#de-duplicate 
for file_n in files('./'):
    if file_n.endswith(".fa"):
        if os.stat(file_n).st_size != 0:
            with open(os.path.join('./', infile.replace('.fasta','.var')), 'a') as outFile:
                record_ids = list()
                for record in SeqIO.parse(file_n, 'fasta'):
                    if record.id not in record_ids:
                        record_ids.append(record.id)
                        SeqIO.write(record, outFile, 'fasta')

os.remove(os.path.join('./', infile.replace('.fasta','.fa')))
#domain screening
dom_D={}
for file_n in files('./'):
    if file_n.endswith(".var"):
        if os.stat(file_n).st_size != 0:
            dom_in = os.path.join('./', file_n)       
            dom_out = os.path.join('./', file_n.replace('.var','.domhits')) 
            os.system('hmmscan --tblout %s --noali --cpu %s %s %s'%(dom_out,threads,os.path.join(db_path, 'Pf_Sm'),dom_in))
            for qresult in SearchIO.parse(dom_out, 'hmmer3-tab'):
                query_id = qresult.id
                query_id = query_id.replace('_H3_1','').replace('_H3_3','')
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
                    print('H3 homologous NOT FOUND')
                    sys.exit(1)          
            os.remove(dom_out)
            os.rename(os.path.join('./', file_n), os.path.join(out_path, file_n))
        else:
            os.remove(os.path.join('./', file_n))
            os.remove(new_input)
            print('H3 homologous NOT FOUND')
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
    os.rename(os.path.join('./', infile.replace('.fasta','_' + hv +'.faa')), os.path.join(out_path, infile.replace('.fasta','_' + hv +'.faa')))
    f.close()     
    os.remove(new_input)


result('H3_1_3','Pfam_Histone')