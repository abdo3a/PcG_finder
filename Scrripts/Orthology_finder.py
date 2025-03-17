#!/usr/bin/env python3

import os, sys
import argparse
from Bio import SearchIO
from Bio import SeqIO
import re
import csv

parser = argparse.ArgumentParser(description='finding proteins orthologs')
parser.add_argument('-i', '--infile', help='.fasta file with the organism total pridicted proteins', required=True)
parser.add_argument('-q', '--quary', help='.fsa file with the querry protein sequence', required=True)
parser.add_argument('-dom', '--domain', help='conserved/catalatic protein domain', required=True)
parser.add_argument('-NOG', '--eggNOG', help='eggNOG protein cluster accession based eggNOG database', required=True)

parser.add_argument('-o', '--outfolder', help='outfolder', required= True)
parser.add_argument('-t', '--threads', help='no_threads', default= 10)


args = parser.parse_args()
infile = args.infile
output = args.outfolder
threads = args.threads
q_file = args.quary
dom_n = args.domain
NOG_n = args.eggNOG


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

#jack and parse results
def jack(quary):
    DB_l =[]
    jack_out = os.path.join('./', infile.replace('.fasta','.jack'))
    os.system('jackhmmer --tblout %s --cpu %s -N 10 --noali %s %s'%(jack_out, threads, q_file, new_input))
    for qresult in SearchIO.parse(jack_out, 'hmmer3-tab'):
        for hit in qresult.hits:
            DB_l.append(hit.id)
    os.remove(jack_out)
    with open(os.path.join('./', infile.replace('.fasta','.fa')), "w") as f:
        seqs_DBs = SeqIO.parse(open(new_input),'fasta')
        for seq in seqs_DBs:
            if seq.id in DB_l:
                SeqIO.write([seq], f, "fasta")
    f.close()

jack(q_file)
print ('jackhmmer done')


#emapper
for file_n in files('./'):
    if file_n.endswith(".fa"):
        if os.stat(file_n).st_size != 0:
            if not os.path.exists(os.path.join(out_path, 'annotations')):
                os.mkdir(os.path.join(out_path, 'annotations'))
            anno_path = os.path.join(out_path, 'annotations/')
            eg_in = os.path.join('./', file_n)       
            eg_out = os.path.join('./', file_n.replace('.fa','')) 
            seed_orthologs = os.path.join('./', file_n.replace('.fa','.emapper.seed_orthologs'))
            NOG_hits = os.path.join('./', file_n.replace('.fa','.emapper.hits'))
            os.system('emapper.py -m diamond --cpu %s -i %s -o %s'%(threads,eg_in,eg_out))
            os.remove(seed_orthologs)
            os.remove(NOG_hits)
            os.remove(eg_in)
        else:
            os.remove(os.path.join('./', file_n))
            os.remove(new_input)
            print('Orthologs NOT FOUND')
            sys.exit(1)


for file_n in files('./'):
    if file_n.endswith(".annotations"):
        with open(os.path.join('./', file_n)) as f:
            pattern = re.compile('.*'+NOG_n+'.*')
            kog_l = re.findall(pattern, f.read())
            os.rename(os.path.join('./', file_n), os.path.join(anno_path, file_n.replace('.emapper.annotations','.annotations')))
kog_select = []
for id in kog_l:
    kog_select.append(id.split('\t')[0])
with open(os.path.join('./', infile.replace('.fasta','.faa')), "w") as f:
    seqs_DBs = SeqIO.parse(open(new_input),'fasta')
    for seq in seqs_DBs:
        if seq.id in kog_select:
            SeqIO.write([seq], f, "fasta")
f.close()
print ('eggNOG done')


#domain screening
dom_D={}
for file_n in files('./'):
    if file_n.endswith(".faa"):
        if os.stat(file_n).st_size != 0:
            dom_in = os.path.join('./', file_n)       
            dom_out = os.path.join('./', file_n.replace('.faa','.domhits')) 
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
                    print('Orthologs NOT FOUND')
                    sys.exit(1)          
            os.remove(dom_out)
        else:
            os.remove(os.path.join('./', file_n))
            os.remove(new_input)
            print('Orthologs NOT FOUND')
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
        seqs_DBs = SeqIO.parse(open(new_input),'fasta')
        for seq in seqs_DBs:
            if seq.id in dic.keys():
                SeqIO.write([seq], f, "fasta")
    f.close()
    os.rename(os.path.join('./', infile.replace('.fasta','.faa')), os.path.join(out_path, infile.replace('.fasta','.faa')))
    os.remove(new_input)

result(dom_n)
