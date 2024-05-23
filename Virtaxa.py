import pandas as pd
import math
from tqdm import tqdm
import warnings
warnings.filterwarnings("ignore")

import sys
import os
import subprocess
import argparse

import copy
import time
from Bio import SeqIO

import Virtaxa_fragment as func_fragment

parser = argparse.ArgumentParser(description='parameters of virtaxa', formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('--query', '-q', help='proteins of query contigs')
parser.add_argument('--db', '-d', help='reference database folder')
parser.add_argument('--output', '-o', help='output folder')
args = parser.parse_args()

test_faa = args.query
db_folder = args.db
result_folder = args.output


test_csv = test_faa.split('.faa')[0] + '.csv'
test_fasta = test_faa.split('.faa')[0] + '.fasta'
train2test_file = '%s/blast_train2test.tsv'%result_folder
singleton2test_file = '%s/blast_singleton2test.tsv'%result_folder

ref_fasta_file = '%s/train_fasta.csv'%db_folder
csv_faa_file = '%s/train_faa.csv'%db_folder

def time_stamp():
    now = int(round(time.time()*1000))
    time_stamp_1 = time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(now/1000))
    
    return time_stamp_1

def if_exist_folder(folder):
    if os.path.exists('%s'%folder)==False:
        os.mkdir('%s'%folder)
    else:
        print('\n---------Folder %s exists-----------\n'%folder)

def if_exist_file(file,cmd):
    if os.path.exists('%s'%file)==False:
        print(cmd)
        res = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        res.communicate()
        end = time.time()
        print('time: ' + time_stamp())
    else:
        print('\n---------File %s exists-----------\n'%file)

def test_faa_csv(result_folder):
    test_faa = SeqIO.parse('%s/tmp/test_1.faa'%result_folder,'fasta')
    tb = pd.DataFrame()
    ind=0
    for rec in tqdm(test_faa):
        new_row = {'Accession':rec.id,'fasta':str(rec.seq)}
        new1 = pd.DataFrame(new_row,index=[ind])
        tb = pd.concat([tb,new1])
        ind = ind + 1 
    return tb

def cmd_run(cmd):
    print(cmd)
    _ = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    _.communicate()

def make_fasta_csv(result_folder,ref_fasta_file):
    df_ref_fasta = pd.read_csv(ref_fasta_file)[['Accession','fasta']]
    test_fasta = SeqIO.parse('%s/test.fasta'%result_folder,'fasta')
    df_test = pd.read_csv('%s/test.csv'%result_folder)
    
    for record in test_fasta:
        df_test.loc[df_test['Accession']==record.id,'fasta']=str(record.seq)
    
    df_test_ref_fasta = pd.concat([df_ref_fasta,df_test])
    
    return df_test_ref_fasta

begin_all = time.time()
print('begin_time:' + time_stamp())

#主程序开始
    
if_exist_folder(result_folder)
if_exist_folder('%s/tmp'%result_folder)

cmd_run("cp %s %s/test.faa"%(test_faa,result_folder))
cmd_run("cp %s %s/test.fasta"%(test_fasta,result_folder))

df_test = pd.DataFrame()
record = SeqIO.parse(test_fasta,'fasta')
test_list = []
for rec in record:
    test_list.append(rec.id)
df_test['Accession']=test_list
df_test.to_csv('%s/test.csv'%result_folder,index=False)


df_test_ref_fasta = make_fasta_csv(result_folder,ref_fasta_file)
df_test_ref_fasta.to_csv('%s/tmp/test_ref_fasta.csv'%result_folder,index=False)
    

#开始 strict diamond prediction

cmd1 ='diamond makedb -p 52 --in %s/train.faa --db %s/train.dmnd'%(db_folder,db_folder)

cmd2 = 'diamond blastp --db %s/train.dmnd --query %s/test.faa --evalue 10 --out %s/tmp/diamond_blast_train2test.tsv --outfmt 6 qseqid sseqid pident qcovhsp length mismatch gapopen qstart qend sstart send evalue bitscore' % (db_folder,result_folder,result_folder)


if os.path.exists('%s/tmp/diamond_blast_train2test.tsv'%result_folder)==False:
    cmd_run(cmd1)
    cmd_run(cmd2)


if os.path.exists('%s/tmp/test_1.faa'%result_folder)==False:
    func_fragment.diamond_blast(db_folder,result_folder)


#开始 hmm prediction
if os.path.exists('%s/tmp/test_1_faa.csv'%result_folder)==False:
    test_faa_csv(result_folder).to_csv('%s/tmp/test_1_faa.csv'%result_folder,index=False)

if os.path.exists('%s/tmp/hmmsearch_result_domain.abc'%result_folder)==False:
    begin = time.time()
    func_fragment.hmm_search(db_folder,result_folder)
    end = time.time()
    print('time: ' + str(math.floor((end-begin)/60))+ ' min ' + str(math.floor((end-begin)%60)) + ' s \n')
else:
    print('\n---------File %s/hmmsearch_result.out exists-----------\n'%result_folder)

    
if os.path.exists('%s/tmp/pred_hmm_info.csv'%result_folder)==False:
    print('\n---------Calculating hmm prediction result using Seq2Seq-----------\n')
    hmm_res = func_fragment.seq2seq_nofaa(db_folder,result_folder,csv_faa_file)
    hmm_res.to_csv('%s/tmp/pred_hmm.csv'%result_folder,index=False)


#开始 singleton prediction

cmd = 'blastp -subject %s/singleton.faa -query %s/tmp/test_2.faa -out %s/tmp/blast_singleton.tsv -outfmt "6 qseqid sseqid pident qcovs length mismatch gapopen qstart qend sstart send evalue bitscore"' % (db_folder,result_folder,result_folder)

if os.path.exists('%s/tmp/blast_singleton.tsv'%result_folder)==False and os.path.getsize('%s/tmp/test_2.faa'%result_folder)>0:
    cmd_run(cmd)
    singleton_result = func_fragment.singleton_res(db_folder,result_folder)
    singleton_result.to_csv('%s/tmp/pred_singleton.csv'%result_folder,index=False)
elif os.path.getsize('%s/tmp/test_2.faa'%result_folder)==0:
    singleton_result = pd.DataFrame()
    singleton_result.to_csv('%s/tmp/pred_singleton.csv'%result_folder,index=False)
    singleton_result = singleton_result.drop(columns=['qcovs','evalue','pident'])


#计算最终结果
final_res = func_fragment.all_res(db_folder,result_folder)
print('\n---------%s/%s sequences have taxonomy------------\n'%(str(len(final_res)),str(len(df_test))))
final_res.to_csv('%s/final_result.csv'%result_folder,index=False)

final_all = time.time()

print('\nfinish_time:' + time_stamp())

print('\ntotal_time:' + str(math.floor((final_all-begin_all)/60))+ ' min ' + str(math.floor((final_all-begin_all)%60)) + ' s \n')


