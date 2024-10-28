import pandas as pd
from tqdm import tqdm
import json
import warnings
warnings.filterwarnings("ignore")

import sys
import os
import subprocess
import argparse
import shutil

import time
from Bio import SeqIO

import Virtaxa_fragment as func_fragment

parser = argparse.ArgumentParser(description='parameters of virtaxa', formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('--query', '-q', help='proteins of query contigs')
parser.add_argument('--db', '-d', help='reference database folder')
parser.add_argument('--output', '-o', help='output folder')
parser.add_argument('--new', '-n', action='store_true', default=False, help='new genera prediction option')
parser.add_argument('--fast', '-f', action='store_true', default=False, help='speed VirTAXA using the pre-built alignments')
args = parser.parse_args()

test_faa = args.query
db_folder = args.db
result_folder = args.output
new_option = args.new
fast_option = args.fast


test_csv = test_faa.split('.faa')[0] + '.csv'
test_fasta = test_faa.split('.faa')[0] + '.fasta'
train2test_file = '%s/blast_train2test.tsv'%result_folder
singleton2test_file = '%s/blast_singleton2test.tsv'%result_folder

train_file = '%s/train.csv'%db_folder
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

def cmd_run(cmd,print_or_not):
    if print_or_not:
        print(cmd)
    _ = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    _.communicate()

def make_fasta_csv(result_folder,train_file):
    df_ref_fasta = pd.read_csv(train_file)[['Accession','fasta']]
    test_fasta = SeqIO.parse('%s/test.fasta'%result_folder,'fasta')
    df_test = pd.read_csv('%s/test.csv'%result_folder)
    
    for record in test_fasta:
        df_test.loc[df_test['Accession']==record.id,'fasta']=str(record.seq)
    
    df_test_ref_fasta = pd.concat([df_ref_fasta,df_test])
    
    return df_test_ref_fasta

begin_all = time.time()

#主程序开始
    
if_exist_folder(result_folder)
if_exist_folder('%s/tmp'%result_folder)

cmd_run("cp %s %s/test.faa"%(test_faa,result_folder),True)
cmd_run("cp %s %s/test.fasta"%(test_fasta,result_folder),True)

df_test = pd.DataFrame()
record = SeqIO.parse(test_fasta,'fasta')
test_list = []
for rec in record:
    test_list.append(rec.id)
df_test['Accession']=test_list
df_test.to_csv('%s/test.csv'%result_folder,index=False)


df_test_ref_fasta = make_fasta_csv(result_folder,train_file)
df_test_ref_fasta.to_csv('%s/tmp/test_ref_fasta.csv'%result_folder,index=False)
    

#strict diamond prediction

cmd1 ='diamond makedb -p 52 --in %s/train.faa --db %s/train.dmnd'%(db_folder,db_folder)

cmd2 = 'diamond blastp --db %s/train.dmnd --query %s/test.faa --evalue 10 --out %s/tmp/diamond_blast_train2test.tsv --outfmt 6 qseqid sseqid pident qcovhsp length mismatch gapopen qstart qend sstart send evalue bitscore' % (db_folder,result_folder,result_folder)


if os.path.exists('%s/tmp/diamond_blast_train2test.tsv'%result_folder)==False:
    print('\n---------Calculating Strict BLASTP result-----------\n')
    cmd_run(cmd1,True)
    cmd_run(cmd2,True)


if os.path.exists('%s/tmp/test_1.faa'%result_folder)==False:
    pred_blast,test_1 = func_fragment.diamond_blast(db_folder,result_folder)
    
if os.path.getsize('%s/tmp/test_1.faa'%result_folder)==0:
    pred_blast['Step'] = 'Strict blast'
    pred_blast.to_csv('%s/final_result.csv'%result_folder,index=False)
else:
    if os.path.exists('%s/tmp/test_1_faa.csv'%result_folder)==False:
        test_faa_csv(result_folder).to_csv('%s/tmp/test_1_faa.csv'%result_folder,index=False)

    if os.path.exists('%s/tmp/hmmsearch_result_domain.abc'%result_folder)==False:
        begin = func_fragment.hmm_search(db_folder,result_folder)
    else:
        print('\n---------File %s/hmmsearch_result.out exists-----------\n'%result_folder)

    if os.path.exists('%s/tmp/pred_hmm_info.csv'%result_folder)==False:
        print('\n---------Calculating hmm prediction result-----------\n')
        hmm_res = func_fragment.hmm_result_calculation(db_folder,result_folder,fast_option)
        hmm_res.to_csv('%s/tmp/pred_hmm.csv'%result_folder,index=False)


    # singleton prediction

    cmd = 'blastp -subject %s/singleton.faa -query %s/tmp/test_2.faa -out %s/tmp/blast_singleton.tsv -outfmt "6 qseqid sseqid pident qcovs length mismatch gapopen qstart qend sstart send evalue bitscore"' % (db_folder,result_folder,result_folder)

    if os.path.exists('%s/tmp/blast_singleton.tsv'%result_folder)==False and os.path.getsize('%s/tmp/test_2.faa'%result_folder)>0:
        cmd_run(cmd,True)
        singleton_result = func_fragment.singleton_res(db_folder,result_folder)
        singleton_result.to_csv('%s/tmp/pred_singleton.csv'%result_folder,index=False)
    elif os.path.getsize('%s/tmp/test_2.faa'%result_folder)==0:
        singleton_result = pd.DataFrame()
        singleton_result.to_csv('%s/tmp/pred_singleton.csv'%result_folder,index=False)
        singleton_result = singleton_result.drop(columns=['qcovs','evalue','pident'])

    final_res = func_fragment.all_res(db_folder,result_folder)
    print('\n---------%s/%s sequences have taxonomy------------\n'%(str(len(final_res)),str(len(df_test))))
    final_res.to_csv('%s/final_result.csv'%result_folder,index=False)


df_other_level = func_fragment.other_level(result_folder, db_folder)
df_other_level.to_csv('%s/high_level.csv'%result_folder,index=False)

if new_option:
    print('predict new genera')
    fm_new_cluster = func_fragment.family_genus(db_folder,result_folder)

    tf = open("%s/new_cluster.json"%result_folder, "w")
    json.dump(fm_new_cluster,tf)
    tf.close()

shutil.rmtree('%s/tmp/'%result_folder)

final_all = time.time()

# print('\nfinish_time:' + time_stamp())

print('VirTAXA is done! Thank you for using')