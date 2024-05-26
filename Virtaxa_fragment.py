import pandas as pd
import numpy as np

from tqdm import tqdm
import os
import subprocess
import copy
import time

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import Phylo

#phmm cutoff
H_score = 0.1
H_cov = 0.95

#step 1 cutoff
B_evalue = 1e-10
B_cov = 95
B_ident = 40
B_multi = 50

#step 4 cutoff
S_evalue = 1e-5
S_cov = 75

def cmd_run(cmd):
    print(cmd)
    _ = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    _.communicate()

def diamond_blast(db_folder,result_folder):
    df_diamond_blast = pd.read_csv('%s/tmp/diamond_blast_train2test.tsv'%result_folder,sep='\t',names=['qseqid','sseqid','pident', 'qcovs', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])
    df_diamond_blast['qseqid'] = df_diamond_blast['qseqid'].str[0:-2].str.strip('_')
    df_diamond_blast['sseqid'] = df_diamond_blast['sseqid'].str[0:-2].str.strip('_')
    df_diamond_blast = df_diamond_blast[['qseqid','sseqid','pident', 'qcovs', 'evalue']]
    df_train = pd.read_csv('%s/train.csv'%db_folder)
    df_test = pd.read_csv('%s/test.csv'%result_folder)
    
    df_diamond_blast = df_diamond_blast[df_diamond_blast['evalue']<B_evalue]
    
    df_diamond_blast['multi'] = round(df_diamond_blast['qcovs']*df_diamond_blast['pident']/100,1)
    
    df_diamond_blast = df_diamond_blast[df_diamond_blast['multi']>B_multi]
    
    for seq in tqdm(df_diamond_blast['qseqid'].unique().tolist()):
        df1 = df_diamond_blast[df_diamond_blast['qseqid'] == seq]
        df1 = df1[df1['evalue'] == df1['evalue'].min()]
        df1 = df1[df1['pident'] == df1['pident'].max()][0:1]
        pred_genus = df_train[df_train['Accession']==df1['sseqid'].values[0]]['Genus'].values[0]
        df_test.loc[df_test['Accession']==seq,'pred_genus']=pred_genus
        df_test.loc[df_test['Accession']==seq,'qcovs']=df1['qcovs'].values[0]
        df_test.loc[df_test['Accession']==seq,'evalue']=df1['evalue'].values[0]
        df_test.loc[df_test['Accession']==seq,'pident']=df1['pident'].values[0]
    df_test = df_test.dropna(subset='pred_genus')
    df_test_list = df_test['Accession'].tolist()

    next_test = []
    for rec in SeqIO.parse('%s/test.faa'%result_folder,'fasta'):
        recid = rec.id[0:-2].strip('_')
        if recid not in df_test_list:
            next_test.append(rec)
    df_test.to_csv('%s/tmp/pred_diamond_blast.csv'%result_folder,index=False)
    SeqIO.write(next_test,'%s/tmp/test_1.faa'%result_folder,'fasta')
    return df_test,next_test

def hmm_search(db_folder,result_folder):
    print('\n---------begin hmmsearch-----------\n')
    cmd =  'hmmsearch --tblout %s/tmp/hmmsearch_result.out --domtblout %s/tmp/hmmsearch_result_domain.out %s/hmm_train.hmm %s/tmp/test_1.faa '%(result_folder,result_folder,db_folder,result_folder)
    print(cmd)
    res = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    res.communicate()
    
    print('\n----begin reformat hmmsearch-----\n')
    
    hmmsearch_out_filename = '%s/tmp/hmmsearch_result_domain.out'%result_folder
    lines = open(hmmsearch_out_filename).readlines()
    with open(hmmsearch_out_filename, mode='a') as filename:
        filename.truncate(0)
        for line in lines:
            if line[0]!='#':
                filename.write(line)
    
    cmd = 'awk \'{print $1,$3,$4,$7,$8,$18,$19}\' %s/tmp/hmmsearch_result_domain.out > %s/tmp/hmmsearch_result_domain.abc'%(result_folder,result_folder)
    res = subprocess.Popen(cmd, shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE) # 使用管道
    res.communicate()

def seq_fa(seq,result_folder):
    rna_yp_1 = pd.read_csv('%s/tmp/test_1_faa.csv'%result_folder)
    rna_yp = copy.deepcopy(rna_yp_1)
    rna_yp['Accession'] = rna_yp['Accession'].str[0:-2].str.strip('_')
    dfp = rna_yp[rna_yp['Accession']==seq]
    dfp = dfp.rename(columns={'Accession':'Nucl'})
    new_seq = []
    if len(dfp)>0:
        df0 = pd.merge(rna_yp_1,dfp,on=['fasta'])
        for index,row in df0.iterrows():
            yp = row['Accession']
            fa = row['fasta']
            des = 'rna virus'
            seq_rec = SeqRecord(Seq(fa), description=des, id=yp, name=yp)
            new_seq.append(seq_rec)
    else:
        print("can not find sequences in csv_fasta:"+seq)
    SeqIO.write(new_seq, '%s/tmp/seq.faa'%(result_folder), "fasta")

def pc_(mcl_file):
    f = open(mcl_file)
    lines = f.readlines()
    pc = {}
    i = 1
    for line in lines:
        h = line.strip().split('\t')
        h1 = []
        for hh in h:
            j = hh
            h1.append(j)
        if(len(h1)>=2):
            pc.update({'cluster%s'%str(i):h1})
            i=i+1
    return pc

def seq2seq_nofaa(db_folder,result_folder,csv_faa_file):
    linux_file=[]
    cmd_run("echo > %s/tmp/job_hmmer_decision.sh"%(result_folder))
    pc = pc_('%s/mcl_train.out'%db_folder)
    df_cl_info = pd.read_csv('%s/mcl_info.csv'%db_folder)
    df_cl_info['cutoff'] = df_cl_info['negative'] + (df_cl_info['positive'] - df_cl_info['negative'])*H_score
    
    blast_train2test = pd.read_csv('%s/tmp/diamond_blast_train2test.tsv'%result_folder,sep='\t',names=['qseqid','sseqid','pident', 'qcovs', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])
    blast_train2test['qseqid'] = blast_train2test['qseqid'].str[0:-2].str.strip('_')
    blast_train2test['sseqid'] = blast_train2test['sseqid'].str[0:-2].str.strip('_')
    blast_train2test = blast_train2test[['qseqid','sseqid','pident', 'qcovs', 'evalue']]
    
    df_train = pd.read_csv('%s/train.csv'%db_folder)
    df_test = pd.read_csv('%s/test.csv'%result_folder)
    df_test['reject_or_not'] = None
    df_test['tree_validation'] = None

    df_hms = pd.read_csv('%s/tmp/hmmsearch_result_domain.abc'%result_folder,sep='\s+',
                     names=['Accession','Length','Cluster','evalue','score','from','to'])
    df_hms['coverage']=round((df_hms['to']-df_hms['from'])/df_hms['Length'],2)
    seq_list = df_hms['Accession'].unique()
    df_final = pd.DataFrame()
    for seq in tqdm(seq_list):
        dfp = df_hms[df_hms['Accession']==seq]   
        cl_list = dfp['Cluster'].unique()
        for cl in cl_list:
            dfpp = dfp[dfp['Cluster']==cl]
            df_new = dfpp[0:1]
            df_new['coverage']=dfpp['coverage'].sum()
            df_final = pd.concat([df_final,df_new])
    df_final.Accession = df_final.Accession.str[0:-2].str.strip('_')
    df_final = df_final.sort_values(by=["Accession","score"], ascending=False)
    df_final.to_csv('%s/tmp/hmmsearch_result_all.csv'%result_folder,index=False)
    
    df_result = pd.DataFrame()
    
    for seq in tqdm(df_final['Accession'].unique()):
        dfd = df_final[df_final['Accession']==seq]
        dfd = dfd.sort_values(by=["coverage"], ascending=False)[0:1]
        hmm_coverage = dfd['coverage'].values[0]
        cluster1 = dfd['Cluster'].values[0]
        df_b = blast_train2test[blast_train2test['qseqid']==seq]
        seq_res = pd.DataFrame()
        score =  dfd['score'].values[0]
        cutoff = df_cl_info[df_cl_info['cluster']==cluster1]['cutoff'].values[0]
        cluster_negative = df_cl_info[df_cl_info['cluster']==cluster1]['negative'].values[0]

        pred_genus = df_cl_info[df_cl_info['cluster']==cluster1]['gn_rate'].values[0].split('\'')[1]
        
        if pred_genus:
            df_test.loc[df_test['Accession']==seq,'pred_genus']=pred_genus

            # df_test.loc[df_test['Accession']==seq,'hit_score']=score
            
            ## decide reject or not
            if dfd['score'].values[0]>cutoff:
                df_test.loc[df_test['Accession']==seq,'reason']='High score'
            elif cluster_negative == 0:
                df_test.loc[df_test['Accession']==seq,'reason']='Specific cluster'
            elif hmm_coverage > H_cov:
                df_test.loc[df_test['Accession']==seq,'reason']='Coverage'
            else:
                # print(seq + ' unsure')
                df_test.loc[df_test['Accession']==seq,'tree_validation']='1'
                df_pred = df_test[df_test['Accession']==seq]
                reject_decision = hmm_reject_decision(dfd,df_pred,pc,db_folder,result_folder,blast_train2test,csv_faa_file)
                if reject_decision:
                    linux_file.append(reject_decision)
                else:
                    df_test.loc[df_test['Accession']==seq,'reject_or_not']='lack enough sequences'
                with open("%s/tmp/decision_file.txt"%result_folder, "w") as tf:
                    tf.write(str(linux_file))

   
    if len(linux_file)>0:
        print('number of sequences need check: ' + str(len(linux_file)))
        last_tree_file_name = linux_file[-1][1]
        print('last tree file name: %s'%last_tree_file_name)

        if os.path.exists(last_tree_file_name)==False:
            cmd0 = "sh %s/tmp/job_hmmer_decision.sh"%result_folder
            cmd_run(cmd0)

        while True:
            if os.path.exists(last_tree_file_name):
                file_size = os.path.getsize(last_tree_file_name)
                if file_size > 0:
                    break
            else:
                time.sleep(5)
                continue  

        for info in tqdm(linux_file):
            print('tree validation:%s'%info[2])
            des_result = treeshrink(info[1],info[2])
            df_test.loc[df_test['Accession']==info[0],'reject_or_not']=des_result
    
    df_test = df_test.dropna(subset='pred_genus')
    
    df_test.to_csv('%s/tmp/pred_hmm_info.csv'%result_folder,index=False)
    
    df_test_1 = df_test[df_test['tree_validation']!=df_test['tree_validation']]
    df_test_2 = df_test[df_test['reject_or_not']=='keep']
    
    df_test = pd.concat([df_test_1,df_test_2])
    
    df_test_list = df_test['Accession'].unique().tolist()
    
    next_test = []
    for rec in SeqIO.parse('%s/tmp/test_1.faa'%result_folder,'fasta'):
        recid = rec.id[0:-2].strip('_')
        if recid not in df_test_list:
            next_test.append(rec)
    SeqIO.write(next_test,'%s/tmp/test_2.faa'%result_folder,'fasta')
    
    return df_test

def hmm_reject_decision(df_hmm, df_pred, pc,db_folder,result_folder,blast_train2test,csv_faa_file):
    if os.path.exists('%s/phy_tree'%result_folder)==False:
        os.mkdir('%s/phy_tree'%result_folder)
        
    all_faa = pd.read_csv('%s/tmp/test_ref_fasta.csv'%result_folder)
    contig = df_pred['Accession'].values[0]    
    df_cl_info = pd.read_csv('%s/mcl_info.csv'%db_folder)        
    df_train = pd.read_csv('%s/train.csv'%db_folder)
    df_train = df_train[['Accession','Family','Genus','Segment']]

    new_file = []
    best_align_cluster = df_hmm['Cluster'].values[0]
    pred_genus = df_pred['pred_genus'].values[0]
    pred_segment = df_cl_info[df_cl_info['cluster']==best_align_cluster]['Segment'].values[0]
    print(contig)
    
    if pred_segment == 'complete':##if align to complete genus，build tree with sequences from this genus
        print('%s,complete genus'%pred_genus)
        df_g = df_train[df_train['Genus']==pred_genus]
        for seq in df_g['Accession'].unique():
            new_id = seq + '_%s'%pred_genus
            des = pred_genus
            fa = all_faa[all_faa['Accession']==seq]['fasta'].values[0]
            new_rec = SeqRecord(Seq(fa), description=des, id=new_id, name=seq)
            new_file.append(new_rec)
        file_name = '%s/phy_tree/fasta_%s_%s.fasta'%(result_folder,contig,pred_genus)
    else:
        print('%s,segment genus'%pred_genus)
        if pred_segment.find(',')!=-1: ##if aligned to the segmented genus, buld tree with the sequences from this specific segment of this genus
            print('isconfused segment')
            list_seg = []
            pred_segment = pred_segment.strip('{')
            pred_segment = pred_segment.strip('}')
            pred_segment = pred_segment.split(',')
            df_g = pd.DataFrame()
            for i in pred_segment:
                i = i.strip('\'')
                list_seg.append(i)
            print(list_seg)
            for seg_i in list_seg:
                df_ = df_train[np.logical_and(df_train['Genus']==pred_genus,df_train['Segment']==seg_i)]
                df_g = pd.concat([df_g,df_])
                for seq in df_['Accession'].unique():
                    if contig.find(seq)==-1:
                        new_id = seq + '_%s_segment_%s'%(pred_genus,seg_i)
                        des = pred_genus
                        fa = all_faa[all_faa['Accession']==seq]['fasta'].values[0]
                        new_rec = SeqRecord(Seq(fa), description=des, id=new_id, name=seq)
                        new_file.append(new_rec)
            file_name = '%s/phy_tree/fasta_%s_%s_segment_confused.fasta'%(result_folder,contig,pred_genus)
        else:
            # print('unconfused segment')
            print(pred_genus,pred_segment)
            df_g = df_train[np.logical_and(df_train['Genus']==pred_genus,df_train['Segment']==pred_segment)]
            for seq in df_g['Accession'].unique():
                if contig.find(seq)==-1:
                    new_id = seq + '_%s_segment%s'%(pred_genus,pred_segment)
                    des = pred_genus
                    fa = all_faa[all_faa['Accession']==seq]['fasta'].values[0]
                    new_rec = SeqRecord(Seq(fa), description=des, id=new_id, name=seq)
                    new_file.append(new_rec)
            file_name = '%s/phy_tree/fasta_%s_%s_segment%s.fasta'%(result_folder,contig,pred_genus,pred_segment)
    if len(new_file)<=1:
        file_name = ''
        des_result =  'database lack enough segment sequences'
        print('Part Reject:do not have enough segment sequences—%s'%(contig))
    else:
        new_id = contig + '_novel'
        des = 'des'
        fa = all_faa[all_faa['Accession']==contig]['fasta'].values[0]
        new_rec = SeqRecord(Seq(fa), description=des, id=new_id, name=contig)
        new_file.append(new_rec)
        
        SeqIO.write(new_file,'%s'%file_name,'fasta')

        aln_file_name = file_name.split('.')[0] + '_aln.fasta'
        tree_file_name = file_name.split('.')[0] + '_tree.treefile'
        cmd1 = 'mafft --thread 56 --auto %s > %s'%(file_name,aln_file_name)
        cmd2 = 'FastTree -nt %s > %s'%(aln_file_name,tree_file_name)
        if os.path.exists(tree_file_name)==False:
            with open("%s/tmp/job_hmmer_decision.sh"%result_folder, mode='a') as filename:
                filename.write('\n' + cmd1)
                filename.write('\n' + cmd2)
        return contig,tree_file_name,new_id

def treeshrink(treefile,contig):
    print(contig)
    out_folder = treefile.split('.treefile')[0]
    # cmd1 = 'python ../TreeShrink-master/run_treeshrink.py  -t %s -o %s'%(treefile,out_folder)
    cmd1 = 'run_treeshrink.py  -t %s -o %s'%(treefile,out_folder)
    cmd_run(cmd1)
    if os.path.exists('%s/output_summary.txt'%(out_folder)):
        df_treeshrink=pd.read_csv('%s/output_summary.txt'%(out_folder),sep='\s+')
        max_score = df_treeshrink['Signature'].max()
        df_novel = df_treeshrink[df_treeshrink['Species']==contig]
        if len(df_novel)>0 and df_novel['Signature'].values[0]==max_score:
            result = 'reject'
        else:
            print('tree_shrink keep')
            result='keep'
    else:
        cmd_run('rm -r %s'%out_folder)
        print('self check')
        result = dist(treefile)
        print(result)
    
    return result

def dist(treefile):
    tree = Phylo.read(treefile,'newick')
    dis_list = []
    dis_dict = {}

    for leaf in tree.get_terminals():
        dis = tree.distance({"name": "%s"%leaf})
        if leaf.name.find('novel')==-1:
            dis_list.append(dis)
        else:
            novel_name = leaf.name
            novel_seq = leaf
        dis_dict.update({leaf:dis})
    mean_dis = round(sum(dis_list) / len(dis_list),2)
    max_dis = round(max(dis_list),2)
    novel_dis = round(dis_dict[novel_seq],2)
    if novel_dis > max_dis:
        partner = if_partner(tree,novel_name)
        # print(novel_dis,mean_dis)
        if partner == 1 and novel_dis <= 2*mean_dis:
            result = 'keep'
        else:
            result = 'reject'
    else:
        result = 'keep'
    
    return result

def if_partner(tree,novel_name):
    result = 0
    if len(tree.get_path(novel_name))>1:
        last_cl = tree.get_path(novel_name)[-2]
        if len(last_cl.get_nonterminals())==1:
            print('have Sibling node：%s'%novel_name)
            result = 1
    return result


def singleton_res(db_folder,result_folder):    
    blast_sing = pd.read_csv('%s/tmp/blast_singleton.tsv'%result_folder,sep='\t',names=['qseqid','sseqid','pident', 'qcovs', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])
    blast_sing['qseqid'] = blast_sing['qseqid'].str[0:-2].str.strip('_')
    blast_sing['sseqid'] = blast_sing['sseqid'].str[0:-2].str.strip('_')
    blast_sing = blast_sing[['qseqid','sseqid','pident', 'qcovs', 'evalue']]
    
    blast_sing = blast_sing[blast_sing['qcovs']>S_cov]
    blast_sing = blast_sing[blast_sing['evalue']<S_evalue]
    
    df_test = pd.read_csv('%s/test.csv'%result_folder)
    df_train = pd.read_csv('%s/train.csv'%db_folder)
    df_test['pred_genus'] = None
    
    singleton_faa = SeqIO.parse('%s/singleton.faa'%db_folder,'fasta')
    blast_singleton = pd.DataFrame()
    singleton_list = []
    for seq in singleton_faa:
        seqid = seq.id[0:-2].strip('_')
        df_new = blast_sing[blast_sing['sseqid']==seqid]
        blast_singleton = pd.concat([blast_singleton,df_new])
    for seq in df_test['Accession'].tolist():
        dfp = blast_singleton[blast_singleton['qseqid']==seq]
        if dfp.empty==False:
            max_sseqid = dfp[dfp.evalue == dfp.evalue.min()][0:1]['sseqid'].values[0]          
            pred_genus = df_train[df_train['Accession']==max_sseqid]['Genus'].values[0]
            df_test.loc[df_test['Accession']==seq,'pred_genus']=pred_genus
    df_test = df_test.dropna(subset='pred_genus')
    return df_test
        
    
def all_res(db_folder,result_folder):
    strict_blast_result = pd.read_csv('%s/tmp/pred_diamond_blast.csv'%result_folder)
    strict_blast_result['Step'] = 'Strict blast'
    strict_blast_result = strict_blast_result[['Accession','Step','pred_genus']]
    
    
    hmm_result = pd.read_csv('%s/tmp/pred_hmm.csv'%result_folder)  
    hmm_result = hmm_result[hmm_result['pred_genus']==hmm_result['pred_genus']]
    
    hmm_result_1 = hmm_result[hmm_result['tree_validation']!=1]
    hmm_result_2 = hmm_result[hmm_result['tree_validation']==1]
    
    hmm_result_1['Step'] = 'pHMM'
    hmm_result_2['Step'] = 'Phylo Tree'
    hmm_result_1 = hmm_result_1[['Accession','Step','pred_genus']]
    hmm_result_2 = hmm_result_2[['Accession','Step','pred_genus']]
    
    
    singleton_result = pd.read_csv('%s/tmp/pred_singleton.csv'%result_folder)
    singleton_result['Step'] = 'Singleton'
    singleton_result = singleton_result.dropna(subset='pred_genus')
    
    singleton_result = singleton_result[['Accession','Step','pred_genus']]
    
    
    df_test = pd.read_csv('%s/test.csv'%result_folder)
    
    print('strict blast result:' + str(len(strict_blast_result)))
    print('pHMMs result:' + str(len(hmm_result_1)))
    print('Phylogenetic tree result:' + str(len(hmm_result_2)))
    print('singleton result:' + str(len(singleton_result)))
    
    df_test_genus = pd.concat([strict_blast_result,hmm_result_1,hmm_result_2,singleton_result])

    print('final result:' + str(len(df_test_genus)))
    
    return df_test_genus
   
def build_fm_tree(res_folder,db_folder):
    res_list = pd.read_csv('%s/final_result.csv'%res_folder)['Accession'].tolist()
    test_list = pd.read_csv('%s/test.csv'%res_folder)['Accession'].tolist()
    no_res = list(set(test_list) - set(res_list))
    df_hmm = pd.read_csv('%s/tmp/hmmsearch_result_all.csv'%res_folder)
    df_f = pd.DataFrame()
    mcl_info = pd.read_csv('%s/mcl_info.csv'%db_folder)
    for seq in no_res:
        df_ = df_hmm[df_hmm['Accession']==seq]
        df_f = pd.concat([df_f,df_])
    for cl in df_f['Cluster'].unique():
        sc = mcl_info[mcl_info['cluster']==cl]['cutoff'].values[0]
        fm = mcl_info[mcl_info['cluster']==cl]['major_family'].values[0]
        df_f.loc[df_f['Cluster']==cl,'cutoff']=sc
        df_f.loc[df_f['Cluster']==cl,'pred_family']=fm
    df_f = df_f[df_f['score']>0.5*df_f['cutoff']]
    df_f = df_f[['Accession','pred_family']]
    df_f.to_csv('%s/final_result_family.csv'%res_folder,index=False)
    fm_new_cluster = family_genus(df_f,db_folder,res_folder)
    
    return fm_new_cluster

def family_genus(df_fm,db_folder,res_folder):
    all_faa = pd.read_csv('%s/tmp/test_ref_fasta.csv'%res_folder)
    df_train = pd.read_csv('%s/train.csv'%db_folder)
    new_file = []
    cmd_run("echo > %s/tmp/job_hmmer_decision_2.sh"%(res_folder))
    if os.path.exists('%s/phy_tree_family/'%(res_folder))==False:
        os.mkdir('%s/phy_tree_family/'%(res_folder))
    list1 = []
    family_list = []
    for fm in tqdm(df_fm['pred_family'].unique()):
        new_file = []
        df_ref_family = df_train[df_train['Family']==fm]
        dfp = df_fm[df_fm['pred_family']==fm]
        if len(dfp)>=3:
            print(fm,len(dfp))
            family_list.append(fm)
            for seq in dfp['Accession'].unique():
                new_id = seq + '_novel'
                des = 'novel'
                fa = all_faa[all_faa['Accession']==seq]['fasta'].values[0]
                new_rec = SeqRecord(Seq(fa), description=des, id=new_id, name=seq)
                new_file.append(new_rec)
            for seq in df_ref_family['Accession'].unique():
                gn = df_train[df_train['Accession']==seq]['Genus'].values[0]
                new_id = seq + '_%s_%s'%(fm,gn)
                des = fm
                fa = all_faa[all_faa['Accession']==seq]['fasta'].values[0]
                new_rec = SeqRecord(Seq(fa), description=des, id=new_id, name=seq)
                new_file.append(new_rec)
            file_name = '%s/phy_tree_family/family_%s.fasta'%(res_folder,fm)
            aln_file_name = '%s/phy_tree_family/family_%s_aln.fasta'%(res_folder,fm)
            tree_file_name = '%s/phy_tree_family/family_%s_aln.treefile'%(res_folder,fm)
            SeqIO.write(new_file,file_name,'fasta')
            cmd1 = 'mafft --thread 56 --auto %s > %s'%(file_name,aln_file_name)
            cmd2 = 'FastTree -nt %s > %s'%(aln_file_name,tree_file_name)
            if os.path.exists(tree_file_name)==False:
                list1.append(tree_file_name)
                with open("%s/tmp/job_hmmer_decision_2.sh"%res_folder, mode='a') as filename:
                    filename.write('\n' + cmd1)
                    filename.write('\n' + cmd2)
    fm_new_cluster = {} 
    if len(family_list)>0:
        if len(list1)>0:
            cmd_run('sh %s/tmp/job_hmmer_decision_2.sh'%res_folder)
            while True:
                print(os.path.exists(list1[-1]))
                if os.path.exists(list1[-1]):
                    file_size = os.path.getsize(list1[-1])
                    if file_size > 0:
                        break
                else:
                    time.sleep(5)
                    continue
        for fm in family_list:
            tree_file = '%s/phy_tree_family/family_%s_aln.treefile'%(res_folder, fm)
            if os.path.exists(tree_file):
                cluster_list = find_new_cluster(tree_file)
                if len(cluster_list)>0:
                    fm_new_cluster.update({fm:cluster_list})
    return fm_new_cluster
   

def find_new_cluster(tree_file):
    tree = Phylo.read(tree_file, "newick")
    clade = tree.get_nonterminals()
    novel_clade = []
    for cl in clade:
        novel_num = 0
        tm_num = len(cl.get_terminals())
        for tm in cl.get_terminals():
            if tm.name.find('novel')!=-1:
                novel_num = novel_num + 1
        if novel_num >= 0.8 * tm_num and tm_num>2:
            print(novel_num, tm_num)
            novel_clade.append(cl)
    # print(novel_clade)

    max_clade = []
    cluster_list = []
    for cl in novel_clade:
        node_list = []
        cl_max = 1
        for cl1 in novel_clade:
            if cl1!=cl and cl in cl1.get_nonterminals():
                cl_max = 0
                break
        if cl_max == 1:
            if len(cl.get_nonterminals())>1:
                max_clade.append(cl)
                cl.color = 'red'
                for node in cl.get_terminals():
                    node_list.append(node.name)
        cluster_list.append(node_list)
    # Phylo.draw(tree)
    
    return cluster_list