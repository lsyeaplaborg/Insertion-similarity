#!/usr/bin/env python
# coding: utf-8

# In[91]:
import os
from os import path
from Bio.Seq import Seq 
from Bio import pairwise2
import re
import sys
# In[92]:
#%%
data_dir=sys.argv[1] 
out_dir=sys.argv[2] 
ref_fa=sys.argv[3]
file_suffix=sys.argv[4]
#%%

def align_func(tmp_seq1,tmp_seq2):
    #seq1 = Seq("CAAGGCGCGAGAAGTTCAAGAC") 
    #seq2 = Seq("AGTACAATGAGAAGTTCAAGAA")
    if len(tmp_seq1)==len(tmp_seq2):
        seq1 = Seq(tmp_seq1) 
        seq2 = Seq(tmp_seq2)
        alignments = pairwise2.align.globalxx(seq1, seq2)
        tmp_score=[]
        for alignment in alignments: 
            tmp_score.append(alignment.score)
        max_score=max(tmp_score)
    else:
        max_score=0
    return(max_score)

def match_two_seq(tmp_seq1,tmp_seq2):
    #seq1 = Seq("CAAGGCGCGAGAAGTTCAAGAC") 
    #seq2 = Seq("AGTACAATGAGAAGTTCAAGAA")
    tmp_score=0
    if len(tmp_seq1) == len(tmp_seq2):
        for i in range(len(tmp_seq1)): 
            if tmp_seq1[i]==tmp_seq2[i]:
                tmp_score+=1
    return(tmp_score)

#%%
def revise_reads(outdir,ref_seq,file_suffix):
    files=os.listdir(outdir)
    for i in files:
        if i.startswith('HQ') and i.endswith(file_suffix): 
            revised_line=[]
            f=open(path.join(outdir,i),'r')
            h1=f.readline().strip().split('\t')
            h1.append('revise_reads')
            revised_line.append(('\t').join(h1)+'\n')
            for line in f:
                if line.strip() != "":
                    line_ls=line.strip().split('\t')
                    score_same=line_ls[36]
                    score_split=line_ls[37]
                    score_local=line_ls[39]
                    # if score_split > score_same and score_split > score_local: # consider the local alignment socre
                    if score_split > score_same:
                        read_seq=line_ls[9]
                        cigar=line_ls[5]
                        ins=line_ls[25]
                        start_pos=line_ls[20]
                        align_pos=line_ls[1]
                        # reads_left=line_ls[29]
                        # reads_right=line_ls[30]
                        # ref_left=line_ls[34]
                        # ref_right=line_ls[35]
                        split_pos=line_ls[38]
                        
                        # revise cigar
                        cigar_num=[m for m in re.split('[MDI]',cigar) if m != '']
                        cigar_type=[m for m in re.split('[\d+]',cigar) if m != '']
                        n=0
                        for o in range(len(cigar_type)):
                            if cigar_type[o]=='M':
                                if n==0:
                                    cigar_num[o]=str(int(cigar_num[o])+int(split_pos))
                                elif n==1:
                                    cigar_num[o]=str(int(cigar_num[o])-int(split_pos))
                                n+=1
                        cigar_ls=[]
                        for j in range(len(cigar_type)):
                            cigar_ls.append(cigar_num[j])
                            cigar_ls.append(cigar_type[j])
                        cigar_v2=('').join(cigar_ls)
                        # revise start position 
                        start_pos_v2=str(int(start_pos)+int(split_pos))
                        # revise insertion
                        ins_v2=ins[(int(split_pos)):]+ins[:(int(split_pos))]
                        # revise left and right sequences
                        reads_left_v2=read_seq[(int(start_pos_v2)-int(align_pos)-len(ins)-1):(int(start_pos_v2)-int(align_pos)-1)]
                        reads_right_v2=read_seq[(int(start_pos_v2)-int(align_pos)+len(ins)-1):(int(start_pos_v2)-int(align_pos)+len(ins)*2-2)]
                        # revise left and right seq of reference 
                        ref_left_v2=ref_seq[(int(start_pos_v2)-len(ins)-1):(int(start_pos_v2)-1)]
                        ref_right_v2=ref_seq[(int(start_pos_v2)-1):(int(start_pos_v2)+len(ins)-1)]
                        # replace the items in the raw reads
                        line_ls[5]=cigar_v2
                        line_ls[25]=ins_v2
                        line_ls[20]=start_pos_v2
                        line_ls[29]=reads_left_v2
                        line_ls[30]=reads_right_v2
                        line_ls[34]=ref_left_v2
                        line_ls[35]=ref_right_v2
                        # add a marker col
                        line_ls.append('yes')
                    else:
                        line_ls.append('no')
                    revised_line.append(('\t').join(line_ls)+'\n')
            f.close()
            ## output 
            f_out=open(outdir+'/'+i.split('.')[0]+'_revise_reads.txt','w')
            for t in revised_line:
                f_out.write(t)
            f_out.close()


# In[95]:

def main(data_dir,out_dir,file_suffix):
    if not path.exists(out_dir):
        os.makedirs(out_dir)
    files=os.listdir(data_dir)
    for file in files:
        if file.startswith('HQ') and file.endswith(file_suffix): 
            print(file)
            f=open(path.join(data_dir,file),'r')
            #f=open('../data/HQ3997_A_12D_kedup.txt','r')
            h1=f.readline()
            ### splite insertion and calculate match score with left and right
            tmp_ls=[]
            for line in f:
                tmp=line.strip('\n').split('\t')
                if len(tmp)<=34:
                    line=line.strip('\n')+'\t'*(34-len(tmp))

                tmp=line.strip('\n').split('\t')
                # sequence of reads
                ins_seq=tmp[25]
                left_seq=tmp[29]
                right_seq=tmp[30]
                # sequence of reference
                start_pos=int(tmp[20])
                if start_pos-len(ins_seq)-1 >=0:
                    left_ref=ref_seq[start_pos-len(ins_seq)-1:start_pos-1]
                else:
                    left_ref=ref_seq[0:start_pos-1]
                if start_pos+len(ins_seq)-1 <=  len(ref_seq):
                    right_ref=ref_seq[start_pos-1:start_pos+len(ins_seq)-1]
                else:
                    right_ref=ref_seq[start_pos-1:len(ref_seq)]

                ## similiraty between insertions and reads/reference 
                same_score=[match_two_seq(ins_seq,left_seq),match_two_seq(ins_seq,right_seq),match_two_seq(ins_seq,left_ref),match_two_seq(ins_seq,right_ref)]
                max_same_score=max(same_score)
                max_same_ratio=max_same_score/len(ins_seq)
                tmp.append(left_ref)
                tmp.append(right_ref)
                tmp.append(str(max_same_ratio))
                ## split score
                if len(ins_seq) >= 6:
                    score_amount=[]
                    for n in range(3,(len(ins_seq)-2)): # split start at 3 bp of insertion
                        tmp_ins_l=ins_seq[:n]
                        tmp_ins_r=ins_seq[n:]
                        tmp_left=left_seq[n:]
                        tmp_right=right_seq[:n]
                        score_left=match_two_seq(tmp_ins_l,tmp_right)
                        score_right=match_two_seq(tmp_ins_r,tmp_left)
                        tmp_score_amount=score_left+score_right
                        score_amount.append(tmp_score_amount)
                    max_score=max(score_amount)
                    split_pos=score_amount.index(max_score)
                    tmp.append(str(max_score/len(ins_seq)))
                    tmp.append(str(split_pos+3))
                else:
                    tmp.append("0")
                    tmp.append("0")
                ## local alianment 
                same_score_local=[align_func(ins_seq,left_seq),align_func(ins_seq,right_seq),align_func(ins_seq,left_ref),align_func(ins_seq,right_ref)]
                max_local_score=max(same_score_local)
                max_score_ratio=max_local_score/len(ins_seq)
                tmp.append(str(max_score_ratio))
                tmp_str=('\t').join(tmp)
                tmp_ls.append(tmp_str)
        ### re-classify reads
            '''
            1: >3bp, same ratio>=60%
            2: 1bp, same ratio=100%
            3: 3bp, same ratio>=60%
            11: 2bp, same ratio=100%
            7: >=3bp, same ratio<60%
            12: 2bp, same ratio<100%
            13: 1bp, same ratio=0
            '''
            tmp_ls_revised=[]
            for k in tmp_ls:
                tmp_line=k.split('\t')
                tmp_line.append(str(max(float(tmp_line[36]),float(tmp_line[37]))))
                if tmp_line[28] in ['1','2','3','11','7','12','13']:
                    ins_len=len(tmp_line[25])
                    same_ratio=max(float(tmp_line[36]),float(tmp_line[37]))
                    if ins_len==1:
                        if same_ratio == 1:
                            tmp_line[28]='2'
                        else:
                            tmp_line[28]='13'
                    elif ins_len==2:
                        if same_ratio == 1:
                            tmp_line[28]='11'
                        else:
                            tmp_line[28]='12'                       
                    elif ins_len==3:
                        if same_ratio >= 0.6:
                            tmp_line[28]='3'
                        else:
                            tmp_line[28]='7'
                    elif ins_len>3:
                        if same_ratio >= 0.6:
                            tmp_line[28]='1'
                        else:
                            tmp_line[28]='7' 
                    tmp_ls_revised.append(('\t').join(tmp_line))
                else:
                    tmp_ls_revised.append(('\t').join(tmp_line))
            ### output
            f=open(path.join(out_dir,file.split('.')[0]+'_new.txt'),'w')
            f.write(h1.strip('\n')+'\t'+'ref_left'+'\t'+'ref_right'+'\t'+'same_ratio_v2'+'\t'+'split_ratio'+'\t'+'split_pos'+'\t'+'local_align_ratio'+'\t'+'merged_ratio'+'\n')
            for i in tmp_ls_revised:
                f.write(i+'\n')
            f.close()


# In[96]:
if __name__ == '__main__':
    # outdir='/Users/lengsiewyeap/Desktop/haoqian/duplicate_insertion/pipeline/results/mouse'
    # data_dir='/Users/lengsiewyeap/Desktop/haoqian/duplicate_insertion/pipeline/results/mouse'

    f=open(ref_fa,'r')
    h2=f.readline()
    ref_seq=f.readline().strip().upper()
    f.close()

    main(data_dir,out_dir,file_suffix)
    # revise cigar
    revise_reads(out_dir,ref_seq,file_suffix.split('.')[0]+'_new.txt')
# %%
