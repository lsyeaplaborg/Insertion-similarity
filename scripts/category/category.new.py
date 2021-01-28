#Python3
import re
import sys

'''
1. cate.xls
2. fasta
3. ins
4. out_kedup
5. out_nokedup

'''

#in_fa="VB18_F3_short.fa"
in_fa = sys.argv[2]

ref = ''.join([t.strip().upper() for t in open(in_fa).readlines() if not t.startswith('>')])

##the number of different base
def diff_base(seq1,seq2):
	num = 0
	if len(seq1)==len(seq2):
		for i in range(len(seq1)):
			if seq1[i] != seq2[i]:
				num += 1
	return(num)

## the ratio of same base
def same_ratio(seq1,seq2):
	num = 0
	if len(seq1)==len(seq2):
		for i in range(len(seq1)):
			if seq1[i] == seq2[i]:
				num += 1
	ratio=num/len(seq1)
	return(ratio)

line_num = 0

out_kedup=open(sys.argv[4],'w')
out_nokedup=open(sys.argv[5],'w')

dict1={}

#for line in open("HQ4664_A_12D.cate.xls"):
for line in open(sys.argv[1]):
	cigar=line.strip().split('\t')[5]
	seq=line.strip().split('\t')[9]
	dst=line.strip().split('\t')[20]
	insize=line.strip().split('\t')[24]
	inseq=line.strip().split('\t')[25]
	category=line.strip().split('\t')[28]

	p = r'(\d*)(I|D|M)'
	m = re.findall(p,cigar)

	before_seq=""
	after_seq=""
	same_ratio_result=""
	basefrom=""
	baseto=""
#	print(start3)
	line_num += 1
	if line_num > 1:
		if int(category) == 5:
			#*M*I*M*I*M
			if int(m[2][0]) <= 5:
				start1=int(m[0][0])
				start2=int(m[1][0])
				start3=int(m[2][0])
				start4=int(m[3][0])
				start5=int(m[4][0])
				other=""
				if len(m) > 5:
					for i in range(5,len(m)):
						other += "".join(m[i])
				else:
					other=""
				m_seq=seq[int(start1)+int(start2):int(start1)+int(start2)+int(start3)]
				m_ins_left=seq[start1:start1+start3]
				m_ins_right=seq[start1+start2+start3+start4-start3:start1+start2+start3+start4]
	
				#M in middle convert to left
				ins_left_new=seq[start1+start3:start1+start3+start2+start4]
				ins_left_new_left=seq[start1+start3-start2-start4:start1+start3]
				ins_left_new_right=seq[start1+start3+start2+start4:start1+start3+start2+start4+start2+start4]
				same_ratio_seq_left_left=same_ratio(ins_left_new,ins_left_new_left)
				same_ratio_seq_left_right=same_ratio(ins_left_new,ins_left_new_left)
				
				ins_left_ref_left=ref[start1+start3-start2-start4:start1+start3]
				ins_left_ref_right=ref[start1+start3:start1+start3+start2+start4]
				same_ratio_ref_left_left=same_ratio(ins_left_new,ins_left_ref_left)
				same_ratio_ref_left_right=same_ratio(ins_left_new,ins_left_ref_right)
	
				#M in middle convert to right
				ins_right_new=seq[start1:start1+start2+start4]
				ins_right_new_left=seq[start1-start2-start4:start1]
				ins_right_new_right=seq[start1+start2+start4:start1+start2+start4+start2+start4]
				same_ratio_seq_right_left=same_ratio(ins_right_new,ins_right_new_left)
				same_ratio_seq_right_right=same_ratio(ins_right_new,ins_right_new_right)
				
				ins_right_ref_left=ref[start1-start2-start4:start1]
				ins_right_ref_right=ref[start1:start1+start2+start4]
				same_ratio_ref_right_left=same_ratio(ins_right_new,ins_right_ref_left)
				same_ratio_ref_right_right=same_ratio(ins_right_new,ins_right_ref_right)
	
				if diff_base(m_seq,m_ins_left) > diff_base(m_seq,m_ins_right) and diff_base(m_seq,m_ins_right) <= 1:
					insize=start2+start4
					inseq=ins_right_new
					dst=start1+1
					cigar=str(start1)+'M'+str(start2+start4)+'I'+str(start5+start3)+'M'+other
					same_ratio_result=max(same_ratio_seq_right_left,same_ratio_seq_right_right,same_ratio_ref_right_left,same_ratio_ref_right_right)
					before_seq=ins_right_new_left
					after_seq=ins_right_new_right
					if max(same_ratio_seq_right_left,same_ratio_ref_right_left) >= max(same_ratio_seq_right_right,same_ratio_ref_right_right):
						basefrom = start1-insize+1
						baseto = start1
					else:
						basefrom = start1+1
						baseto = start1+insize
				elif diff_base(m_seq,m_ins_left) < diff_base(m_seq,m_ins_right) and diff_base(m_seq,m_ins_left) <= 1:
					insize=start2+start4
					inseq=ins_left_new
					dst=start1+start3+1
					cigar=str(start1+start3)+'M'+str(start2+start4)+'I'+str(start5)+'M'+other
					same_ratio_result=max(same_ratio_seq_left_left,same_ratio_seq_left_right,same_ratio_ref_left_left,same_ratio_ref_left_right)
					before_seq=ins_left_new_left
					after_seq=ins_left_new_right
					if max(same_ratio_seq_left_left,same_ratio_ref_left_left) >= max(same_ratio_seq_left_right,same_ratio_ref_left_right):
						basefrom = start1+start3-insize+1
						baseto = start1+start3
					else:
						basefrom = start1+start3+1
						baseto = start1+start3+insize

				elif diff_base(m_seq,m_ins_left) == diff_base(m_seq,m_ins_right) and diff_base(m_seq,m_ins_left) <= 1:
					insize=start2+start4
					same_ratio_left=max(same_ratio_seq_left_left,same_ratio_seq_left_right,same_ratio_ref_left_left,same_ratio_ref_left_right)
					same_ratio_right=max(same_ratio_seq_right_left,same_ratio_seq_right_right,same_ratio_ref_right_left,same_ratio_ref_right_right)
					if same_ratio_left >= same_ratio_right:
						same_ratio_result = same_ratio_left
						inseq=ins_left_new
						dst=start1+start3+1
						cigar=str(start1+start3)+'M'+str(start2+start4)+'I'+str(start5)+'M'+other
						before_seq=ins_left_new_left
						after_seq=ins_left_new_right
						if max(same_ratio_seq_left_left,same_ratio_ref_left_left) >= max(same_ratio_seq_left_right,same_ratio_ref_left_right):
							basefrom = start1+start3-insize+1
							baseto = start1+start3
						else:
							basefrom = start1+start3+1
							baseto = start1+start3+insize
					else:
						same_ratio_result = same_ratio_right
						inseq=ins_right_new
						dst=start1+1
						cigar=str(start1)+'M'+str(start2+start4)+'I'+str(start5+start3)+'M'+other
						before_seq=ins_right_new_left
						after_seq=ins_right_new_right
						if max(same_ratio_seq_right_left,same_ratio_ref_right_left) >= max(same_ratio_seq_right_right,same_ratio_ref_right_right):
							basefrom = start1-insize+1
							baseto = start1
						else:
							basefrom = start1+1
							baseto = start1+insize

#					insize=start2+start4
				else:
					cigar=cigar
#				print(str(diff_base(m_seq,m_ins_left))+'\t'+str(diff_base(m_seq,m_ins_right)))
				if cigar.count("I")==1:				
					if start2+start4 >3 and same_ratio_result >= 0.6:
						category = 1
					elif start2+start4 >3 and same_ratio_result < 0.6:
						category = 7
					elif start2+start4 == 3 and same_ratio_result >= 0.6:
						category = 3
					elif start2+start4 == 3 and same_ratio_result < 0.6:
						category = 7
					elif start2+start4 == 1 and same_ratio_result >= 0.6:
						category = 2
					elif start2+start4 == 1 and same_ratio_result < 0.6:
						category = 13
					elif start2+start4 == 2 and same_ratio_result >= 0.6:
						category = 11
					elif start2+start4 == 2 and same_ratio_result < 0.6:
						category = 12
					##uniq dst and inseq
					
					dict1[line.strip().split('\t')[0]] = "\t".join(line.strip().split('\t')[0:5])+'\t'+cigar+'\t'+"\t".join(line.strip().split('\t')[6:20])+'\t'+str(dst)+'\t'+"\t".join(line.strip().split('\t')[21:24])+'\t'+str(insize)+'\t'+inseq
#				print(line.strip()+'\traw',file=out_kedup)
				print("\t".join(line.strip().split('\t')[0:5])+'\t'+cigar+'\t'+"\t".join(line.strip().split('\t')[6:20])+'\t'+str(dst)+'\t'+"\t".join(line.strip().split('\t')[21:24])+'\t'+str(insize)+'\t'+inseq+'\t'+"\t".join(line.strip().split('\t')[26:28])+'\t'+str(category)+'\t'+before_seq+'\t'+after_seq+'\t'+str(same_ratio_result)+'\t'+str(basefrom)+'\t'+str(baseto),file=out_kedup)
			else:
				print(line.strip(),file=out_kedup)
		else:
			print(line.strip(),file=out_kedup)
	else:
		print(line.strip(),file=out_kedup)



#for line in open("HQ4664_A_Uins2.txt"):
for line in open(sys.argv[3]):
	header = line.strip().split('\t')[0]
	if header in dict1.keys():
#		print(line.strip()+'\tchange_raw',file=out_nokedup)
		print(dict1[header]+'\tchange_new',file=out_nokedup)
	else:
		print(line.strip()+'\tnochange',file=out_nokedup)
