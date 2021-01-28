#python3
import sys

dict_uniq={}
dict_all={}
#for line in open("HQ4664_A_Uins2_new.txt"):
for line in open(sys.argv[1]):
	dst=line.strip().split('\t')[20]
	inseq=line.strip().split('\t')[25]
	t=line.strip().split('\t')[26]
	if t == "nochange":
		dict_all[dst+'\t'+inseq]=1
		print(line.strip())
	else:
		dict_uniq[dst+'\t'+inseq]=line.strip()


for k in dict_uniq.keys():
	if k not in dict_all.keys():
		print(dict_uniq[k])
#	else:
#		print(dict_uniq[k]+'\tfilter_out')
