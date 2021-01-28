#!/usr/bin/perl -w
use List::Util qw/sum max min/;
use List::MoreUtils;

open IN,"$ARGV[0]";
#open OUT,">$ARGV[1]";

$linenum=0;
while (<IN>) {
	chomp;
	my @a=split(/\s+/);
	if ($linenum==0){
		print "$_\n";
	}else {
	##if the ref is long, then short the seqfrom
	if ($a[3]>10){
		splice(@a,3,1,$a[3]-169);
		splice(@a,20,1,$a[20]-169);
	}

	##use the seq that dst is bigger than 0
	if ($a[20]>0){
        $cigar=$a[5];
        my $countI=($cigar=~s/I/#/g);
        my $countM=($cigar=~s/M/#/g);
        my $countD=($cigar=~s/D/#/g);
        @mid=split/\d+/,$a[5];
        @mid_new=@mid[1..$#mid];
        @num=split/[a-zA-Z]/,$a[5];
        @num_new=@num[0..$#num];

	$pos_sum=$a[3]-1;$len_sum=0;
	for ($i=0;$i<=$#mid_new;$i++){
		if ($mid_new[$i] =~m/M/) {
			$pos_sum=$pos_sum+$num_new[$i];
			$len_sum=$len_sum+$num_new[$i];
			if ($pos_sum>=144){
				$pos_144=144-($pos_sum-$num_new[$i])+($len_sum-$num_new[$i])-1;
                                $str_144=substr($a[9],$pos_144,1);
				if ($str_144 ne "T"){
					print join"\t",@a[0..21];print "\t\t\t$a[22]\t$a[23]\n";
					last;
				}else {
					last;
				}
			}
		}elsif ($mid_new[$i] =~m/D/){
			$pos_sum=$pos_sum+$num_new[$i];
			if ($pos_sum>=144){
				#print join"\t",@a;print "\tD\n";
				print join"\t",@a[0..21];print "\t\t\t$a[22]\t$a[23]\n";
				last;
			}
		}elsif ($mid_new[$i] =~m/I/) {
			$len_sum+=$num_new[$i];
		}
	}
	}
	}
	$linenum++;
}
