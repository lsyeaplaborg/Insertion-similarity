#!/usr/bin/perl -w
use List::Util qw/sum max min/;
use List::MoreUtils;

open IN,"$ARGV[0]";
open TWO,"$ARGV[1]";
open OUT,">$ARGV[2]";

#open TWO,"/Users/guituantuan/Desktop/scripts/hao_qian_script/VB18_F3_short.fa";
local $/=">";
$ref="";
while (<TWO>) {
	chomp;
	my @a=split(/\n/);
	for ($i=1;$i<=$#a;$i++){
              $ref="$ref$a[$i]";
        }	
}

$line=0;
local $/="\n";
while (<IN>) {
	chomp;
	my @array=split(/\t/);
	$line++;

	if ($line<2) {
		print OUT "$_\tnum_I\tnum_D\tcategory\tbefore_seq\tafter_seq\tsame_ratio\tbasefrom\tbaseto\n";
	}else {
		my @array=split(/\t/);
		$cigar=$array[5];
		$ins=$array[20];
		$insize=$array[24];

		@mid=split(/\d+/,$cigar);
                @mid_new=@mid[1..$#mid];
                @num=split(/[a-zA-Z]/,$cigar);
                @num_new=@num[0..$#num];

		$sum=$array[3]-1;@num_new1=();@mid_new1=();
		for ($i=0;$i<=$#num_new;$i++){
			if ($mid_new[$i] ne "I") {
				if ($sum<=201){
					push @num_new1,$num_new[$i];
					push @mid_new1,$mid_new[$i];
				}
				$sum=$sum+$num_new[$i];
			}else {
				if ($sum<=201){
					push @num_new1,$num_new[$i];
	              	        	 push @mid_new1,$mid_new[$i];
				}
			}
	
		}

		$countI=0;$countD=0;
		foreach $a(@mid_new1){
			if ($a eq "I") {
				$countI++;
			}
			if ($a eq "D") {
				$countD++;
			}
		}
		
		##if countD=0 and countI>1,category is 5;if countD=1 and countI>1, or countD>1,category is 6
        if ($countD==0 and $countI>1){
                print OUT "$_\tcountI_$countI\tcountD_$countD\t5\n";
        }elsif ($countD==1 and $countI>1) {
                print OUT "$_\tcountI_$countI\tcountD_$countD\t6\n";
        }elsif ($countD>1){
                print OUT "$_\tcountI_$countI\tcountD_$countD\t6\n";	
		}else {
			@mid_new=@mid_new1;@num_new=@num_new1;

			if ($countD==0) {
				$seq=$array[9];
				$seq_forward=substr($seq,($num_new[0]-$num_new[1]),$num_new[1]);
				$seq_backward=substr($seq,($num_new[0]+$num_new[1]),$num_new[1]);

				$ref_forward=substr($ref,($num_new[0]-$num_new[1]+$array[3]-1),$num_new[1]);
				$ref_backward=substr($ref,($num_new[0]+$array[3]-1),$num_new[1]);
				
				$str_forward_num=0;$str_backward_num=0;
				for ($seq_i=0;$seq_i<length($seq_forward);$seq_i++) {
					$str_forward=substr($seq_forward,$seq_i,1);
					$str_backward=substr($seq_backward,$seq_i,1);
					$str_ins=substr($array[25],$seq_i,1);
					if ($str_forward eq $str_ins) {
						$str_forward_num++;
					}
					if ($str_backward eq $str_ins) {
						$str_backward_num++;
					}
				}

				$str_forward_ref_num=0;$str_backward_ref_num=0;
                                for ($seq_i=0;$seq_i<length($ref_forward);$seq_i++) {
                                        $str_forward=substr($ref_forward,$seq_i,1);
                                        $str_backward=substr($ref_backward,$seq_i,1);
                                        $str_ins=substr($array[25],$seq_i,1);
                                        if ($str_forward eq $str_ins) {
                                                $str_forward_ref_num++;
                                        }
                                        if ($str_backward eq $str_ins) {
                                                $str_backward_ref_num++;
                                        }
                                }
				$ref_forward_ratio=$str_forward_ref_num/$array[24];
				$ref_backward_ratio=$str_backward_ref_num/$array[24];

				$str_forward_ratio=$str_forward_num/$array[24];
				$str_backward_ratio=$str_backward_num/$array[24];

				$str_forward_ratio=max ($ref_forward_ratio,$str_forward_ratio);
				$str_backward_ratio=max ($ref_backward_ratio,$str_backward_ratio);

				if ($num_new[1]==1 or $num_new[1]==2) {
				
					if ($num_new[1]==1){
						$category_1_1=2;
						$category_1_2=13;
					}elsif ($num_new[1]==2){
						$category_1_1=11;
						$category_1_2=12;
					}
					if ($str_forward_ratio>=1 or $str_backward_ratio>=1) {
						if ($str_forward_ratio>$str_backward_ratio) {
							$baseto=$ins-1;
							$basefrom=$baseto-$insize+1;
							print OUT "$_\t$num_new[1]\t0\t$category_1_1\t$seq_forward\t$seq_backward\t$str_forward_ratio\t$basefrom\t$baseto\n";
						}else {
							$basefrom=$ins;
							$baseto=$basefrom+$insize-1;
							print OUT "$_\t$num_new[1]\t0\t$category_1_1\t$seq_forward\t$seq_backward\t$str_backward_ratio\t$basefrom\t$baseto\n";
						}
					}else {
						if ($str_forward_ratio>$str_backward_ratio) {
							$baseto=$ins-1;
                                                        $basefrom=$baseto-$insize+1;
                                                        print OUT "$_\t$num_new[1]\t0\t$category_1_2\t$seq_forward\t$seq_backward\t$str_forward_ratio\n";
                                                }else {
							$basefrom=$ins;
                                                        $baseto=$basefrom+$insize-1;
                                                        print OUT "$_\t$num_new[1]\t0\t$category_1_2\t$seq_forward\t$seq_backward\t$str_backward_ratio\n";
                                                }
					}
				}elsif ($num_new[1]==3) {
					if ($str_forward_ratio>=0.6 or $str_backward_ratio>=0.6) {
						if ($str_forward_ratio>$str_backward_ratio) {
							$baseto=$ins-1;
                                                        $basefrom=$baseto-$insize+1;
                                                        print OUT "$_\t$num_new[1]\t0\t3\t$seq_forward\t$seq_backward\t$str_forward_ratio\t$basefrom\t$baseto\n";
                                                }else {
							$basefrom=$ins;
                                                        $baseto=$basefrom+$insize-1;
                                                        print OUT "$_\t$num_new[1]\t0\t3\t$seq_forward\t$seq_backward\t$str_backward_ratio\t$basefrom\t$baseto\n";
                                                }
                                        }else {
						if ($str_forward_ratio>$str_backward_ratio) {
							$baseto=$ins-1;
                                                        $basefrom=$baseto-$insize+1;
                                                        print OUT "$_\t$num_new[1]\t0\t7\t$seq_forward\t$seq_backward\t$str_forward_ratio\n";
                                                }else {
							$basefrom=$ins;
                                                        $baseto=$basefrom+$insize-1;
                                                        print OUT "$_\t$num_new[1]\t0\t7\t$seq_forward\t$seq_backward\t$str_backward_ratio\n";
                                                }
                                        }
				}elsif ($num_new[1]>3) {
					if ($str_forward_ratio>=0.6 or $str_backward_ratio>=0.6) {
                                        	if ($str_forward_ratio>$str_backward_ratio) {
							$baseto=$ins-1;
                                                        $basefrom=$baseto-$insize+1;
                                                        print OUT "$_\t$num_new[1]\t0\t1\t$seq_forward\t$seq_backward\t$str_forward_ratio\t$basefrom\t$baseto\n";
                                                }else {
							$basefrom=$ins;
                                                        $baseto=$basefrom+$insize-1;
                                                        print OUT "$_\t$num_new[1]\t0\t1\t$seq_forward\t$seq_backward\t$str_backward_ratio\t$basefrom\t$baseto\n";
                                                }
                                        }else {
						if ($str_forward_ratio>$str_backward_ratio) {
							$baseto=$ins-1;
                                                        $basefrom=$baseto-$insize+1;
                                                        print OUT "$_\t$num_new[1]\t0\t7\t$seq_forward\t$seq_backward\t$str_forward_ratio\n";
                                                }else {
							$basefrom=$ins;
                                                        $baseto=$basefrom+$insize-1;
                                                        print OUT "$_\t$num_new[1]\t0\t7\t$seq_forward\t$seq_backward\t$str_backward_ratio\n";
                                                } 
                                        }
				}
			}else {
				$index_I=List::MoreUtils::firstidx{/I/} @mid_new;
				$index_D=List::MoreUtils::firstidx{/D/} @mid_new;
				if ($num_new[$index_I] == $num_new[$index_D] and $num_new[2]<=5){
					print OUT "$_\t$num_new[$index_I]\t$num_new[$index_D]\t0\n";
				}elsif ($num_new[$index_I] == $num_new[$index_D] and $num_new[2]>5){
					print OUT "$_\t$num_new[$index_I]\t$num_new[$index_D]\t4\n";
				}elsif ($num_new[$index_I] != $num_new[$index_D]) {
					if ($num_new[2]>5) {
						print OUT "$_\t$num_new[$index_I]\t$num_new[$index_D]\t8\n";
					}else {
						if ($num_new[$index_I]>$num_new[$index_D]) {
							print OUT "$_\t$num_new[$index_I]\t$num_new[$index_D]\t9\n";
						}elsif ($num_new[$index_I]<$num_new[$index_D]) {
							print OUT "$_\t$num_new[$index_I]\t$num_new[$index_D]\t10\n";
						}
					}
				}
			}

		}
	}
}

