#!/usr/bin/perl -w

use List::Util qw/sum max min/;
use List::MoreUtils;

open IN,"$ARGV[0]";
open OUT1,">$ARGV[1]";
open OUT2,">$ARGV[2]";

$line=0;
while (<IN>) {
	chomp;
	my @array=split(/\s+/);
	if ($line>0) {
		$cigar=$array[5];
#                $ins=$array[20];
 #               $insize=$array[24];

                @mid=split(/\d+/,$cigar);
                @mid_new=@mid[1..$#mid];
                @num=split(/[a-zA-Z]/,$cigar);
                @num_new=@num[0..$#num];

		$pos_sum=$array[3]-1;$i=0;
		while ($i<=$#num_new){
			if ($mid_new[$i] ne "I") {
				$pos_sum+=$num_new[$i];
				if ($pos_sum>201){
					if ($mid_new[$i+1] eq "D" and $num_new[$i+1] eq 12){
						print OUT1 "$_\n";
						last;
					}else {
						$i++;
					}
				}else {
					$i++;
				}
			}else {
				$i++;
			}
		}
		$b=$#num_new+1;
		if ($i==$b){
			print OUT2 "$_\n";
		}
	}else {
		print OUT1 "$_\n";
		print OUT2 "$_\n";
	}
	$line++;
}
close IN;
