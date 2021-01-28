#!/usr/bin/perl -w

open IN,"$ARGV[0]";
open OUT,">$ARGV[1]";

my %count;
$line=0;
while (<IN>) {
	chomp;
	if ($line>0) {
		my @a=split(/\t/);
		$count{$a[28]}++;
	}
	$line++;
}
close IN;

@num=(0..13);
@key=keys %count;
#my @a = (1,2,3,4,5,6,7,8,);
#my @b = (1,9,0,4,15,6,12,8);
my %hash_num = map{$_=>1} @num;
my %hash_key = map{$_=>1} @key;
my @num_only = grep {!$hash_key{$_}} @num;
	
foreach $num_only(@num_only){
	$count{$num_only}=0;
}

print OUT"\t$ARGV[0]\n";
foreach $key (sort {$a<=>$b} keys %count) {
                        print OUT "$key\t$count{$key}\n";
}

$sum=0;
foreach $key (keys %count) {
        if ($key==1 or $key==3 or $key==11){
                $sum=$sum+$count{$key};
        }
}

$sum_7_12=$count{7}+$count{12};
$dup_all=$sum+$count{2};
$nondup_all=$sum_7_12+$count{13};
$total_ins=$sum+$count{2}+$sum_7_12+$count{13}+$count{5};

print OUT "dup_1bp\t$count{2}\n";
print OUT "dup>1bp\t$sum\n";
print OUT "nondup_1bp\t$count{13}\n";
print OUT "nondup>1bp\t$sum_7_12\n";


print OUT "dup_all\t$dup_all\n";
print OUT "nondup_all\t$nondup_all\n";
print OUT "gt1I\t$count{5}\n";

print OUT "total_ins\t$total_ins\n";