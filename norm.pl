#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my $in;
my $sizes;


GetOptions("in=s"=>\$in, "sizes=s"=>\$sizes);

open(IN, "<$sizes"); 

my $szmax=0;

while(<IN>){
chomp; 
my @arr=split(/\t/, $_); $h{$arr[0]}=$arr[1]; 
if($arr[1]>$szmax){$szmax=$arr[1]}
} 
close(IN); 

open(IN, "<$in");
while(<IN>){
chomp; 
my @arr=split(); 
if(defined($h{$arr[0]})&&defined($h{$arr[1]})){
my $val = ($szmax*$szmax)*$arr[2]/($h{$arr[0]}*$h{$arr[1]}); 
print “$_\t$val\n";
}
}
close(IN);





