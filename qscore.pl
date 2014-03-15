#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my $cl;
my $pl;
my $statsfile;

GetOptions("clustering=s"=>\$cl,"placements=s"=>\$pl, "statsfile=s"=>\$statsfile);

use POSIX;
my %h;
open(IN, "<$cl");
my $clid = 0;
while(<IN>){
$clid++;
chomp;
my @arr=split();
# For each contig, record the cluster number it was placed in
for(my $i=0; $i<@arr; $i++){ 
$h{$arr[$i]}{cluster_id} = $clid;
}
}
close(IN);

open(IN, "<$pl");
while(<IN>){
chomp;
my @arr=split();
# Record the id of where each contig should be placed
$h{$arr[0]}{org_id} = $arr[1];
}
close(IN);

# How many of those placed together belong together
my $tp = 0;

# How many of those separated should be separated
my $tn = 0;

# How many of those separated should not be separated
my $fn = 0;

# How many of those separated should be separated.
my $fp = 0;

my $numtosample=100000;
my @keyarr = (keys %h);
my $len = @keyarr;

for(my $i=0; $i<$numtosample; $i++){
# Randomly choose a pair
# Determine whether they are placed together.
# Determine whether they should be placed together.
my $val1 = floor(rand($len));
my $val2 = floor(rand($len));
my $id1 = $keyarr[$val1];
my $id2 = $keyarr[$val2];

if(defined($h{$id1}{cluster_id})&&defined($h{$id2}{cluster_id})&&defined($h{$id1}{org_id})&&defined($h{$id2}{org_id})){
if(($h{$id1}{cluster_id}==$h{$id2}{cluster_id})&&($h{$id1}{org_id}==$h{$id2}{org_id})){

$tp++;

}elsif(($h{$id1}{cluster_id}!=$h{$id2}{cluster_id})&&($h{$id1}{org_id}!=$h{$id2}{org_id})){

$tn++;

}elsif(($h{$id1}{cluster_id}!=$h{$id2}{cluster_id})&&($h{$id1}{org_id}==$h{$id2}{org_id})){

$fn++;

}elsif(($h{$id1}{cluster_id}==$h{$id2}{cluster_id})&&($h{$id1}{org_id}!=$h{$id2}{org_id})){

$fp++;

}

}

}

my $tpr = 0;
my $fpr = 0;
my $ppv = 0;
my $npv = 0;

if(($tp+$fn)!=0){$tpr = $tp/($tp+$fn);}
if(($fp+$tn)!=0){$fpr = $fp/($fp+$tn);}
if(($tp+$fp)!=0){$ppv = $tp/($tp+$fp);}
if(($tn+$fn)!=0){$npv = $tn/($tn+$fn);}

open(STATS, ">$statsfile");

print STATS "TPR: $tpr\n";
print STATS "FPR: $fpr\n";
print STATS "PPV: $ppv\n";
print STATS "NPV: $npv\n";

close(STATS);





