#!/usr/bin/perl
#author         :SHOUCHENG LIU
#email          :286596224@QQ.com
#last modified  :

#-------------------------------------------------------------------------
#please consult DOCUMENTATION for detailed information...
use strict;
use diagnostics;
use warnings;
use Getopt::Long;
# perl /home/liusc/proj/wheat/code/agp.pl -i ragtag.patch.ctg.agp -i1 ragtag.patch.debug.filtered.paf -start seq00000000 -end seq00000001 -o test.agp
my %opts=();
GetOptions(\%opts,"i:s","i1:s","o:s","start:s","end:s");

if (!$opts{i} or !$opts{o}){
    print "----------------------------------------------------------------------
    USAGE: perl $0
        -i input  file
        -o out file
----------------------------------------------------------------------\n";
    exit;
    
}
my %h=();
open IN, $opts{i} or die "Cannot open file: $opts{i}!\n";
open IN1, $opts{i1} or die "Cannot open file: $opts{i1}!\n";
open OUT,">$opts{o}" or die "Cannot create file: $opts{o}!\n";
L:while(<IN>){
#chr8	1	1000000	1	W	seq00000000	1	1000000	+
#chr8	1000001	1000100	2	N	100	scaffold	yes	align_genus
#chr8	1000101	4999101	3	W	seq00000001	1	3999001	+
#scf0	1000001	1000999	2	W	qseq00000000	10002	11000	+
 chomp;
 if ($_=~ $opts{start}) {
	  my @f = split /\t/, $_;

  $h{$f[5]}=[@f];
  @f = split /\t/, <IN>;
  $h{"scaffold"}=[@f];
  @f = split /\t/, <IN>;
  $h{$f[5]}=[@f];
  last L;
 }
  }
    my $start;
    my $end;
	my $len2;
while(<IN1>){
#qseq00000269	12474565	2150891	12471563	+	seq00000001	79845467	259191	10579680	9480999	10320756	60
#qseq00000269	12474565	75	2194153	+	seq00000000	673241388	671105857	673241383	2029712	2194088	60
my $l1=chomp($_);
  my $number=0;
  my @f = split /\t/, $_;

  if ($f[4] eq "+") {
	  if ($f[5] eq $opts{start}) {
		my $offset=$f[6]-$f[8];
		$h{$f[5]}[2]-=$offset;
		$h{$f[5]}[7]-=$offset;
		
		$h{"scaffold"}[1]-=$offset;
		$start=$f[3];
		$h{"scaffold"}[1]=$h{$opts{start}}[2]+1;
		$h{"scaffold"}[6]=$f[3]+1;
	  }elsif ($f[5] eq $opts{end}){
		$end=$f[2];
		  $h{"scaffold"}[4]="W";
		  $h{"scaffold"}[5]=$f[0];
		  
		  $h{"scaffold"}[8]="+";
		  $h{$opts{end}}[6]=$f[7]+1;
	  }
	  
  }elsif($f[4] eq "-") {
	  if ($f[5] eq $opts{start}) {
		  my $offset=$f[6]-$f[8];
		  $h{$f[5]}[2]-=$offset;
		  $h{$f[5]}[7]-=$offset;
		  $h{"scaffold"}[1]-=$offset;
		  $end=$f[2];
		  $h{"scaffold"}[1]=$h{$opts{start}}[2]+1;
		  $h{"scaffold"}[7]=$f[2];


  }elsif ($f[5] eq $opts{end}){
		$start=$f[3];
		  $h{"scaffold"}[4]="W";
		  $h{"scaffold"}[5]=$f[0];
		  
		  $h{"scaffold"}[8]="-";
		  $h{$f[5]}[6]=$f[7]+1;
		  $h{"scaffold"}[6]=$f[3]+1;
	  }
  }
}
  		  my $gapsize=$end-$start;
		  print $gapsize."\n";
		  $h{"scaffold"}[2]=$h{"scaffold"}[1]+$gapsize-1;
		  $h{"scaffold"}[7]=$h{"scaffold"}[6]+$gapsize-1;
if ($h{"scaffold"}[2]>$h{$opts{start}}[2]) {

		$h{$opts{end}}[1]=$h{"scaffold"}[2]+1;
		$h{$opts{end}}[2]=$h{$opts{end}}[1]+($h{$opts{end}}[7]-$h{$opts{end}}[6]);
		my $string = join "\t", @{$h{$opts{start}}};
		print OUT "## agp-version 2.1\n";
		print OUT "# AGP created by RagTag v2.1.0\n";
		print OUT $string."\n";
		$string = join "\t", @{$h{"scaffold"}};
		print OUT $string."\n";
		$string = join "\t", @{$h{$opts{end}}};
		print OUT $string;
}else{
		$h{$opts{end}}[1]=$h{$opts{start}}[2]+1;
		$h{$opts{end}}[2]=$h{$opts{end}}[1]+($h{$opts{end}}[7]-$h{$opts{end}}[6]);
		$h{$opts{end}}[3]--;
		my $string = join "\t", @{$h{$opts{start}}};
		print OUT "## agp-version 2.1\n";
		print OUT "# AGP created by RagTag v2.1.0\n";
		print OUT $string."\n";
		$string = join "\t", @{$h{$opts{end}}};
		print OUT $string;
		}