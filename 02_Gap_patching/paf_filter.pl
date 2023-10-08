#!/usr/bin/perl
#author         :SHOUCHENG LIU
#email          :286596224@QQ.com
#last modified  :

#-------------------------------------------------------------------------
#please consult DOCUMENTATION for detailed information...
#perl paf_filter.pl -i ragtag.patch.debug.filtered.paf -minlen 500000 -iden 0.8
use strict;
use diagnostics;
use warnings;
use Getopt::Long;

my %opts=();
GetOptions(\%opts,"i:s","minlen:s","iden:s");

my %h=();
open IN, $opts{i} or die "Cannot open file: $opts{i}!\n";
my $o=$opts{i}."_fiter.paf";
open OUT,">$o" or die "Cannot create file: $o!\n";
#ragtag.patch.debug.filtered.paf
while(<IN>){
 chomp;
  my $number=0;
  my @f = split /\t/, $_;
  if ($f[8]>($f[6]-$opts{minlen}) or $f[7]<$opts{minlen}) {
	  print $f[8]."\t".$f[9]/$f[10]."\n";
	if ($f[9]/$f[10]>$opts{iden}) {
$h{$f[0]}{$f[5]}=1
	}}
  }
close IN;
open IN, $opts{i} or die "Cannot open file: $opts{i}!\n";
while(<IN>){
 chomp;
  my $number=0;
  my @f = split /\t/, $_;
  if ($f[8]>($f[6]-$opts{minlen}) or $f[7]<$opts{minlen}) {
	  print $f[8]."\t".$f[9]/$f[10]."\n";
	if ($f[9]/$f[10]>$opts{iden}) {
		my $le=keys %{$h{$f[0]}};
		print $le."\n";
		if ($le>1) {

			print OUT $_."\n";
		}
	}
  }
  }