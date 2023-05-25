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
# perl renameagp.pl -i ragtag.patch.rename.agp -i1 ragtag.patch.ctg.agp -i2 gap.agp -i3 ragtag.patch.agp -o gap.rename.agp
my %opts=();
GetOptions(\%opts,"i:s","i1:s","o:s","i2:s","i3:s");

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
open IN2, $opts{i2} or die "Cannot open file: $opts{i2}!\n";
open IN3, $opts{i3} or die "Cannot open file: $opts{i3}!\n";
open OUT,">$opts{o}" or die "Cannot create file: $opts{o}!\n";
while(<IN>){
#qseq00000000    1       209209621       1       W       contig_10       1       209209621       +
 chomp;
	  my @f = split /\t/, $_;

  $h{$f[0]}=$f[5];
 }
my $n=0;
my %h2;
while(<IN1>){
#chr1A_RagTag_MOD_MOD    1       2046621 1       W       seq00000001     1       2046621 +
 chomp;
	  my @f = split /\t/, $_;
    if ($f[0]=~/hr/){

  if (!exists $h2{$f[0]}){
                $h{$f[5]}=$f[0];

                $h2{$f[0]}=$f[5];

        }

    }

 }
while(<IN3>){
#scf00000000     1       2046621 1       W       seq00000001     1       2046621 +
 chomp;
	  my @f = split /\t/, $_;
    if (exists $h2{$h{$f[5]}}){
         $h{$f[0]}=$h{$f[5]};
          print $h{$f[5]}."\n";

    }

 }
while(<IN2>){
#scf00000001     457330265       457811306       10      W       qseq00000182
 chomp;
	  my @f = split /\t/, $_;

  $f[0]=$h{$f[0]};
     $f[5]=$h{$f[5]};
    my $string = join "\t", @f;
    print OUT $string."\n";
 }
