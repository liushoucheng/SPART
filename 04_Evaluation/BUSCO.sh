#!/bin/sh

geno=$1
profix=$2

busco -m geno -i $geno -l poales_odb10 -o $profix -c 52

