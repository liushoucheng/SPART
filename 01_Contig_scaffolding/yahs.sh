#!/bin/sh

enzyme=$1
ref=$2
bed=$3
profix=$4

yahs -e $enzyme $ref $bed -o $4

