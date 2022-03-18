#!/bin/bash

raw_fq1=
raw_fq2=
strain=
adapters=
phix=

#trim adapters (shouldn't be any?)
bbduk.sh in1=$raw_fq1 in2=$raw_fq2 out1=${strain}_trim1.fq out2=${strain}_trim2.fq ref=$adapters ktrim=r k=23 mink=11 hdist=1 tpe tbo

#contaminant (phix) filtering
bbduk.sh in1=${strain}_trim1.fq in2=${strain}_trim2.fq out1=${strain}_unmatched1.fq out2=${strain}_unmatched2.fq outm1=${strain}_matched1.fq outm2=${strain}_matched2.fq ref=$phix k=31 hdist=1 stats=phistats.txt

#quality trimming (bbduk user guide recommends this as separate step from adapter trimming)
bbduk.sh in1=${strain}_unmatched1.fq in2=${strain}_unmatched2.fq out1=${strain}_bbduk1.fq out2=${strain}_bbduk2.fq qtrim=rl trimq=10

#ready for alignment!
