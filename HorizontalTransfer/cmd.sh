#!/bin/bash


# create directories for the plots

for DIR in 1 2 A B C1 C2 D E ; do mkdir plots_T$DIR ;  done

# make a detailed transfer plot for each clone

for i in `seq 4663 4722`; do R --vanilla --args $i 1  290.0 < plotContigs.R ; done | grep "^\[1\]" > T1.log
for i in `seq 4723 4782`; do R --vanilla --args $i 2  290.0 < plotContigs.R ; done | grep "^\[1\]" > T2.log
for i in `seq 5175 5214`; do R --vanilla --args $i A  290.0 < plotContigs.R ; done | grep "^\[1\]" > TA.log
for i in `seq 5217 5256`; do R --vanilla --args $i B  290.0 < plotContigs.R ; done | grep "^\[1\]" > TB.log
for i in `seq 5257 5296`; do R --vanilla --args $i C1 290.0 < plotContigs.R ; done | grep "^\[1\]" > TC1.log
for i in `seq 5385 5424`; do R --vanilla --args $i C2 290.0 < plotContigs.R ; done | grep "^\[1\]" > TC2.log
for i in `seq 5301 5340`; do R --vanilla --args $i D  290.0 < plotContigs.R ; done | grep "^\[1\]" > TD.log
for i in `seq 5343 5382`; do R --vanilla --args $i E  290.0 < plotContigs.R ; done | grep "^\[1\]" > TE.log

# make a summary plot for each experiment

for DIR in 1 2 A B C1 C2 D E  ; do  R --vanilla --args $DIR < plotContigsSummary.R ; done

# collect transfer boundaries in one file per experiment for further analysis

for i in T1 T2 TA TB TC1 TC2 TD TE ; do cd plots_$i ; grep -h ^sample *.txt | head -1 > ../$i\_details.txt ; grep -h -v ^sample *.txt >> ../$i\_details.txt ; cd .. ; done

# end
