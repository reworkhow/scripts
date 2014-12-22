#!/bin/bash

#reform data
sed s/,,/,.,/ ped.csv > ped1
sed s/,$/,./ ped1 > ped2
sed 's/,/ /g' ped2 > ped3

#1st run to generate nopedrec;animals with no parents
#"Sires and dams without pedigree record found.  Placed in file nopedrec."
stack_ped ped3 stacked_ped

#add one extra column of .
awk '{$(NF+1)=".";}1' nopedrec > nopedrec1

#add one extra column of .
awk '{$(NF+1)=".";}1' nopedrec1 > nopedrec2

#paste vertically
cat nopedrec2 ped3 | paste > ped_done
stack_ped ped_done stacked_ped

invnrm -a -v Ainverse -i stacked_ped -o inbreedings 

rm ped1 ped2 ped3 ped_done
rm nopedrec1 nopedrec2
