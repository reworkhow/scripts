#!/bin/bash

OLDDATA="/home/haocheng/Hereford/13Traits_small/data_original" #original data, pedigree with A inverse, G matrix
NEWDATA="/home/haocheng/Hereford/13Traits_small/data_new" #modified data

#replace spaces in original dataset with .
sed s/,,/,.,/ $OLDDATA"/ped100.csv" | sed s/,$/,./ | sed 's/,/ /g' > ped1

#1st run to generate nopedrec;animals with no parents
#"Sires and dams without pedigree record found.  Placed in file nopedrec."
stack_ped ped1 stacked_ped

#add two extra column of .
awk '{$(NF+1)=".";}1' nopedrec | awk '{$(NF+1)=".";}1' > nopedrec1

#paste vertically
cat nopedrec1 ped1 > ped_done
stack_ped ped_done stacked_ped

invnrm -a -v Ainverse -i stacked_ped -o inbreedings 

mv stacked_ped Ainverse inbreedings $NEWDATA
rm ped1  ped_done nopedrec1
