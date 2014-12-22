#!/bin/bash

sort -n phe.mcw > phe.mcw.sort
awk '{if(a[$1]) a[$1]=a[$1]" "$0; else a[$1]=$0;}END{for(i in a) print a[i];}' phe.mcw.sort > phe.mcw.done

sort -n phe.ss > phe.ss.sort
awk '{if(b[$1]) b[$1]=b[$1]" "$0; else b[$1]=$0;}END{for(i in b) print b[i];}' phe.ss.sort > phe.ss.done

#now just keep the first obervation for all individuals
awk '{print $1,$2,$3,$4}' phe.mcw.done > phe.mcw.single
awk '{print $1,$2,$3,$4}' phe.ss.done > phe.ss.single

rm phe.mcw.sort phe.mcw.done
rm phe.ss.sort phe.ss.done
