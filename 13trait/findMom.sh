#!/bin/bash

sort phe.bw >phe.bw.sorted
sort phe.ww >phe.ww.sorted
sort stacked_ped >stacked.ped.sorted

awk '{print $1}' phe.bw.sorted > id.bw.sorted
awk '{print $1}' phe.ww.sorted > id.ww.sorted

join id.bw.sorted stacked.ped.sorted| awk '{print $1,$3}' > mom.bw 
join id.ww.sorted stacked.ped.sorted| awk '{print $1,$3}' > mom.ww 
