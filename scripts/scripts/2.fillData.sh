#!/bin/bash

OLDDATA="/home/haocheng/Hereford/13Traits_small/data_original" #original data, pedigree with A inverse, G matrix
NEWDATA="/home/haocheng/Hereford/13Traits_small/data_new" #modified data

cd $OLDDATA

awk '{print $1}' phe100.bw phe100.ww phe100.ss phe100.yw phe100.fw phe100.carc_ema phe100.carc_imf phe100.carc_rib phe100.cwt phe100.ema phe100.imf phe100.rib phe100.mcw.single | sort -u > IDs.sorted

sort phe100.bw >phe100.bw.sorted
sort phe100.ww >phe100.ww.sorted
sort phe100.yw >phe100.yw.sorted
sort phe100.fw >phe100.fw.sorted
sort phe100.carc_ema >phe100.carc_ema.sorted
sort phe100.carc_imf >phe100.carc_imf.sorted
sort phe100.carc_rib >phe100.carc_rib.sorted
sort phe100.cwt >phe100.cwt.sorted
sort phe100.ema >phe100.ema.sorted
sort phe100.imf >phe100.imf.sorted
sort phe100.rib >phe100.rib.sorted
sort phe100.ss >phe100.ss.sorted
sort phe100.mcw.single >phe100.mcw.single.sorted

join -a 1 -a 2 phe100.yw.sorted IDs.sorted > phe100.yw.fill_IDs
join -a 1 -a 2 phe100.ww.sorted IDs.sorted > phe100.ww.fill_IDs
join -a 1 -a 2 phe100.fw.sorted IDs.sorted > phe100.fw.fill_IDs
join -a 1 -a 2 phe100.bw.sorted IDs.sorted > phe100.bw.fill_IDs
join -a 1 -a 2 phe100.carc_ema.sorted IDs.sorted > phe100.carc_ema.fill_IDs
join -a 1 -a 2 phe100.carc_imf.sorted IDs.sorted > phe100.carc_imf.fill_IDs
join -a 1 -a 2 phe100.carc_rib.sorted IDs.sorted > phe100.carc_rib.fill_IDs
join -a 1 -a 2 phe100.cwt.sorted IDs.sorted > phe100.cwt.fill_IDs
join -a 1 -a 2 phe100.ema.sorted IDs.sorted > phe100.ema.fill_IDs
join -a 1 -a 2 phe100.imf.sorted IDs.sorted > phe100.imf.fill_IDs
join -a 1 -a 2 phe100.rib.sorted IDs.sorted > phe100.rib.fill_IDs
join -a 1 -a 2 phe100.mcw.single.sorted IDs.sorted > phe100.mcw.single.fill_IDs
join -a 1 -a 2 phe100.ss.sorted IDs.sorted > phe100.ss.fill_IDs

awk '{for(i=1;i<=4;i++) if($i=="") $i="\."}1' phe100.yw.fill_IDs > phe100.yw.fill
awk '{for(i=1;i<=4;i++) if($i=="") $i="\."}1' phe100.ww.fill_IDs > phe100.ww.fill
awk '{for(i=1;i<=4;i++) if($i=="") $i="\."}1' phe100.fw.fill_IDs > phe100.fw.fill
awk '{for(i=1;i<=4;i++) if($i=="") $i="\."}1' phe100.bw.fill_IDs > phe100.bw.fill
awk '{for(i=1;i<=4;i++) if($i=="") $i="\."}1' phe100.carc_ema.fill_IDs > phe100.carc_ema.fill
awk '{for(i=1;i<=4;i++) if($i=="") $i="\."}1' phe100.carc_imf.fill_IDs > phe100.carc_imf.fill
awk '{for(i=1;i<=4;i++) if($i=="") $i="\."}1' phe100.carc_rib.fill_IDs > phe100.carc_rib.fill
awk '{for(i=1;i<=4;i++) if($i=="") $i="\."}1' phe100.cwt.fill_IDs > phe100.cwt.fill
awk '{for(i=1;i<=4;i++) if($i=="") $i="\."}1' phe100.ema.fill_IDs > phe100.ema.fill
awk '{for(i=1;i<=4;i++) if($i=="") $i="\."}1' phe100.imf.fill_IDs > phe100.imf.fill
awk '{for(i=1;i<=4;i++) if($i=="") $i="\."}1' phe100.rib.fill_IDs > phe100.rib.fill
awk '{for(i=1;i<=4;i++) if($i=="") $i="\."}1' phe100.mcw.single.fill_IDs > phe100.mcw.single.fill
awk '{for(i=1;i<=4;i++) if($i=="") $i="\."}1' phe100.ss.fill_IDs > phe100.ss.fill


#add pedgree to the last two columns
sort $NEWDATA"/stacked_ped" > peds
join phe100.bw.fill peds > phe100.bw.fill.ped
join phe100.yw.fill peds > phe100.yw.fill.ped
join phe100.ww.fill peds > phe100.ww.fill.ped
join phe100.fw.fill peds > phe100.fw.fill.ped
join phe100.carc_ema.fill peds > phe100.carc_ema.fill.ped
join phe100.carc_imf.fill peds > phe100.carc_imf.fill.ped
join phe100.carc_rib.fill peds > phe100.carc_rib.fill.ped
join phe100.cwt.fill peds > phe100.cwt.fill.ped
join phe100.ema.fill peds > phe100.ema.fill.ped
join phe100.imf.fill peds > phe100.imf.fill.ped
join phe100.rib.fill peds > phe100.rib.fill.ped
join phe100.mcw.single.fill peds > phe100.mcw.fill.ped
join phe100.ss.fill peds > phe100.ss.fill.ped

rm phe100*sorted
rm phe100*fill_IDs
rm phe100*fill
rm peds
rm IDs.sorted

mv phe100.bw.fill.ped phe100.yw.fill.ped phe100.ww.fill.ped phe100.fw.fill.ped phe100.carc_ema.fill.ped phe100.carc_imf.fill.ped phe100.carc_rib.fill.ped phe100.cwt.fill.ped phe100.ema.fill.ped phe100.imf.fill.ped phe100.rib.fill.ped phe100.mcw.fill.ped phe100.ss.fill.ped $NEWDATA

