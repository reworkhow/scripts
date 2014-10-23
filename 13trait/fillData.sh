#!/bin/bash

awk '{print $1}' phe.bw phe.ww phe.ss phe.yw phe.fw phe.carc_ema phe.carc_imf phe.carc_rib phe.cwt phe.ema phe.imf phe.rib phe.mcw.single | sort -u > IDs.sorted

sort phe.bw >phe.bw.sorted
sort phe.ww >phe.ww.sorted
sort phe.yw >phe.yw.sorted
sort phe.fw >phe.fw.sorted
sort phe.carc_ema >phe.carc_ema.sorted
sort phe.carc_imf >phe.carc_imf.sorted
sort phe.carc_rib >phe.carc_rib.sorted
sort phe.cwt >phe.cwt.sorted
sort phe.ema >phe.ema.sorted
sort phe.imf >phe.imf.sorted
sort phe.rib >phe.rib.sorted
sort phe.ss.single >phe.ss.single.sorted
sort phe.mcw.single >phe.mcw.single.sorted

join -a 1 -a 2 phe.yw.sorted IDs.sorted > phe.yw.fill_IDs
join -a 1 -a 2 phe.ww.sorted IDs.sorted > phe.ww.fill_IDs
join -a 1 -a 2 phe.fw.sorted IDs.sorted > phe.fw.fill_IDs
join -a 1 -a 2 phe.bw.sorted IDs.sorted > phe.bw.fill_IDs
join -a 1 -a 2 phe.carc_ema.sorted IDs.sorted > phe.carc_ema.fill_IDs
join -a 1 -a 2 phe.carc_imf.sorted IDs.sorted > phe.carc_imf.fill_IDs
join -a 1 -a 2 phe.carc_rib.sorted IDs.sorted > phe.carc_rib.fill_IDs
join -a 1 -a 2 phe.cwt.sorted IDs.sorted > phe.cwt.fill_IDs
join -a 1 -a 2 phe.ema.sorted IDs.sorted > phe.ema.fill_IDs
join -a 1 -a 2 phe.imf.sorted IDs.sorted > phe.imf.fill_IDs
join -a 1 -a 2 phe.rib.sorted IDs.sorted > phe.rib.fill_IDs
join -a 1 -a 2 phe.mcw.single.sorted IDs.sorted > phe.mcw.single.fill_IDs
join -a 1 -a 2 phe.ss.single.sorted IDs.sorted > phe.ss.single.fill_IDs

awk '{for(i=1;i<=4;i++) if($i=="") $i="\."}1' phe.yw.fill_IDs > phe.yw.fill
awk '{for(i=1;i<=4;i++) if($i=="") $i="\."}1' phe.ww.fill_IDs > phe.ww.fill
awk '{for(i=1;i<=4;i++) if($i=="") $i="\."}1' phe.fw.fill_IDs > phe.fw.fill
awk '{for(i=1;i<=4;i++) if($i=="") $i="\."}1' phe.bw.fill_IDs > phe.bw.fill
awk '{for(i=1;i<=4;i++) if($i=="") $i="\."}1' phe.carc_ema.fill_IDs > phe.carc_ema.fill
awk '{for(i=1;i<=4;i++) if($i=="") $i="\."}1' phe.carc_imf.fill_IDs > phe.carc_imf.fill
awk '{for(i=1;i<=4;i++) if($i=="") $i="\."}1' phe.carc_rib.fill_IDs > phe.carc_rib.fill
awk '{for(i=1;i<=4;i++) if($i=="") $i="\."}1' phe.cwt.fill_IDs > phe.cwt.fill
awk '{for(i=1;i<=4;i++) if($i=="") $i="\."}1' phe.ema.fill_IDs > phe.ema.fill
awk '{for(i=1;i<=4;i++) if($i=="") $i="\."}1' phe.imf.fill_IDs > phe.imf.fill
awk '{for(i=1;i<=4;i++) if($i=="") $i="\."}1' phe.rib.fill_IDs > phe.rib.fill
awk '{for(i=1;i<=4;i++) if($i=="") $i="\."}1' phe.mcw.single.fill_IDs > phe.mcw.single.fill
awk '{for(i=1;i<=4;i++) if($i=="") $i="\."}1' phe.ss.single.fill_IDs > phe.ss.single.fill


#add pedgree to the last two columns
sort stacked_ped > peds
join phe.bw.fill peds > phe.bw.fill.ped
join phe.yw.fill peds > phe.yw.fill.ped
join phe.ww.fill peds > phe.ww.fill.ped
join phe.fw.fill peds > phe.fw.fill.ped
join phe.carc_ema.fill peds > phe.carc_ema.fill.ped
join phe.carc_imf.fill peds > phe.carc_imf.fill.ped
join phe.carc_rib.fill peds > phe.carc_rib.fill.ped
join phe.cwt.fill peds > phe.cwt.fill.ped
join phe.ema.fill peds > phe.ema.fill.ped
join phe.imf.fill peds > phe.imf.fill.ped
join phe.rib.fill peds > phe.rib.fill.ped
join phe.mcw.single.fill peds > phe.mcw.fill.ped
join phe.ss.single.fill peds > phe.ss.fill.ped

rm phe*sorted
rm phe*fill_IDs
rm phe*fill
rm peds
rm IDs.sorted
