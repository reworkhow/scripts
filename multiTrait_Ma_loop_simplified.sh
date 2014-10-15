#!/bin/bash


##########################################
###make matrix for trait with maternal effects
##########################################
array=(bw ww ss)

for i in "${array[@]}"
do

#####make X and Z matrix
awk '{print $2}' phe.$i".fill.ped" > y.$i &
awk '{print $1}' phe.$i".fill.ped" > id.dat.$i &
awk '{print $1}' stacked_ped > id.eff &
awk '{print $4}' phe.$i".fill.ped" > cg1.dat.$i &
awk '{print $4}' phe.$i".fill.ped" | sort -u > cg1.eff.$i &
wait
cgen_z -d cg1.dat.$i -e cg1.eff.$i -r y.$i -o X.$i
cgen_z -d id.dat.$i -e id.eff -r y.$i -o Z.$i
wait

#####maternal effect
awk '{print $1}' stacked_ped > mat.eff.$i &
awk '{print $6}' phe.bw.fill.ped > mat.dat.$i &
wait
cgen_z -d mat.dat.$i -e mat.eff.$i -r y.$i -o M.$i

######permanent env effect
awk '$6!="\."{print $6}' phe.$i".fill.ped" | sort -u > perm.eff.$i &
awk '{print $6}' phe.$i".fill.ped" > perm.dat.$i &
wait
cgen_z -d perm.dat.$i -e perm.eff.$i -r y.$i -o P.$i

done

###########################################
###make matrix for trait without maternal effects
###########################################

array=(fw yw carc_ema carc_imf carc_rib cwt ema imf rib)

for i in "${array[@]}"
do

#####make X and Z matrix
awk '{print $2}' phe.$i".fill.ped" > y.$i &
awk '{print $1}' phe.$i".fill.ped" > id.dat.$i &
awk '{print $1}' stacked_ped > id.eff &
awk '{print $4}' phe.$i".fill.ped" > cg1.dat.$i &
awk '{print $4}' phe.$i".fill.ped" | sort -u > cg1.eff.$i &
wait
cgen_z -d cg1.dat.$i -e cg1.eff.$i -r y.$i -o X.$i
cgen_z -d id.dat.$i -e id.eff -r y.$i -o Z.$i
wait

done

############################################
####make matrix for trait with repeatability
############################################
#
#array=(mcw)
#
#for i in "${array[@]}"
#do
#
######make X and Z matrix
#awk '{print $2}' phe.$i".fill.ped" > y.$i &
#awk '{print $1}' phe.$i".fill.ped" > id.dat.$i &
#awk '{print $1}' stacked_ped > id.eff &
#awk '{print $4}' phe.$i".fill.ped" > cg1.dat.$i &
#awk '{print $4}' phe.$i".fill.ped" | sort -u > cg1.eff.$i &
#wait
#cgen_z -d cg1.dat.$i -e cg1.eff.$i -r y.$i -o X.$i
#cgen_z -d id.dat.$i -e id.eff -r y.$i -o Z.$i
#wait


###########################################
## rename design matrix
###########################################
n=1
array=(bw ww yw fw carc_ema carc_imf carc_rib cwt ema imf rib mcw ss)

for i in "${array[@]}"
do

cat y.$i >y$n &
cat Z.$i >Z$n &
cat X.$i >X$n &
wait

((n+=1))

done

n=1
array=(bw ww ss)

for i in "${array[@]}"
do
cat M.$i > M$n & 
cat P.$i > P$n &
wait

((n+=1))

done

###########################################
###########################################

cnewr -R R -r resid.14x14 y1 y2 y3 y4 y5 y6 y7 y9 y9 y10 y11 y12 y13 y14
wait

###########################################
##construct X'RX; Z'RZ; X'RZ'; Z'RX 
###########################################
array=(1 2 3 4 5 6 7 8 9 10 11 12)
array2=(1 2 3 4 5 6 7 8 9 10 11 12)

for i in "${array[@]}"
do
for j in "${array2[@]}"
do

if[ $j < $i ]{
cat R$j$i > R$i$j
}

#construct X'RX
cmult -t -a X$i -R R$i$j -b X$j -c X$i"RX"$j &

#construct Z'RZ
cmult -t -a Z$i -R R$i$j -b Z$j -c Z$i"RZ"$j &

#construct X'RZ
cmult -t -a X$i -R R$i$j -b Z$j -c X$i"RZ"$j &

#construct Z'RX
cmult -t -a Z$i -R R$i$j -b X$j -c Z$i"RX"$j &

done
done
###########################################
##construct Zm'RZm ; Zp'RZp ; Zm'RZp ; Zp'RZm 
###########################################
array=(1 2 3)
array2=(1 2 3)

for i in "${array[@]}"
do
for j in "${array2[@]}"
do

if[ $j < $i ]{
cat R$j$i > R$i$j
}
#construct Zm'RZm
cmult -t -a M$i -R R$i$j -b M$j -c M$i"RM"$j &

#construct Zp'RZp
cmult -t -a P$i -R R$i$j -b P$j -c P$i"RP"$j &

#construct Zm'RZp
cmult -t -a M$i -R R$i$j -b P$j -c M$i"RP"$j &

#construct Zp'RZm
cmult -t -a P$i -R R$i$j -b M$j -c P$i"RM"$j &

done
done

###########################################
##construct X'RZm ; X'RZp; Z'RZm ; Z'RZp  X 2 
###########################################
array=(1 2 3 4 5 6 7 8 9 10 11 12)
array2=(1 2 3) 

for i in "${array[@]}"
do
for j in "${array2[@]}"
do


if[ $j < $i ]{
cat R$j$i > R$i$j
}

#construct X'RZm
cmult -t -a X$i -R R$i$j -b M$j -c X$i"RM"$j &

#construct Zm'RX
cmult -t -a M$j -R R$i$j -b X$i -c M$j"RX"$i &

#construct Z'RZm
cmult -t -a Z$i -R R$i$j -b M$j -c Z$i"RM"$j &

#construct Zm'RZ
cmult -t -a M$j -R R$i$j -b Z$i -c M$j"RZ"$i &

#construct X'RZp
cmult -t -a X$i -R R$i$j -b P$j -c X$i"RP"$j &

#construct Zp'RX
cmult -t -a P$j -R R$i$j -b X$i -c P$j"RX"$i &

#construct Z'RZp
cmult -t -a Z$i -R R$i$j -b P$j -c Z$i"RP"$j &
#construct Zp'RZ
cmult -t -a P$j -R R$i$j -b Z$i -c P$j"RZ"$i &

done
done

wait

###########################################
##construct 
###########################################

array=(1 2 3 4 5 6 7 8 9 10 11 12)
array2=(1 2 3 4 5 6 7 8 9 10 11 12) 

for i in "${array[@]}"
do
for j in "${array2[@]}"
do

gInvert =$(awk '$1==$i&&$2=$j{print $3}' Ginvert)

cadd -a Z$i"RZ"$j -r gInvert Ainverse -c Z$i"RZ"$j".v" &

done
done

wait


array=(1 2 3)
array2=(1 2 3) 

for i in "${array[@]}"
do
for j in "${array2[@]}"
do

mInvert =$(awk '$1==$i&&$2=$j{print $3}' Minvert)

cadd -a M$i"RM"$j -r mInvert -b Ainverse -c M$i"RM"$j".v" &
done
done
wait


indnum=$(awk '$6!="\."{print $6}' phe.bw.fill.ped|sort -u|wc| awk '{print $1}')  

array=(1 2 3)
array2=(1 2 3) 

for i in "${array[@]}"
do
for j in "${array2[@]}"
do

ident 1611519 > myIdent
cadd -a P$i"RP"$j -r $p11 -b myIndent -c P$i"RP"$j".v" &

done
done

wait

###########################################
##rhs 
###########################################

array=(1 2 3 4 5 6 7 8 9 10 11 12)
array2=(1 2 3 4 5 6 7 8 9 10 11 12)

for i in "${array[@]}"
do
for j in "${array2[@]}"
do

cmult -t -a X$i -R R$i$j -b y$j -c X$i"Ry"$j &
cmult -t -a Z$i -R R$i$j -b y$j -c Z$i"Ry"$j &

done
done

array=(1 2)
array2=(1 2 3 4)

for i in "${array[@]}"
do
for j in "${array2[@]}"
do

cmult -t -a M$i -R R$i$j -b y$j -c M$i"Ry"$j &
cmult -t -a P$i -R R$i$j -b y$j -c P$i"Ry"$j &

done
done

wait

###########################################
##concat rhs 
###########################################

array=(X1 X2 X3 X4 Z1 Z2 Z3 Z4 M1 M2 P1 P2)

for i in "${array[@]}" #X,Z..
do

cadd -a $i"Ry1" -b $i"Ry2" -c rhs.1.1 &
cadd -a $i"Ry3" -b $i"Ry4" -c rhs.1.2 &
cadd -a $i"Ry5" -b $i"Ry6" -c rhs.1.3 &
cadd -a $i"Ry7" -b $i"Ry8" -c rhs.1.4 &
cadd -a $i"Ry9" -b $i"Ry10" -c rhs.1.5 &
cadd -a $i"Ry11" -b $i"Ry12" -c rhs.1.6 &
wait
cadd -a rhs.1.1 -b rhs.1.2 -c rhs.1.7 &
cadd -a rhs.1.3 -b rhs.1.4 -c rhs.1.8 &
cadd -a rhs.1.5 -b rhs.1.6 -c rhs.1.9 &
wait
cadd -a rhs.1.7 -b rhs.1.8 -c rhs.1.10 &
cadd -a rhs.1.9 -b $i"Ry13" -c rhs.1.11 &
wait
cadd -a rhs.1.10 -b rhs.1.11 -c rhs.$i

done

#pcgmgpu -A lhs -x solve_gpu.sln -b rhs -t 1e-15 > solve_gpu.log
