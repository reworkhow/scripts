#!/bin/bash


##########################################
###make matrix for trait with maternal effects
##########################################
array=(bw ww)

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
##########################################
###make matrix for trait with maternal effect but just genetic part
##########################################
array=(ss)

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

done

###########################################
###make matrix for trait without maternal effects
###########################################

array=(fw yw carc_ema carc_imf carc_rib cwt ema imf rib mcw)

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
array=(bw ww yw fw mcw imf ema rib ss cwt carc_ema carc_rib carc_imf)

for i in "${array[@]}"
do

cat y.$i >y$n &
cat Z.$i >Z$n &
cat X.$i >X$n &
wait

((n+=1))

done

n=1
array=(bw ww)

for i in "${array[@]}"
do
cat M.$i > M$n & 
cat P.$i > P$n &
wait

((n+=1))

done


cat M.ss > M3

wait 

###########################################
###########################################

cnewr -R R -r resid.13x13 y1 y2 y3 y4 y5 y6 y7 y9 y9 y10 y11 y12 y13
wait

###########################################
##construct X'RX; Z'RZ; X'RZ'; Z'RX 
###########################################
array=(1 2 3 4 5 6 7 8 9 10 11 12 13)
array2=(1 2 3 4 5 6 7 8 9 10 11 12 13) 

for i in "${array[@]}"
do
for j in "${array2[@]}"
do

if(($j < $i ));then
{
cat R$j$i > R$i$j
}
fi
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
##construct Zm'RZm ; Zp'RZp ; Zm'RZp ; Zp'RZm for traits with both m and p 
###########################################
array=(1 2)
array2=(1 2)

for i in "${array[@]}"
do
for j in "${array2[@]}"
do

if(($j < $i ));then
{
cat R$j$i > R$i$j
}
fi
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
##construct Zm'RZm for trait with only m
###########################################
array=(3)
array2=(3)

for i in "${array[@]}"
do
for j in "${array2[@]}"
do

if(($j < $i ));then
{
cat R$j$i > R$i$j
}
fi
#construct Zm'RZm
cmult -t -a M$i -R R$i$j -b M$j -c M$i"RM"$j &

done
done
###########################################
##construct X'RZm ; X'RZp; Z'RZm ; Z'RZp  X 2 for trait with m and p 
###########################################
array=(1 2 3 4 5 6 7 8 9 10 11 12 13)
array2=(1 2) 

for i in "${array[@]}"
do
for j in "${array2[@]}"
do

if(($j < $i ));then
{
cat R$j$i > R$i$j
}
fi

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
##construct X'RZm and Z'RZm for trit with only m 
###########################################
array=(1 2 3 4 5 6 7 8 9 10 11 12 13)
array2=(3) 

for i in "${array[@]}"
do
for j in "${array2[@]}"
do


if(($j < $i ));then
{
cat R$j$i > R$i$j
}
fi

#construct X'RZm
cmult -t -a X$i -R R$i$j -b M$j -c X$i"RM"$j &

#construct Zm'RX
cmult -t -a M$j -R R$i$j -b X$i -c M$j"RX"$i &

#construct Z'RZm
cmult -t -a Z$i -R R$i$j -b M$j -c Z$i"RM"$j &

#construct Zm'RZ
cmult -t -a M$j -R R$i$j -b Z$i -c M$j"RZ"$i &


done
done

wait
############################################
###construct 
############################################
#
#array=(1 2 3 4 5 6 7 8 9 10 11 12 13)
#array2=(1 2 3 4 5 6 7 8 9 10 11 12 13) 
#
#for i in "${array[@]}"
#do
#for j in "${array2[@]}"
#do
#
#gInvert =$(awk '$1==$i&&$2=$j{print $3}' genetic.13x13.inverted)
#
#cadd -a Z$i"RZ"$j -r $gInvert -b Ainverse -c Z$i"RZ"$j".v" &
#
#done
#done
#
#wait
#
#
#array=(1 2 3)
#array2=(1 2 3) 
#
#for i in "${array[@]}"
#do
#for j in "${array2[@]}"
#do
#
#mInvert =$(awk '$1==$i&&$2=$j{print $3}' maternal.3x3.inverted)
#
#cadd -a M$i"RM"$j -r $mInvert -b Ainverse -c M$i"RM"$j".v" &
#done
#done
#wait
#
#
#indnum=$(awk '$6!="\."{print $6}' phe.bw.fill.ped|sort -u|wc| awk '{print $1}')  
#
#array=(1 2)
#array2=(1 2) 
#
#for i in "${array[@]}"
#do
#for j in "${array2[@]}"
#do
#
#pInvert =$(awk '$1==$i&&$2=$j{print $3}' permanent.2x2.inverted)
#
#ident $indnum  > myIdent
#cadd -a P$i"RP"$j -r $pInvert -b $myIndent -c P$i"RP"$j".v" &
#
#done
#done
#
#wait
#
############################################
###rhs 
############################################
#
#array=(1 2 3 4 5 6 7 8 9 10 11 12 13)
#array2=(1 2 3 4 5 6 7 8 9 10 11 12 13)
#
#for i in "${array[@]}"
#do
#for j in "${array2[@]}"
#do
#
#cmult -t -a X$i -R R$i$j -b y$j -c X$i"Ry"$j &
#cmult -t -a Z$i -R R$i$j -b y$j -c Z$i"Ry"$j &
#
#done
#done
#
#array=(1 2)
#array2=(1 2 3 4 5 6 7 8 9 10 11 12 13)
#
#for i in "${array[@]}"
#do
#for j in "${array2[@]}"
#do
#
#cmult -t -a M$i -R R$i$j -b y$j -c M$i"Ry"$j &
#cmult -t -a P$i -R R$i$j -b y$j -c P$i"Ry"$j &
#
#done
#done
#
#array=(3)
#array2=(1 2 3 4 5 6 7 8 9 10 11 12 13)
#
#for i in "${array[@]}"
#do
#for j in "${array2[@]}"
#do
#
#cmult -t -a M$i -R R$i$j -b y$j -c M$i"Ry"$j &
#
#done
#done
#
#wait
#
############################################
###concat rhs 
############################################
#
#array=(X1 X2 X3 X4 X5 X6 X7 X8 X9 X10 X11 X12 X13 Z1 Z2 Z3 Z4 Z5 Z6 Z7 Z8 Z9 Z10 Z11 Z12 Z13  M1 M2 M3 P1 P2)
#
#for i in "${array[@]}" #X,Z..
#do
#
#cadd -a $i"Ry1" -b $i"Ry2" -c rhs.1.1 &
#cadd -a $i"Ry3" -b $i"Ry4" -c rhs.1.2 &
#cadd -a $i"Ry5" -b $i"Ry6" -c rhs.1.3 &
#cadd -a $i"Ry7" -b $i"Ry8" -c rhs.1.4 &
#cadd -a $i"Ry9" -b $i"Ry10" -c rhs.1.5 &
#cadd -a $i"Ry11" -b $i"Ry12" -c rhs.1.6 &
#wait
#cadd -a rhs.1.1 -b rhs.1.2 -c rhs.1.7 &
#cadd -a rhs.1.3 -b rhs.1.4 -c rhs.1.8 &
#cadd -a rhs.1.5 -b rhs.1.6 -c rhs.1.9 &
#wait
#cadd -a rhs.1.7 -b rhs.1.8 -c rhs.1.10 &
#cadd -a rhs.1.9 -b $i"Ry13" -c rhs.1.11 &
#wait
#cadd -a rhs.1.10 -b rhs.1.11 -c rhs.$i
#
#done
#
#wait
#
#cvcat rhs.X1 rhs.X2 rhs.X3 rhs.X4 rhs.X5 rhs.X6 rhs.X7 rhs.X8 rhs.X9 rhs.X10 rhs.X11 rhs.X12 rhs.X13 rhs.Z1 rhs.Z2 rhs.Z3 rhs.Z4 rhs.Z5 rhs.Z6 rhs.Z7 rhs.Z8 rhs.Z9 rhs.Z10 rhs.Z11 rhs.Z12 rhs.Z13 rhs.M1 rhs.M2 rhs.M3 rhs.P1 rhs.p2 rhs
#
#
#pcgmgpu -A lhs.MAP -x solve_gpu.sln -b rhs -t 1e-15 > solve_gpu.log
