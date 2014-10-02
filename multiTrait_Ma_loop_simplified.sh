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

###########################################
###make matrix for trait without maternal effects
###########################################

array=(fw yw)

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


###########################################
## rename design matrix
###########################################
n=1
array=(bw ww yw fw)

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

###########################################
###########################################

cnewr -R R -r resid.4x4 y1 y2 y3 y4
wait

###########################################
##construct X'RX; Z'RZ; X'RZ'; Z'RX 
###########################################
array=(X1 X2 X3 X4)
array=(X1 X2 X3 X4)

cat R12 > R21

for i in "${array[@]}"
do
for j in "${array2[@]}"
do

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
array=(1 2)
array=(1 2)

for i in "${array[@]}"
do
for j in "${array2[@]}"
do
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
array=(1 2 3 4)
array=(1 2)

for i in "${array[@]}"
do
for j in "${array2[@]}"
do

#construct X'RZm
cmult -t -a X$i -R R$i$j -b M$j -c X$i"RM"$j &

#construct Zm'RX
cmult -t -a M$j -R R$j$i -b X$i -c M$j"RX"$i &

#construct Z'RZm
cmult -t -a Z$i -R R$i$j -b M$j -c Z$i"RM"$j &

#construct Zm'RZ
cmult -t -a M$j -R R$j$i -b Z$i -c M$j"RZ"$i &

#construct X'RZp
cmult -t -a X$i -R R$i$j -b P$j -c X$i"RP"$j &

#construct Zp'RX
cmult -t -a P$j -R R$j$i -b X$i -c P$j"RX"$i &

#construct Z'RZp
cmult -t -a Z$i -R R$i$j -b P$j -c Z$i"RP"$j &
#construct Zp'RZ
cmult -t -a P$j -R R$j$i -b Z$i -c P$j"RZ"$i &

done
done

wait

###########################################
##construct 
###########################################

g11=0.0505
g12=-0.0033
g13=-0.0031
g14=0.0015
g21=-0.0033
g22=0.0044
g23=-0.0015
g24=-0.0005
g31=-0.0031
g32=-0.0015
g33=0.0023
g34=-0.0011
g41=0.0015
g42=-0.0005
g43=-0.0011
g44=0.0015



cadd -a Z$i"RZ"$j -r $g11 Ainverse -c Z$i"RZ"$j".v" &

cadd -a Z1RZ1 -r $g11 -b Ainverse -c Z1RZ1.v &
cadd -a Z1RZ2 -r $g12 -b Ainverse -c Z1RZ2.v &
cadd -a Z1RZ3 -r $g13 -b Ainverse -c Z1RZ3.v &
cadd -a Z1RZ4 -r $g14 -b Ainverse -c Z1RZ4.v &
cadd -a Z2RZ1 -r $g21 -b Ainverse -c Z2RZ1.v &
cadd -a Z2RZ2 -r $g22 -b Ainverse -c Z2RZ2.v &
cadd -a Z2RZ3 -r $g23 -b Ainverse -c Z2RZ3.v &
cadd -a Z2RZ4 -r $g24 -b Ainverse -c Z2RZ4.v &
cadd -a Z3RZ1 -r $g31 -b Ainverse -c Z3RZ1.v &
cadd -a Z3RZ2 -r $g32 -b Ainverse -c Z3RZ2.v &
cadd -a Z3RZ3 -r $g33 -b Ainverse -c Z3RZ3.v &
cadd -a Z3RZ4 -r $g34 -b Ainverse -c Z3RZ4.v &
cadd -a Z4RZ1 -r $g41 -b Ainverse -c Z4RZ1.v &
cadd -a Z4RZ2 -r $g42 -b Ainverse -c Z4RZ2.v &
cadd -a Z4RZ3 -r $g43 -b Ainverse -c Z4RZ3.v &
cadd -a Z4RZ4 -r $g44 -b Ainverse -c Z4RZ4.v &
wait

m11=0.2009
m12=-0.0069
m21=-0.0069
m22=0.0031
cadd -a M1RM1 -r $m11 -b Ainverse -c M1RM1.v &
cadd -a M1RM2 -r $m12 -b Ainverse -c M1RM2.v &
cadd -a M2RM1 -r $m21 -b Ainverse -c M2RM1.v &
cadd -a M2RM2 -r $m22 -b Ainverse -c M2RM2.v &
wait


p11=0.4994
p12=-0.0092
p21=-0.0092
p22=0.0019
indnum=$(awk '$6!="\."{print $6}' phe.bw.fill.ped|sort -u|wc| awk '{print $1}')  
ident 1611519 > ident.1611519
cadd -a P1RP1 -r $p11 -b ident.1611519 -c P1RP1.v &
cadd -a P1RP2 -r $p12 -b ident.1611519 -c P1RP2.v &
cadd -a P2RP1 -r $p21 -b ident.1611519 -c P2RP1.v &
cadd -a P2RP2 -r $p22 -b ident.1611519 -c P2RP2.v &
wait

###########################################
##rhs 
###########################################

array=(1 2 3 4)
array2=(1 2 3 4)

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


for i in "${array[@]}" #X,Z..
do


cvcat rhs.X1 rhs.X2 rhs.X3 rhs.X4 rhs.Z1 rhs.Z2 rhs.Z3 rhs.Z4 rhs.M1 rhs.M2 rhs.P1 rhs.P2 rhs

#csolve -A lhs.MAP -x solve.Ma,sln -b rhs -h -t 1e-15 > solve.Ma.log
#pcgmgpu -A lhs -x solve_gpu.sln -b rhs -t 1e-15 > solve_gpu.log
