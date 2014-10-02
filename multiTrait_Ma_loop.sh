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
##construct X'RX 
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
cmult -t -a X$i -R R11 -b Z$j -c X$i"RZ"$j &

#construct Z'RX
cmult -t -a Z$i -R R11 -b X$j -c Z$i"RX"$j &




#construct Zm'RZm
cmult -t -a M1 -R R11 -b M1 -c M1RM1 &
cmult -t -a M1 -R R12 -b M2 -c M1RM2 &
cmult -t -a M2 -R R12 -b M1 -c M2RM1 &
cmult -t -a M2 -R R22 -b M2 -c M2RM2 &

#construct Zp'RZp
cmult -t -a P1 -R R11 -b P1 -c P1RP1 &
cmult -t -a P2 -R R12 -b P1 -c P2RP1 &
cmult -t -a P1 -R R12 -b P2 -c P1RP2 &
cmult -t -a P2 -R R22 -b P2 -c P2RP2 &

#construct X'RZm
cmult -t -a X1 -R R11 -b M1 -c X1RM1 &
cmult -t -a X1 -R R12 -b M2 -c X1RM2 &
cmult -t -a X2 -R R12 -b M1 -c X2RM1 &
cmult -t -a X2 -R R22 -b M2 -c X2RM2 &
cmult -t -a X3 -R R13 -b M1 -c X3RM1 &
cmult -t -a X3 -R R23 -b M2 -c X3RM2 &
cmult -t -a X4 -R R14 -b M1 -c X4RM1 &
cmult -t -a X4 -R R24 -b M2 -c X4RM2 &
#construct Zm'RX
cmult -t -a M1 -R R11 -b X1 -c M1RX1 &
cmult -t -a M2 -R R12 -b X1 -c M2RX1 &
cmult -t -a M1 -R R12 -b X2 -c M1RX2 &
cmult -t -a M2 -R R22 -b X2 -c M2RX2 &
cmult -t -a M1 -R R13 -b X3 -c M1RX3 &
cmult -t -a M2 -R R23 -b X3 -c M2RX3 &
cmult -t -a M1 -R R14 -b X4 -c M1RX4 &
cmult -t -a M2 -R R24 -b X4 -c M2RX4 &

#construct Z'RZm
cmult -t -a Z1 -R R11 -b M1 -c Z1RM1 &
cmult -t -a Z1 -R R12 -b M2 -c Z1RM2 &
cmult -t -a Z2 -R R12 -b M1 -c Z2RM1 &
cmult -t -a Z2 -R R22 -b M2 -c Z2RM2 &
cmult -t -a Z3 -R R13 -b M1 -c Z3RM1 &
cmult -t -a Z3 -R R23 -b M2 -c Z3RM2 &
cmult -t -a Z4 -R R14 -b M1 -c Z4RM1 &
cmult -t -a Z4 -R R24 -b M2 -c Z4RM2 &
#construct Zm'RZ
cmult -t -a M1 -R R11 -b Z1 -c M1RZ1 &
cmult -t -a M2 -R R12 -b Z1 -c M2RZ1 &
cmult -t -a M1 -R R12 -b Z2 -c M1RZ2 &
cmult -t -a M2 -R R22 -b Z2 -c M2RZ2 &
cmult -t -a M1 -R R13 -b Z3 -c M1RZ3 &
cmult -t -a M2 -R R23 -b Z3 -c M2RZ3 &
cmult -t -a M1 -R R14 -b Z4 -c M1RZ4 &
cmult -t -a M2 -R R24 -b Z4 -c M2RZ4 &

#construct X'RZp
cmult -t -a X1 -R R11 -b P1 -c X1RP1 &
cmult -t -a X1 -R R12 -b P2 -c X1RP2 &
cmult -t -a X2 -R R12 -b P1 -c X2RP1 &
cmult -t -a X2 -R R22 -b P2 -c X2RP2 &
cmult -t -a X3 -R R13 -b P1 -c X3RP1 &
cmult -t -a X3 -R R23 -b P2 -c X3RP2 &
cmult -t -a X4 -R R14 -b P1 -c X4RP1 &
cmult -t -a X4 -R R24 -b P2 -c X4RP2 &
#construct Zp'RX
cmult -t -a P1 -R R11 -b X1 -c P1RX1 &
cmult -t -a P2 -R R12 -b X1 -c P2RX1 &
cmult -t -a P1 -R R12 -b X2 -c P1RX2 &
cmult -t -a P2 -R R22 -b X2 -c P2RX2 &
cmult -t -a P1 -R R13 -b X3 -c P1RX3 &
cmult -t -a P2 -R R23 -b X3 -c P2RX3 &
cmult -t -a P1 -R R14 -b X4 -c P1RX4 &
cmult -t -a P2 -R R24 -b X4 -c P2RX4 &


#construct Z'RZp
cmult -t -a Z1 -R R11 -b P1 -c Z1RP1 &
cmult -t -a Z1 -R R12 -b P2 -c Z1RP2 &
cmult -t -a Z2 -R R12 -b P1 -c Z2RP1 &
cmult -t -a Z2 -R R22 -b P2 -c Z2RP2 &
cmult -t -a Z3 -R R13 -b P1 -c Z3RP1 &
cmult -t -a Z3 -R R23 -b P2 -c Z3RP2 &
cmult -t -a Z4 -R R14 -b P1 -c Z4RP1 &
cmult -t -a Z4 -R R24 -b P2 -c Z4RP2 &
#construct Zp'RZ
cmult -t -a P1 -R R11 -b Z1 -c P1RZ1 &
cmult -t -a P2 -R R12 -b Z1 -c P2RZ1 &
cmult -t -a P1 -R R12 -b Z2 -c P1RZ2 &
cmult -t -a P2 -R R22 -b Z2 -c P2RZ2 &
cmult -t -a P1 -R R13 -b Z3 -c P1RZ3 &
cmult -t -a P2 -R R23 -b Z3 -c P2RZ3 &
cmult -t -a P1 -R R14 -b Z4 -c P1RZ4 &
cmult -t -a P2 -R R24 -b Z4 -c P2RZ4 &


#construct Zm'RZp
cmult -t -a M1 -R R11 -b P1 -c M1RP1 &
cmult -t -a M1 -R R12 -b P2 -c M1RP2 &
cmult -t -a M2 -R R12 -b P1 -c M2RP1 &
cmult -t -a M2 -R R22 -b P2 -c M2RP2 &
#construct Zp'RZm
cmult -t -a P1 -R R11 -b M1 -c P1RM1 &
cmult -t -a P2 -R R12 -b M1 -c P2RM1 &
cmult -t -a P1 -R R12 -b M2 -c P1RM2 &
cmult -t -a P2 -R R22 -b M2 -c P2RM2 &

wait

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


cmult -t -a X1 -R R11 -b y1 -c X1Ry1 &
cmult -t -a X1 -R R12 -b y2 -c X1Ry2 &
cmult -t -a X1 -R R13 -b y3 -c X1Ry3 &
cmult -t -a X1 -R R14 -b y4 -c X1Ry4 &
cmult -t -a X2 -R R12 -b y1 -c X2Ry1 &
cmult -t -a X2 -R R22 -b y2 -c X2Ry2 &
cmult -t -a X2 -R R23 -b y3 -c X2Ry3 &
cmult -t -a X2 -R R24 -b y4 -c X2Ry4 &
cmult -t -a X3 -R R13 -b y1 -c X3Ry1 &
cmult -t -a X3 -R R23 -b y2 -c X3Ry2 &
cmult -t -a X3 -R R33 -b y3 -c X3Ry3 &
cmult -t -a X3 -R R34 -b y4 -c X3Ry4 &
cmult -t -a X4 -R R14 -b y1 -c X4Ry1 &
cmult -t -a X4 -R R24 -b y2 -c X4Ry2 &
cmult -t -a X4 -R R34 -b y3 -c X4Ry3 &
cmult -t -a X4 -R R44 -b y4 -c X4Ry4 &

cmult -t -a Z1 -R R11 -b y1 -c Z1Ry1 &
cmult -t -a Z1 -R R12 -b y2 -c Z1Ry2 &
cmult -t -a Z1 -R R13 -b y3 -c Z1Ry3 &
cmult -t -a Z1 -R R14 -b y4 -c Z1Ry4 &
cmult -t -a Z2 -R R12 -b y1 -c Z2Ry1 &
cmult -t -a Z2 -R R22 -b y2 -c Z2Ry2 &
cmult -t -a Z2 -R R23 -b y3 -c Z2Ry3 &
cmult -t -a Z2 -R R24 -b y4 -c Z2Ry4 &
cmult -t -a Z3 -R R13 -b y1 -c Z3Ry1 &
cmult -t -a Z3 -R R23 -b y2 -c Z3Ry2 &
cmult -t -a Z3 -R R33 -b y3 -c Z3Ry3 &
cmult -t -a Z3 -R R34 -b y4 -c Z3Ry4 &
cmult -t -a Z4 -R R14 -b y1 -c Z4Ry1 &
cmult -t -a Z4 -R R24 -b y2 -c Z4Ry2 &
cmult -t -a Z4 -R R34 -b y3 -c Z4Ry3 &
cmult -t -a Z4 -R R44 -b y4 -c Z4Ry4 &

cmult -t -a M1 -R R11 -b y1 -c M1Ry1 &
cmult -t -a M1 -R R12 -b y2 -c M1Ry2 &
cmult -t -a M1 -R R13 -b y3 -c M1Ry3 &
cmult -t -a M1 -R R14 -b y4 -c M1Ry4 &
cmult -t -a M2 -R R12 -b y1 -c M2Ry1 &
cmult -t -a M2 -R R22 -b y2 -c M2Ry2 &
cmult -t -a M2 -R R23 -b y3 -c M2Ry3 &
cmult -t -a M2 -R R24 -b y4 -c M2Ry4 &

cmult -t -a P1 -R R11 -b y1 -c P1Ry1 &
cmult -t -a P1 -R R12 -b y2 -c P1Ry2 &
cmult -t -a P1 -R R13 -b y3 -c P1Ry3 &
cmult -t -a P1 -R R14 -b y4 -c P1Ry4 &
cmult -t -a P2 -R R12 -b y1 -c P2Ry1 &
cmult -t -a P2 -R R22 -b y2 -c P2Ry2 &
cmult -t -a P2 -R R23 -b y3 -c P2Ry3 &
cmult -t -a P2 -R R24 -b y4 -c P2Ry4 &

wait

cadd -a X1Ry1 -b X1Ry2 -c rhs.1.1 &
cadd -a X1Ry3 -b X1Ry4 -c rhs.1.2 &
wait
cadd -a rhs.1.1 -b rhs.1.2 -c rhs.1

cadd -a X2Ry1 -b X2Ry2 -c rhs.2.1 &
cadd -a X2Ry3 -b X2Ry4 -c rhs.2.2 &
wait
cadd -a rhs.2.1 -b rhs.2.2 -c rhs.2

cadd -a X3Ry1 -b X3Ry2 -c rhs.3.1 &
cadd -a X3Ry3 -b X3Ry4 -c rhs.3.2 &
wait
cadd -a rhs.3.1 -b rhs.3.2 -c rhs.3

cadd -a X4Ry1 -b X4Ry2 -c rhs.4.1 &
cadd -a X4Ry3 -b X4Ry4 -c rhs.4.2 &
wait
cadd -a rhs.4.1 -b rhs.4.2 -c rhs.4

cadd -a Z1Ry1 -b Z1Ry2 -c rhs.5.1 &
cadd -a Z1Ry3 -b Z1Ry4 -c rhs.5.2 &
wait
cadd -a rhs.5.1 -b rhs.5.2 -c rhs.5

cadd -a Z2Ry1 -b Z2Ry2 -c rhs.6.1 &
cadd -a Z2Ry3 -b Z2Ry4 -c rhs.6.2 &
wait
cadd -a rhs.6.1 -b rhs.6.2 -c rhs.6

cadd -a Z3Ry1 -b Z3Ry2 -c rhs.7.1 &
cadd -a Z3Ry3 -b Z3Ry4 -c rhs.7.2 &
wait
cadd -a rhs.7.1 -b rhs.7.2 -c rhs.7

cadd -a Z4Ry1 -b Z4Ry2 -c rhs.8.1 &
cadd -a Z4Ry3 -b Z4Ry4 -c rhs.8.2 &
wait
cadd -a rhs.8.1 -b rhs.8.2 -c rhs.8

cadd -a M1Ry1 -b M1Ry2 -c rhs.9.1 &
cadd -a M1Ry3 -b M1Ry4 -c rhs.9.2 &
wait
cadd -a rhs.9.1 -b rhs.9.2 -c rhs.9

cadd -a M2Ry1 -b M2Ry2 -c rhs.10.1 &
cadd -a M2Ry3 -b M2Ry4 -c rhs.10.2 &
wait
cadd -a rhs.10.1 -b rhs.10.2 -c rhs.10

cadd -a P1Ry1 -b P1Ry2 -c rhs.11.1 &
cadd -a P1Ry3 -b P1Ry4 -c rhs.11.2 &
wait
cadd -a rhs.11.1 -b rhs.11.2 -c rhs.11

cadd -a P2Ry1 -b P2Ry2 -c rhs.12.1 &
cadd -a P2Ry3 -b P2Ry4 -c rhs.12.2 &
wait
cadd -a rhs.12.1 -b rhs.12.2 -c rhs.12

cvcat rhs.1 rhs.2 rhs.3 rhs.4 rhs.5 rhs.6 rhs.7 rhs.8 rhs.9 rhs.10 rhs.11 rhs.12 rhs

#csolve -A lhs.MAP -x solve.Ma,sln -b rhs -h -t 1e-15 > solve.Ma.log
#pcgmgpu -A lhs -x solve_gpu.sln -b rhs -t 1e-15 > solve_gpu.log
