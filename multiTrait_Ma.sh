#!/bin/bash

#####birth weight
awk '{print $2}' phe.bw.fill.ped > y.bw &
awk '{print $1}' phe.bw.fill.ped > id.dat.bw &
awk '{print $1}' stacked_ped > id.eff &
awk '{print $4}' phe.bw.fill.ped > cg1.dat.bw &
awk '{print $4}' phe.bw.fill.ped | sort -u > cg1.eff.bw &
wait
cgen_z -d cg1.dat.bw -e cg1.eff.bw -r y.bw -o X.bw
cgen_z -d id.dat.bw -e id.eff -r y.bw -o Z.bw
wait

#maternal effect
awk '{print $1}' stacked_ped > mat.eff.bw &
awk '{print $6}' phe.bw.fill.ped > mat.dat.bw &
wait
cgen_z -d mat.dat.bw -e mat.eff.bw -r y.bw -o M.bw

#permanent env effect
awk '$6!="\."{print $6}' phe.bw.fill.ped | sort -u > perm.eff.bw &
awk '{print $6}' phe.bw.fill.ped > perm.dat.bw &
wait
cgen_z -d perm.dat.bw -e perm.eff.bw -r y.bw -o P.bw


#####weaning weight
awk '{print $2}' phe.ww.fill.ped > y.ww &
awk '{print $1}' phe.ww.fill.ped > id.dat.ww &
awk '{print $1}' stacked_ped > id.eff &
awk '{print $4}' phe.ww.fill.ped > cg1.dat.ww &
awk '{print $4}' phe.ww.fill.ped | sort -u > cg1.eff.ww &
wait
cgen_z -d cg1.dat.ww -e cg1.eff.ww -r y.ww -o X.ww
cgen_z -d id.dat.ww -e id.eff -r y.ww -o Z.ww
wait

#maternal effect
awk '{print $1}' stacked_ped > mat.eff.ww &
awk '{print $6}' phe.ww.fill.ped > mat.dat.ww &
wait
cgen_z -d mat.dat.ww -e mat.eff.ww -r y.ww -o M.ww

#permanent env effect
awk '$6!="\."{print $6}' phe.ww.fill.ped | sort -u > perm.eff.ww &
awk '{print $6}' phe.ww.fill.ped > perm.dat.ww &
wait
cgen_z -d perm.dat.ww -e perm.eff.ww -r y.bw -o P.ww


#####final weight
awk '{print $2}' phe.fw.fill.ped > y.fw &
awk '{print $1}' phe.fw.fill.ped > id.dat.fw &
awk '{print $1}' stacked_ped > id.eff &
awk '{print $4}' phe.fw.fill.ped > cg1.dat.fw &
awk '{print $4}' phe.fw.fill.ped | sort -u > cg1.eff.fw &
wait
cgen_z -d cg1.dat.fw -e cg1.eff.fw -r y.fw -o X.fw
cgen_z -d id.dat.fw -e id.eff -r y.fw -o Z.fw
wait

#####yearling weight
awk '{print $2}' phe.yw.fill.ped > y.yw &
awk '{print $1}' phe.yw.fill.ped > id.dat.yw &
awk '{print $1}' stacked_ped > id.eff &
awk '{print $4}' phe.yw.fill.ped > cg1.dat.yw &
awk '{print $4}' phe.yw.fill.ped | sort -u > cg1.eff.yw &
wait
cgen_z -d cg1.dat.yw -e cg1.eff.yw -r y.yw -o X.yw
cgen_z -d id.dat.yw -e id.eff -r y.yw -o Z.yw
wait

cat y.bw >y1 &
cat y.ww >y2 &
cat y.yw >y3 &
cat y.fw >y4 &
cat Z.bw >Z1 &
cat Z.ww >Z2 &
cat Z.yw >Z3 &
cat Z.fw >Z4 &
cat X.bw >X1 &
cat X.ww >X2 &
cat X.yw >X3 &
cat X.fw >X4 &

cat M.bw > M1 
cat P.bw > P1
cat M.ww > M2
cat P.ww > P2
wait

cnewr -R R -r resid.4x4 y1 y2 y3 y4
wait
#construct X'RX
cmult -t -a X1 -R R11 -b X1 -c X1RX1 &
cmult -t -a X1 -R R12 -b X2 -c X1RX2 &
cmult -t -a X1 -R R13 -b X3 -c X1RX3 &
cmult -t -a X1 -R R14 -b X4 -c X1RX4 &
cmult -t -a X2 -R R12 -b X1 -c X2RX1 &
cmult -t -a X2 -R R22 -b X2 -c X2RX2 &
cmult -t -a X2 -R R23 -b X3 -c X2RX3 &
cmult -t -a X2 -R R24 -b X4 -c X2RX4 &
cmult -t -a X3 -R R13 -b X1 -c X3RX1 &
cmult -t -a X3 -R R23 -b X2 -c X3RX2 &
cmult -t -a X3 -R R33 -b X3 -c X3RX3 &
cmult -t -a X3 -R R34 -b X4 -c X3RX4 &
cmult -t -a X4 -R R14 -b X1 -c X4RX1 &
cmult -t -a X4 -R R24 -b X2 -c X4RX2 &
cmult -t -a X4 -R R34 -b X3 -c X4RX3 &
cmult -t -a X4 -R R44 -b X4 -c X4RX4 &

#construct Z'RZ
cmult -t -a Z1 -R R11 -b Z1 -c Z1RZ1 &
cmult -t -a Z1 -R R12 -b Z2 -c Z1RZ2 &
cmult -t -a Z1 -R R13 -b Z3 -c Z1RZ3 &
cmult -t -a Z1 -R R14 -b Z4 -c Z1RZ4 &
cmult -t -a Z2 -R R12 -b Z1 -c Z2RZ1 &
cmult -t -a Z2 -R R22 -b Z2 -c Z2RZ2 &
cmult -t -a Z2 -R R23 -b Z3 -c Z2RZ3 &
cmult -t -a Z2 -R R24 -b Z4 -c Z2RZ4 &
cmult -t -a Z3 -R R13 -b Z1 -c Z3RZ1 &
cmult -t -a Z3 -R R23 -b Z2 -c Z3RZ2 &
cmult -t -a Z3 -R R33 -b Z3 -c Z3RZ3 &
cmult -t -a Z3 -R R34 -b Z4 -c Z3RZ4 &
cmult -t -a Z4 -R R14 -b Z1 -c Z4RZ1 &
cmult -t -a Z4 -R R24 -b Z2 -c Z4RZ2 &
cmult -t -a Z4 -R R34 -b Z3 -c Z4RZ3 &
cmult -t -a Z4 -R R44 -b Z4 -c Z4RZ4 &

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

#construct X'RZ
cmult -t -a X1 -R R11 -b Z1 -c X1RZ1 &
cmult -t -a X1 -R R12 -b Z2 -c X1RZ2 &
cmult -t -a X1 -R R13 -b Z3 -c X1RZ3 &
cmult -t -a X1 -R R14 -b Z4 -c X1RZ4 &
cmult -t -a X2 -R R12 -b Z1 -c X2RZ1 &
cmult -t -a X2 -R R22 -b Z2 -c X2RZ2 &
cmult -t -a X2 -R R23 -b Z3 -c X2RZ3 &
cmult -t -a X2 -R R24 -b Z4 -c X2RZ4 &
cmult -t -a X3 -R R13 -b Z1 -c X3RZ1 &
cmult -t -a X3 -R R23 -b Z2 -c X3RZ2 &
cmult -t -a X3 -R R33 -b Z3 -c X3RZ3 &
cmult -t -a X3 -R R34 -b Z4 -c X3RZ4 &
cmult -t -a X4 -R R14 -b Z1 -c X4RZ1 &
cmult -t -a X4 -R R24 -b Z2 -c X4RZ2 &
cmult -t -a X4 -R R34 -b Z3 -c X4RZ3 &
cmult -t -a X4 -R R44 -b Z4 -c X4RZ4 &

#construct Z'RX
cmult -t -a Z1 -R R11 -b X1 -c Z1RX1 &
cmult -t -a Z1 -R R12 -b X2 -c Z1RX2 &
cmult -t -a Z1 -R R13 -b X3 -c Z1RX3 &
cmult -t -a Z1 -R R14 -b X4 -c Z1RX4 &
cmult -t -a Z2 -R R12 -b X1 -c Z2RX1 &
cmult -t -a Z2 -R R22 -b X2 -c Z2RX2 &
cmult -t -a Z2 -R R23 -b X3 -c Z2RX3 &
cmult -t -a Z2 -R R24 -b X4 -c Z2RX4 &
cmult -t -a Z3 -R R13 -b X1 -c Z3RX1 &
cmult -t -a Z3 -R R23 -b X2 -c Z3RX2 &
cmult -t -a Z3 -R R33 -b X3 -c Z3RX3 &
cmult -t -a Z3 -R R34 -b X4 -c Z3RX4 &
cmult -t -a Z4 -R R14 -b X1 -c Z4RX1 &
cmult -t -a Z4 -R R24 -b X2 -c Z4RX2 &
cmult -t -a Z4 -R R34 -b X3 -c Z4RX3 &
cmult -t -a Z4 -R R44 -b X4 -c Z4RX4 &

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
cadd -a Z4RZ4 -r $g44 -b Ainverse -c Z4RZ4.v &
cadd -a Z4RZ3 -r $g43 -b Ainverse -c Z4RZ3.v &
wait

m11=0.2009
m12=-0.0069
m21=-0.0069
m22=0.0031
cadd -a M1RM1 -r $m11 -b Ainverse -c M1RM1.v &
cadd -a M1RM2 -r $m12 -b Ainverse -c M1RM2.v &
cadd -a M2RM1 -r $m21 -b Ainverse -c M2RM1.v &
cadd -a M2RM2 -r $m22 -b Ainverse -c M2RM1.v &
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

#cadd -a X1Ry1 -b X1Ry2 -c rhs.1.1 &
#cadd -a X1Ry3 -b X1Ry4 -c rhs.1.2 &
#cadd -a rhs.1.1 -b rhs.1.2 -c rhs.1
#
#cadd -a X2Ry1 -b X2Ry2 -c rhs.2.1 &
#cadd -a X2Ry3 -b X2Ry4 -c rhs.2.2 &
#cadd -a rhs.2.1 -b rhs.2.2 -c rhs.1
#
#cadd -a X3Ry1 -b X3Ry2 -c rhs.3.1 &
#cadd -a X3Ry3 -b X3Ry4 -c rhs.3.2 &
#cadd -a rhs.3.1 -b rhs.3.2 -c rhs.3
#
#cadd -a X4Ry1 -b X4Ry2 -c rhs.4.1 &
#cadd -a X4Ry3 -b X4Ry4 -c rhs.4.2 &
#cadd -a rhs.4.1 -b rhs.4.2 -c rhs.4
#wait
#
#
#cadd -a Z1Ry1 -b Z1Ry2 -c rhs.5.1 &
#cadd -a Z1Ry3 -b Z1Ry4 -c rhs.5.2 &
#cadd -a rhs.5.1 -b rhs.5.2 -c rhs.5
#
#cadd -a Z2Ry1 -b Z2Ry2 -c rhs.6.1 &
#cadd -a Z2Ry3 -b Z2Ry4 -c rhs.6.2 &
#cadd -a rhs.6.1 -b rhs.6.2 -c rhs.6
#
#cadd -a Z3Ry1 -b Z3Ry2 -c rhs.7.1 &
#cadd -a Z3Ry3 -b Z3Ry4 -c rhs.7.2 &
#cadd -a rhs.7.1 -b rhs.7.2 -c rhs.7
#
#cadd -a Z4Ry1 -b Z4Ry2 -c rhs.8.1 &
#cadd -a Z4Ry3 -b Z4Ry4 -c rhs.8.2 &
#cadd -a rhs.8.1 -b rhs.8.2 -c rhs.8
#
#cadd -a M1Ry1 -b M1Ry2 -c rhs.9.1 &
#cadd -a M1Ry3 -b M1Ry4 -c rhs.9.2 &
#cadd -a rhs.9.1 -b rhs.9.2 -c rhs.9
#
#cadd -a M2Ry1 -b M2Ry2 -c rhs.10.1 &
#cadd -a M2Ry3 -b M2Ry4 -c rhs.10.2 &
#cadd -a rhs.10.1 -b rhs.10.2 -c rhs.10
#
#cadd -a P1Ry1 -b P1Ry2 -c rhs.11.1 &
#cadd -a P1Ry3 -b P1Ry4 -c rhs.11.2 &
#cadd -a rhs.11.1 -b rhs.11.2 -c rhs.11
#
#cadd -a P2Ry1 -b P2Ry2 -c rhs.12.1 &
#cadd -a P2Ry3 -b P2Ry4 -c rhs.12.2 &
#cadd -a rhs.12.1 -b rhs.12.2 -c rhs.12
#
#cvcat rhs.1 rhs.2 rhs.3 rhs.4 rhs.5 rhs.6 rhs.7 rhs.8 rhs.9 rhs.10 rhs.11 rhs.12 rhs

#csolve -A lhs -x solve.sln -b rhs -h -t 1e-15 > solve.log
#pcgmgpu -A lhs -x solve_gpu.sln -b rhs -t 1e-15 > solve_gpu.log
