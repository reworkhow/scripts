#!/bin/bash

#####birth weight
awk '{print $2}' phe.bw.fill.ped > y.bw &
awk '{print $1}' phe.bw.fill.ped > id.dat.bw &
awk '{print $1}' stacked_ped > id.eff &
awk '{print $4}' phe.bw.fill.ped > cg1.dat.bw &
awk '{print $4}' phe.bw.fill.ped | sort -u > cg1.eff.bw &
wait
cgen_z -d cg1.dat.bw -e cg1.eff.bw -r y.bw -o X.bw
cgen_z -d id.dat.bw -e id.eff -r y.bw -o Z
wait

#maternal effect
awk '{print $1}' stacked_ped > mat.eff.bw &
awk '{print $6}' phe.bw.fill.ped > mat.dat.bw &
wait
cgen_z -d mat.dat.bw -e mat.eff.bw -r y.bw -o M

#permanent env effect
awk '$6!=.{print $6}' phe.bw.fill.ped | sort -u > perm.eff.bw &
awk '{print $6}' phe.bw.fill.ped > perm.dat.bw &
wait
cgen_z -d perm.dat.bw -e perm.eff.bw -r y.bw -o P


#####weaning weight
awk '{print $2}' phe.ww.fill.ped > y.ww &
awk '{print $4}' phe.ww.fill.ped > cg1.dat.ww &
awk '{print $4}' phe.ww.fill.ped | sort -u > cg1.eff.ww &
wait
cgen_z -d cg1.dat.ww -e cg1.eff.ww -r y.ww -o X.ww
wait


#final weight
awk '{print $2}' phe.fw.fill.ped > y.fw &
awk '{print $4}' phe.fw.fill.ped > cg1.dat.fw &
awk '{print $4}' phe.fw.fill.ped | sort -u > cg1.eff.fw &
wait
cgen_z -d cg1.dat.fw -e cg1.eff.fw -r y.fw -o X.fw
wait

#yearling weight
awk '{print $2}' phe.yw.fill.ped > y.yw &
awk '{print $4}' phe.yw.fill.ped > cg1.dat.yw &
awk '{print $4}' phe.yw.fill.ped | sort -u > cg1.eff.yw &
wait
cgen_z -d cg1.dat.yw -e cg1.eff.yw -r y.yw -o X.yw
wait

cat y.bw >y1 &
cat y.ww >y2 &
cat y.yw >y3 &
cat y.fw >y4 &
cat X.bw >X1 &
cat X.ww >X2 &
cat X.yw >X3 &
cat X.fw >X4 &

cat Z.bw > Z &
cat M.bw > M &
cat P.bw > P &
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
cmult -t -a Z -b Z -c ZPZ
ckron -a ZPZ -b Rinverse  -c ZRZ &

#construct Zm'RZm
cmult -t -a M -b M -c MPM
ckron -a MPM -b Rinverse -c MRM &

#construct Zp'RZp
cmult -t -a P -b P -c PPP
ckron -a PPP -b Rinverse -c PRP &

###Multiplication related to X cannot use kroneker, since X1!=X2!=X3
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

###kronecker

#construct Z'RZm
cmult -t -a Z -b M -c ZPM
ckron -a ZPM -b Rinverse -c ZRM
#construct Zm'RZ
cmult -t -a M -b Z -c MPZ
ckron -a MPZ -b Rinverse -c MRZ

#construct Z'RZp
cmult -t -a Z -b P -c ZPP
ckron -a ZPP -b Rinverse -c ZRP
#construct Zp'RZ
cmult -t -a P -b Z -c PPZ
ckron -a PPZ -b Rinverse -c PRZ


#construct Zm'RZp
cmult -t -a M -b P -c MPP
ckron -a MPP -b Rinverse -c MRP
#construct Zp'RZm
cmult -t -a P -b M -c PPM
ckron -a PPM -b Rinverse -c PRM




wait

g11=30.00
g12=75.64
g13=128.45
g14=88.40
g21=75.64
g22=679.00
g23=905.31
g24=817.80
g31=128.45
g32=905.31
g33=1886.00
g34=1557.61
g41=88.40
g42=817.80
g43=1557.61
g44=2010.00


cadd -a Z1RZ1 -r 30.00 -b Ainverse -c Z1RZ1.v
cadd -a Z1RZ2 -r 75.64 -b Ainverse -c Z1RZ2.v
cadd -a Z1RZ3 -r 128.45 -b Ainverse -c Z1RZ3.v
cadd -a Z1RZ4 -r 88.40 -b Ainverse -c Z1RZ4.v
cadd -a Z2RZ1 -r 75.64 -b Ainverse -c Z2RZ1.v
cadd -a Z2RZ2 -r 679.00 -b Ainverse -c Z2RZ2.v
cadd -a Z3RZ1 -r 128.45 -b Ainverse -c Z3RZ1.v
cadd -a Z2RZ3 -r 905.31 -b Ainverse -c Z2RZ3.v
cadd -a Z4RZ1 -r 88.40 -b Ainverse -c Z4RZ1.v
cadd -a Z2RZ4 -r 817.80 -b Ainverse -c Z2RZ4.v
cadd -a Z3RZ2 -r 905.31 -b Ainverse -c Z3RZ2.v
cadd -a Z3RZ3 -r 1886.00 -b Ainverse -c Z3RZ3.v
cadd -a Z3RZ4 -r 1557.61 -b Ainverse -c Z3RZ4.v
cadd -a Z4RZ2 -r 817.80 -b Ainverse -c Z4RZ2.v
cadd -a Z4RZ4 -r 2010.00 -b Ainverse -c Z4RZ4.v
cadd -a Z4RZ3 -r 1557.61 -b Ainverse -c Z4RZ3.v
wait

m11=5.40
m12=1
m21=1
m22=1
cadd -a M1RM1 -r $m11 -b Ainverse -c Z1RZ1.v &
cadd -a M1RM2 -r $m12 -b Ainverse -c Z1RZ2.v &
cadd -a M2RM1 -r $m21 -b Ainverse -c Z1RZ3.v &
cadd -a M2RM2 -r $m22 -b Ainverse -c Z1RZ4.v &
wait


p11=5.40
p12=1
p21=1
p22=1
cadd -a P1RP1 -r $p11 -b Ainverse -c P1RP1.v &
cadd -a P1RP2 -r $p12 -b Ainverse -c P1RP2.v &
cadd -a P2RP1 -r $p21 -b Ainverse -c P1RP3.v &
cadd -a P2RP2 -r $p22 -b Ainverse -c P1RP4.v &
wait


chcat X1RX1 X1RX2 X1RX3 X1RX4 X1RZ1 X1RZ2 X1RZ3 X1RZ4 lhs.1
chcat X2RX1 X2RX2 X2RX3 X2RX4 X2RZ1 X2RZ2 X2RZ3 X2RZ4 lhs.2
chcat X3RX1 X3RX2 X3RX3 X3RX4 X3RZ1 X3RZ2 X3RZ3 X3RZ4 lhs.3
chcat X4RX1 X4RX2 X4RX3 X4RX4 X4RZ1 X4RZ2 X4RZ3 X4RZ4 lhs.4
chcat Z3RX1 Z3RX2 Z3RX3 Z3RX4 Z3RZ1.v Z3RZ2.v Z3RZ3.v Z3RZ4.v lhs.7
chcat Z2RX1 Z2RX2 Z2RX3 Z2RX4 Z2RZ1.v Z2RZ2.v Z2RZ3.v Z2RZ4.v lhs.6
chcat Z1RX1 Z1RX2 Z1RX3 Z1RX4 Z1RZ1.v Z1RZ2.v Z1RZ3.v Z1RZ4.v lhs.5
chcat Z4RX1 Z4RX2 Z4RX3 Z4RX4 Z4RZ1.v Z4RZ2.v Z4RZ3.v Z4RZ4.v lhs.8
wait
cvcat lhs.1 lhs.2 lhs.3 lhs.4 lhs.5 lhs.6 lhs.7 lhs.8 lhs

cmult -t -a X1 -R R11 -b y1 -c X1Ry1
cmult -t -a X1 -R R12 -b y2 -c X1Ry2
cmult -t -a X2 -R R22 -b y2 -c X2Ry2
cmult -t -a X2 -R R23 -b y3 -c X2Ry3
cmult -t -a X1 -R R13 -b y3 -c X1Ry3
cmult -t -a X3 -R R13 -b y1 -c X3Ry1
cmult -t -a X1 -R R14 -b y4 -c X1Ry4
cmult -t -a X4 -R R14 -b y1 -c X4Ry1
cmult -t -a Z4 -R R24 -b y2 -c Z4Ry2
cmult -t -a Z2 -R R23 -b y3 -c Z2Ry3
cmult -t -a X4 -R R24 -b y2 -c X4Ry2
cmult -t -a Z1 -R R14 -b y4 -c Z1Ry4
cmult -t -a X4 -R R34 -b y3 -c X4Ry3
cmult -t -a X4 -R R44 -b y4 -c X4Ry4
cmult -t -a Z3 -R R23 -b y2 -c Z3Ry2
cmult -t -a Z3 -R R13 -b y1 -c Z3Ry1
cmult -t -a Z1 -R R11 -b y1 -c Z1Ry1
cmult -t -a Z1 -R R12 -b y2 -c Z1Ry2
cmult -t -a Z3 -R R33 -b y3 -c Z3Ry3
cmult -t -a X3 -R R34 -b y4 -c X3Ry4
cmult -t -a X3 -R R33 -b y3 -c X3Ry3
cmult -t -a X3 -R R23 -b y2 -c X3Ry2
cmult -t -a X2 -R R24 -b y4 -c X2Ry4
cmult -t -a Z4 -R R44 -b y4 -c Z4Ry4
cmult -t -a X2 -R R12 -b y1 -c X2Ry1
cmult -t -a Z4 -R R34 -b y3 -c Z4Ry3
cmult -t -a Z1 -R R13 -b y3 -c Z1Ry3
cmult -t -a Z2 -R R22 -b y2 -c Z2Ry2
cmult -t -a Z4 -R R14 -b y1 -c Z4Ry1
cmult -t -a Z3 -R R34 -b y4 -c Z3Ry4
cmult -t -a Z2 -R R24 -b y4 -c Z2Ry4
cmult -t -a Z2 -R R12 -b y1 -c Z2Ry1
wait

cadd -a X1Ry1 -b X1Ry2 -c rhs.1.1
cadd -a X1Ry3 -b X1Ry4 -c rhs.1.2
cadd -a rhs.1.1 -b rhs.1.2 -c rhs.1
cadd -a X2Ry1 -b X2Ry2 -c rhs.2.1
cadd -a X2Ry3 -b X2Ry4 -c rhs.2.2
cadd -a rhs.2.1 -b rhs.2.2 -c rhs.2
cadd -a X3Ry1 -b X3Ry2 -c rhs.3.1
cadd -a X3Ry3 -b X3Ry4 -c rhs.3.2
cadd -a rhs.3.1 -b rhs.3.2 -c rhs.3
cadd -a X4Ry1 -b X4Ry2 -c rhs.4.1
cadd -a X4Ry3 -b X4Ry4 -c rhs.4.2
cadd -a rhs.4.1 -b rhs.4.2 -c rhs.4
cadd -a Z1Ry1 -b Z1Ry2 -c rhs.5.1
cadd -a Z1Ry3 -b Z1Ry4 -c rhs.5.2
cadd -a rhs.5.1 -b rhs.5.2 -c rhs.5
cadd -a Z2Ry1 -b Z2Ry2 -c rhs.6.1
cadd -a Z2Ry3 -b Z2Ry4 -c rhs.6.2
cadd -a rhs.6.1 -b rhs.6.2 -c rhs.6
cadd -a Z3Ry1 -b Z3Ry2 -c rhs.7.1
cadd -a Z3Ry3 -b Z3Ry4 -c rhs.7.2
cadd -a rhs.7.1 -b rhs.7.2 -c rhs.7
cadd -a Z4Ry1 -b Z4Ry2 -c rhs.8.1
cadd -a Z4Ry3 -b Z4Ry4 -c rhs.8.2
cadd -a rhs.8.1 -b rhs.8.2 -c rhs.8
wait
cvcat rhs.1 rhs.2 rhs.3 rhs.4 rhs.5 rhs.6 rhs.7 rhs.8 rhs

#csolve -A lhs -x solve.sln -b rhs -h -t 1e-15 > solve.log
#pcgmgpu -A lhs -x solve_gpu.sln -b rhs -t 1e-15 > solve_gpu.log
