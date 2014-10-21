#ID sire dam CG1 CG2 T1 T2.1 T2.2
#1  . .  a1  a2  50  34  32  
#2  1 .   .  a1   .  42   .
#3  1 2  a2   .  65   .   .
#4  1 3  a2  a3  70   .  35
#5  4 2  a2  a4   .  43   .

awk '{print $7}' data.txt > y.1 
awk '{print $1}' data.txt > id.1

awk '{if($7!="."||$8!=".")print 1;else print "."}' data.txt >y.dummy
awk '{print $1}' data.txt > id.2

awk '{print $1;print $1}' data.txt>id.rep
awk '{print $7;print $8}' data.txt>y.t
cgen_z -d id.rep -e id.2 -r y.t -o D
cnewr -R R -r resid.2x2 y.1 y.dummy

cmult -a D -b R22 -c DR22
cmult -a DR22 -t -b D -c DR22D

awk '{print $1}' data.txt > id.eff
cgen_z -d id.rep -e id.eff -r y.t -o Z2

cmult -t -a Z2 -b DR22D -c Z2DR22D
cmult -a Z2DR22D -b Z2 -c Z2DR22DZ2

cgen_z -d id.1 -e id.eff -r y.1 -o Z1
cmult -t -a Z1 -b DR12 -c Z1DR12
cmult -a Z1DR12 -b Z2 -c Z1DR12Z2


