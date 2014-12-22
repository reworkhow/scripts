#!/bin/bash

OLDDATA="/home/haocheng/Hereford/13Traits_small/data_original" #original data, pedigree with A inverse, G matrix 
NEWDATA="/home/haocheng/Hereford/13Traits_small/data_new" #modified data
SOLUTION="/home/haocheng/Hereford/13Traits_small/solutions" #results
SCRIPT="/home/haocheng/Hereford/13Traits_small/scripts" #results
TEMP="/home/haocheng/Hereford/13Traits_small/temp" #temporary files

##################################################################
##
## make coefficient matrices for trait with 
## direct genetic, maternal genetic, permanent environment effects
##
##################################################################

cd $NEWDATA
# These traits names can be in any order
array1=(bw ww ss fw yw carc_ema carc_imf carc_rib cwt ema imf rib mcw) #traits with direct genetic effects
array2=(bw ww ss) #traits with maternal genetic
array3=(bw ww) #traits with permanent environment effects

awk '{print $1}' stacked_ped >id.eff 
for i in "${array1[@]}"
do
awk '{print $2}' phe100.$i".fill.ped" > y.$i &
awk '{print $1}' phe100.$i".fill.ped" > id.dat.$i &
awk '{print $4}' phe100.$i".fill.ped" > cg1.dat.$i &
awk '$2!="."{print $4}' phe100.$i".fill.ped" | sort -u > cg1.eff.$i &
wait
cgen_z -d cg1.dat.$i -e cg1.eff.$i -r y.$i -o X.$i &
cgen_z -d id.dat.$i -e id.eff -r y.$i -o Z.$i &
done

wait
for i in "${array2[@]}"
do
awk '{print $6}' phe100.$i".fill.ped" > mat.dat.$i
cgen_z -d mat.dat.$i -e id.eff -r y.$i -o M.$i &
done

for i in "${array3[@]}"
do
awk '$2!="."&&$6!="."{print $6}' phe100.$i".fill.ped" | sort -u > perm.eff.$i
cgen_z -d mat.dat.$i -e perm.eff.$i -r y.$i -o P.$i &
done
wait
rm cg1.dat.* id.dat.* mat.dat.*  


###########################################
#
# Make G inverse and R inverse
#
###########################################
#  Note the order of these variables must correspond to the rows and columns of R and G
rm Rinv.*.*

cnewr -R Rinv -r $OLDDATA"/resid.13x13" y.bw y.ww y.yw y.fw y.mcw y.imf y.ema y.rib y.ss y.cwt y.carc_ema y.carc_rib y.carc_imf
array=(bw ww yw fw mcw imf ema rib ss cwt carc_ema carc_rib carc_imf)
for i in $(seq 13)
do
namei=${array[$(($i-1))]}
for j in $(seq $i)
do
namej=${array[$(($j-1))]}
#echo $i $j $namei $namej
mv Rinv$j$i Rinv.$namej.$namei
if [ $i -ne $j ] ; then 
ln -s Rinv.$namej.$namei Rinv.$namei.$namej 
fi 
done
done

# Convert dense to sparse matrix
invert -i $OLDDATA"/bigG.18x18" -o bigG.inv
ident 18 > ident.18
cmult -a bigG.inv -b ident.18 -c bigG.inverted
rm bigG.inv ident.18 

###########################################
#
# Make LHS: step 1
#
###########################################
echo "MAP" > lhs.map 

arrayX=(bw ww yw fw mcw imf ema rib ss cwt carc_ema carc_rib carc_imf)
arrayZ=(bw ww yw fw mcw imf ema rib ss cwt carc_ema carc_rib carc_imf)
arrayM=(bw ww ss)
arrayP=(bw ww)
arrayC=(X Z M P)

for i in ${arrayC[@]}
do
array1=array$i[@]

for j in ${arrayC[@]}
do
array2=array$j[@]

for k in ${!array1}
do
for l in ${!array2}
do

Rinv=Rinv.$k.$l

cmult -t -a $i.$k -R $Rinv -b $j.$l -c $i.$k."Rinv."$j.$l
echo -n $i.$k."Rinv."$j.$l" " >> lhs.map

cmult -t -a $i.$k -R $Rinv -b y.$l -c $i.$k."Rinv.y".$l
done

echo >> lhs.map
done

done


done

###########################################
#
# Make LHS: step 2  Add the inverse var-cov matrix for random effects
#
############################################
array1=($(seq 0 15))
array2=($(seq 16 17))
array3=(Z.bw Z.ww Z.yw Z.fw Z.mcw Z.imf Z.ema Z.rib Z.ss Z.cwt Z.carc_ema Z.carc_rib Z.carc_imf M.bw M.ww M.ss P.bw P.ww)

for i in "${array1[@]}"
do
for j in "${array1[@]}"
do

gInvert=`extract_ij $(($i+1)) $(($j+1)) bigG.inverted` #check later to see gInvert = != 0

cadd -a ${array3[i]}".Rinv."${array3[j]}  -r $gInvert -b Ainverse -c ${array3[i]}".Rinv."${array3[j]}".v"

done
done

idnum=$(awk '$6!="\."{print $6}' phe100.bw.fill.ped|sort -u|wc| awk '{print $1}')  
ident $idnum  > myIdent
for i in "${array2[@]}"
do
for j in "${array2[@]}"
do

gInvert=`extract_ij $(($i+1)) $(($j+1)) bigG.inverted` #check later to see gInvert = != 0

cadd -a ${array3[i]}".Rinv."${array3[j]}  -r $gInvert -b myIdent -c ${array3[i]}".Rinv."${array3[j]}".v"

done
done

##########################################
# 
#  concat RHS (MAP file later) 
#
###########################################
echo "MAP" >rhs.map

array=(X.bw X.ww X.yw X.fw X.mcw X.imf X.ema X.rib X.ss X.cwt X.carc_ema X.carc_rib X.carc_imf Z.bw Z.ww Z.yw Z.fw Z.mcw Z.imf Z.ema Z.rib Z.ss Z.cwt Z.carc_ema Z.carc_rib Z.carc_imf M.bw M.ww M.ss P.bw P.ww)

for i in "${array[@]}" #X,Z..
do

cadd -a $i".Rinv.y.bw" -b $i".Rinv.y.ww" -c rhs.1.1 &
cadd -a $i".Rinv.y.yw" -b $i".Rinv.y.fw" -c rhs.1.2 &
cadd -a $i".Rinv.y.mcw" -b $i".Rinv.y.imf" -c rhs.1.3 &
cadd -a $i".Rinv.y.ema" -b $i".Rinv.y.rib" -c rhs.1.4 &
cadd -a $i".Rinv.y.ss" -b $i".Rinv.y.cwt" -c rhs.1.5 &
cadd -a $i".Rinv.y.carc_ema" -b $i".Rinv.y.carc_rib" -c rhs.1.6 &
wait
cadd -a rhs.1.1 -b rhs.1.2 -c rhs.1.7 &
cadd -a rhs.1.3 -b rhs.1.4 -c rhs.1.8 &
cadd -a rhs.1.5 -b rhs.1.6 -c rhs.1.9 &
wait
cadd -a rhs.1.7 -b rhs.1.8 -c rhs.1.10 &
cadd -a rhs.1.9 -b $i".Rinv.y.carc_imf" -c rhs.1.11 &
wait
cadd -a rhs.1.10 -b rhs.1.11 -c rhs.$i
echo "rhs.$i" >> rhs.map
done

rm rhs.1.*

#######################################
#
# Make Map files
#
#######################################


#######################################
#
#solve the equations
#
#######################################

csolve -A $SCRIPT"/lhs.map" -x $SOLUTION"/solve_13trait_1219.sln" -b rhs -h > $SOLUTION"/solve_13trait_1219.log"

mv X* Z* M* P* R* rhs* y* $TEMP
rm myIdent

