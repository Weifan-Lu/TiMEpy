#!/bin/sh
#
#     TidalStrain.2.sh
#
#     2020/06/25 HIROSE F.@MRI
#        Main points of modification since "TidalStrain" are as follows:
#        * Based on gotic2(ver.2004) for ocean tidal loading effect.
#          Note that files with "XXXX_mod2020.f" indicate subroutines that are different from subroutines in gotic2(ver.2004).
#        * We revised "TidalStrain.2/source_mod2/convl3_mod2020.f" because an error case appears depending on computing environment.
#        * We revised "TidalStrain.2/Calc/calcu_idep_ninf.f" & "TidalStrain.2/source_mod2/loadgreenf3_mod.f" because of errors in calculating for an event of shallower than 10 km in deep.
#
#     2024/03/25
#        We replaced select_Lame.f that obtains the PREM parameters and calculates the Lame constant with a more sophisticated code.
#

#date
echo "*****************************************"
echo "**********    TidalStrain.2    **********"
echo "*****************************************"

#-f77.sh

if [[ -f *.all ]]; then
rm *.all
#rm ocean_loading.in.all strn2dVdTdS.prm.all
#rm earthtide.input5.all
fi

data=./event.in

count=1
nevent=3
while [ $count -le $nevent ];
do

#== make fundamental input data ===
ID=`awk -v l=$count 'NR==l {print $1}' $data`
lon=`awk -v l=$count 'NR==l {print $2}' $data`
lat=`awk -v l=$count 'NR==l {print $3}' $data`
dep=`awk -v l=$count 'NR==l {print $4}' $data`
stime=`awk -v l=$count 'NR==l {print $5}' $data`
etime=`awk -v l=$count 'NR==l {print $6}' $data`

echo $ID
echo $dep
sed s/eventID/$ID/ ocean_loading.in.base > buff
sed s/longitude/$lon/ buff > buff2
sed s/latitude/$lat/ buff2 > buff3
sed s/depth/$dep/ buff3 > buff4
sed s/starttime/$stime/ buff4 > buff5
sed s/endtime/$etime/ buff5 > buff6

rm buff buff2 buff3 buff4 buff5
cp buff6 ocean_loading.in

cat buff6 >> ocean_loading.in.all
rm buff6

#== Solid tidal effect ===
#echo "Solid tide                 [by earthtide_mod]"

JA=`awk -v l=$count 'NR==l {print $5}' $data | cut -c1-4`
JB=`awk -v l=$count 'NR==l {print $5}' $data | cut -c5-6`
JC=`awk -v l=$count 'NR==l {print $5}' $data | cut -c7-8`
JD=`awk -v l=$count 'NR==l {print $5}' $data | cut -c9-10`
JE=`awk -v l=$count 'NR==l {print $5}' $data | cut -c11-12`
#PH=`awk -v l=$count 'NR==l {print $3}' $data`
#RM=`awk -v l=$count 'NR==l {print $2}' $data`
#DP=`echo $dep/-1000 | bc -l | cut -c1-4`
DP=`echo $dep/-1000 | bc -l | cut -c1-6`

sed s/JA/$JA/ earthtide.input5.base > buff
sed s/JB/$JB/ buff > buff2
sed s/JC/$JC/ buff2 > buff3
sed s/JD/$JD/ buff3 > buff4
sed s/JE/$JE/ buff4 > buff5
sed s/PH/$lat/ buff5 > buff6
sed s/RM/$lon/ buff6 > buff7
sed s/DP/$DP/ buff7 > buff8

rm buff buff2 buff3 buff4 buff5 buff6 buff7
cp buff8 earthtide.input5
./earthtide_mod < earthtide.input5
# -> solid_earthtide.output9

cat buff8 >> earthtide.input5.all
rm buff8

#== Ocean tidal loading effect ===
#echo "Ocean tide loading effect  [by GOTIC2_mod2]"
#depkm=`echo $dep/-1000 | bc -l | cut -c1-4`
depkm=`echo $dep/-1000 | bc -l | cut -c1-6`
echo $depkm
echo $DP
echo $depkm > depthkm.out

# get idep1, idep2, ninf2
./calcu_idep_ninf
# -> idep_ninf.out

idep1=`awk 'NR==1 {print $1}' idep_ninf.out`
idep2=`awk 'NR==2 {print $1}' idep_ninf.out`
ninf2=`awk 'NR==3 {print $1}' idep_ninf.out`
sed s/IDDDDDD1/$idep1/ loadgreenf3.input5.base > buff
sed s/IDDDDDD2/$idep2/ buff > buff2
sed s/NIIIIII2/$ninf2/ buff2 > buff3

# make input file for loadgreenf3_mod
cp buff3 loadgreenf3.input5
./loadgreenf3_mod < loadgreenf3.input5
cp loadgreenf3_mod.output9 ./data/grn1.data

# calculate strain tensor
#cp buff6 event.in
#./gotic2_mod_2 < event.in
#./gotic2_mod_2_v2 < ocean_loading.in
#./GOTIC2_mod2 < ocean_loading.in
./GOTIC2_mod2020 < ocean_loading.in
sed 's/D/E/g' ocean_loading.out > ocean_loading.out.2


#== memo
#-- time series
#--  Solid: solid_earthtide.output9
#--  Ocean: ocean_loading.out.2

#== Solid + Ocean
./sum_both_list
# -> solid-ocean.out
cp -p solid-ocean.out ./output/solid-ocean.out_$ID

#== delta CFF ===

# select Lame's constant
./select_Lame
# -> select_Lame.out

# strain tensor -> dcff
lambdain=`awk 'NR==1 {print $1}' select_Lame.out`
muin=`awk 'NR==1 {print $2}' select_Lame.out`
strikein=`awk -v l=$count 'NR==l {print $7}' $data`
dipin=`awk -v l=$count 'NR==l {print $8}' $data`
rakein=`awk -v l=$count 'NR==l {print $9}' $data`

sed s/strike1/$strikein/ strn2dVdTdS.prm.base > buff
sed s/dip1/$dipin/ buff > buff2
sed s/rake1/$rakein/ buff2 > buff3
sed s/lambda1/$lambdain/ buff3 > buff4
sed s/mu1/$muin/ buff4 > buff5

rm buff buff2 buff3 buff4

cp buff5 strn2dVdTdS.prm

#fhirose 2015/11/16
#./STRN2DCFF
# -> strn2dcff.out

./STRN2dVdTdS
# -> strn2dVdTdS.out

cat buff5 >> strn2dVdTdS.prm.all
rm buff5

cp -p strn2dVdTdS.out ./output/strn2dVdTdS.out_$ID

count=`expr $count + 1`
done

./conv_dTdS2dCFF.sh
#-> Tidal_response.out

#date
