#!/bin/bash
# select a basin, and use warmed-up hydrological model to conduct the IRC.

basinid=500
spininfo='S2'

estart=2

if [ "$estart" -eq 0 ]; then
# start from scratch, from event #1
starteve=1
elif [ "$estart" -eq 1 ]; then
# start from current event, this is due to power being cut off, computer down, and similar issues etc, as IRC is time consuming
starteve=`cat /shared/dondo/home/ml423/world_mts_runs/IRC1/Basin*_"$basinid".txt`
elif [ "$estart" -eq 2 ]; then
# start from next event, this is due to current event does not meet quality control
starteve=`cat /shared/dondo/home/ml423/world_mts_runs/IRC1/Basin*_"$basinid".txt`
starteve=$((starteve+1))
fi

# basin locations 
# 1-132 466-522 are in the Alps
# 133-434 are in the Andes
# 450-487 are in Asia

# copy files from backup folder to current working folder, so multiple IRC can be executed simultaneously without interfering each other.
utch=0

printf "%d" "$basinid" >basin_code_rep85.txt
printf "%d" "$basinid" >current_basin_rep85.txt
printf "%d" "$utch" >utchh_rep85.txt

cp /shared/dondo/home/ml423/world_mts_runs/backup/gauge_event_info.m -r /shared/dondo/home/ml423/world_mts_runs/
cp /shared/dondo/home/ml423/world_mts_runs/backup/check_events.m -r /shared/dondo/home/ml423/world_mts_runs/
cp /shared/dondo/home/ml423/world_mts_runs/backup/get_events.m -r /shared/dondo/home/ml423/world_mts_runs/

mv check_events.m check_events_rep85.m
mv get_events.m get_events_rep85.m
mv gauge_event_info.m gauge_event_info_rep85.m

sed -i 's/basin_code/basin_code_rep85/g' check_events_rep85.m
sed -i 's/events_number/events_number_rep85/g' check_events_rep85.m
sed -i 's/basin_gauge_id/basin_gauge_id_rep85/g' check_events_rep85.m

sed -i 's/basin_code/basin_code_rep85/g' get_events_rep85.m
sed -i 's/event_code/event_code_rep85/g' get_events_rep85.m
sed -i 's/datnm/datnm_rep85/g' get_events_rep85.m
sed -i 's/icdaten/icdaten_rep85/g' get_events_rep85.m
sed -i 's/event_qc/event_qc_rep85/g' get_events_rep85.m
sed -i 's/basin_code/basin_code_rep85/g' gauge_event_info_rep85.m
sed -i 's/gauge_number/gauge_number_rep85/g' gauge_event_info_rep85.m
sed -i 's/too_short/too_short_rep85/g' gauge_event_info_rep85.m

# pull out basin information and check if it satisfys quality control.
matlab -nodisplay -nosplash -nodesktop -r "run /shared/dondo/home/ml423/world_mts_runs/check_events_rep85.m"
matlab -nodisplay -nosplash -nodesktop -r "run /shared/dondo/home/ml423/world_mts_runs/gauge_event_info_rep85.m"
tooshort=`cat too_short_rep85.txt`

if [ "$tooshort" -eq 1 ]; then
  echo 'The gauge records is less than a year post 19730901, or drainage area too big, skip to next one!'
  exit
  continue
fi

# set up the hydrological model.

totevents=`cat events_number_rep85.txt`
gid=`cat basin_gauge_id_rep85.txt`

echo "$totevents is the total number of events for this basin $gid"
cadjvalue=`cat cadj_"$gid".txt`

cp /shared/dondo/home/ml423/world_mts_runs/backup/DCHM_900m1hr_any_event_IRCICC.f -r /shared/dondo/home/ml423/world_mts_runs/
mv DCHM_900m1hr_any_event_IRCICC.f DCHM_900m1hr.f

sed -i "s/Cadj=0/Cadj=$cadjvalue/g" DCHM_900m1hr.f

# loop through the events.

for evet in `seq $starteve $totevents` 

do
#mail -s "rep85 starts event $evet " mochil@illinois.edu  <<< "basin is $gid, index is $basinid"
printf "%d" "$evet" > /shared/dondo/home/ml423/world_mts_runs/IRC1/Basin"$gid"_"$basinid".txt

printf "%d" "$evet" >event_code_rep85.txt
matlab -nodisplay -nosplash -nodesktop -r "run /shared/dondo/home/ml423/world_mts_runs/get_events_rep85.m"
datenm=`cat datnm_rep85.txt`
lastday=`cat icdaten_rep85.txt`
eqc=`cat event_qc_rep85.txt`


if [[ "$eqc" -ne 1 ]]; then
echo "event quality is no good, peak too ahead or too late within the 30-day window, due to multiple peaks in the window, or IRC already being done"
continue
fi

# load events window points defined by the hydrograph as in Liao and Barros 2022.
evetindex=$(printf %03d $evet)
w2d=`sed -n '1p' /shared/dondo/home/ml423/world_mts_runs/Basins_events_chars/Basin"$gid"_event"$datenm".txt`
w2u=`sed -n '2p' /shared/dondo/home/ml423/world_mts_runs/Basins_events_chars/Basin"$gid"_event"$datenm".txt`
w3u=`sed -n '3p' /shared/dondo/home/ml423/world_mts_runs/Basins_events_chars/Basin"$gid"_event"$datenm".txt`
w4u=`sed -n '4p' /shared/dondo/home/ml423/world_mts_runs/Basins_events_chars/Basin"$gid"_event"$datenm".txt`

# constrained number of event to 150 to avoid large sparse matrix in calculations. This number is sufficiently high. 99% of basins don't have 150 events.
if [ "$evet" -lt 150 ]; then

# update IRC window settings based on Liao and Barros 2022, revise file names for parallel executions

cp /shared/dondo/home/ml423/world_mts_runs/backup/Backtrack_IRC.m -r /shared/dondo/home/ml423/world_mts_runs/
cp /shared/dondo/home/ml423/world_mts_runs/backup/Backtrack_IRC.m -r /shared/dondo/home/ml423/world_mts_runs/
sed -i "s/w2bd1 = 17/w2bd1=$w2d/g" Backtrack_IRC.m
sed -i "s/w2bd = 28/w2bd=$w2u/g" Backtrack_IRC.m
sed -i "s/w3bd1 = 28/w3bd1=$w2u/g" Backtrack_IRC.m
sed -i "s/w3bd = 43/w3bd=$w3u/g" Backtrack_IRC.m
sed -i "s/w4bd1 = 43/w4bd1=$w3u/g" Backtrack_IRC.m
sed -i "s/w4bd = 59/w4bd=$w4u/g" Backtrack_IRC.m

sed -i "s/w2bd1 = 17/w2bd1=$w2d/g" Backtrack_IRC.m
sed -i "s/w2bd = 28/w2bd=$w2u/g" Backtrack_IRC.m
sed -i "s/w3bd1 = 28/w3bd1=$w2u/g" Backtrack_IRC.m
sed -i "s/w3bd = 43/w3bd=$w3u/g" Backtrack_IRC.m
sed -i "s/w4bd1 = 43/w4bd1=$w3u/g" Backtrack_IRC.m
sed -i "s/w4bd = 59/w4bd=$w4u/g" Backtrack_IRC.m



cp /shared/dondo/home/ml423/world_mts_runs/backup/Slow_recession.m -r /shared/dondo/home/ml423/world_mts_runs/
cp /shared/dondo/home/ml423/world_mts_runs/backup/Pre_rising.m -r /shared/dondo/home/ml423/world_mts_runs/
cp /shared/dondo/home/ml423/world_mts_runs/backup/more_iteration_IRC.m -r /shared/dondo/home/ml423/world_mts_runs/
cp /shared/dondo/home/ml423/world_mts_runs/backup/compare_WORLD_30day.m -r /shared/dondo/home/ml423/world_mts_runs/
cp /shared/dondo/home/ml423/world_mts_runs/backup/rainfall_update.m -r /shared/dondo/home/ml423/world_mts_runs/
cp /shared/dondo/home/ml423/world_mts_runs/backup/transition_IRC.m -r /shared/dondo/home/ml423/world_mts_runs/


mv Backtrack_IRC.m backtrack_rep85.m
mv Backtrack_IRC.m backtrackfollow_rep85.m
mv rainfall_update.m rainupv2ori_rep85.m
mv transition_IRC.m transitionb1_rep85.m

mv Slow_recession.m Slow_recession_rep85.m
mv Pre_rising.m Pre_rising_rep85.m
mv more_iteration_IRC.m more_itra_rep85.m
mv compare_WORLD_30day.m compare_rep85.m





sed -i 's/filenum/filenum_rep85/g' backtrack_rep85.m
sed -i 's/filenum/filenum_rep85/g' backtrackfollow_rep85.m
sed -i 's/filenum/filenum_rep85/g' rainupv2ori_rep85.m
sed -i 's/filenum/filenum_rep85/g' transitionb1_rep85.m
sed -i 's/filenum/filenum_rep85/g' Slow_recession_rep85.m
sed -i 's/flow_dis_str/flow_dis_str_rep85/g' Slow_recession_rep85.m
sed -i 's/filenum/filenum_rep85/g' Pre_rising_rep85.m
sed -i 's/datnm/datnm_rep85/g' backtrack_rep85.m
sed -i 's/datnm/datnm_rep85/g' backtrackfollow_rep85.m
sed -i 's/datnm/datnm_rep85/g' rainupv2ori_rep85.m
sed -i 's/datnm/datnm_rep85/g' transitionb1_rep85.m
sed -i 's/datnm/datnm_rep85/g' Slow_recession_rep85.m
sed -i 's/datnm/datnm_rep85/g' Pre_rising_rep85.m
sed -i 's/datnm/datnm_rep85/g' more_itra_rep85.m
sed -i 's/datnm/datnm_rep85/g' compare_rep85.m
sed -i 's/smkge/smkge_rep85/g' compare_rep85.m
sed -i 's/bl/bl_rep85/g' backtrack_rep85.m
sed -i 's/bl/bl_rep85/g' backtrackfollow_rep85.m
sed -i 's/bl/bl_rep85/g' rainupv2ori_rep85.m
sed -i 's/bl/bl_rep85/g' transitionb1_rep85.m
sed -i 's/bl/bl_rep85/g' Slow_recession_rep85.m
sed -i 's/bl/bl_rep85/g' Pre_rising_rep85.m
sed -i 's/xtmpbb1/xtmpbb1_rep85/g' backtrack_rep85.m
sed -i 's/xtmpbb1/xtmpbb1_rep85/g' backtrackfollow_rep85.m
sed -i 's/xtmpbb1/xtmpbb1_rep85/g' rainupv2ori_rep85.m
sed -i 's/xtmpbb1/xtmpbb1_rep85/g' transitionb1_rep85.m
sed -i 's/xtmpbb1/xtmpbb1_rep85/g' Slow_recession_rep85.m
sed -i 's/xtmpbb1/xtmpbb1_rep85/g' Pre_rising_rep85.m
sed -i 's/rainfallb1/rainfallb1_rep85/g' backtrack_rep85.m
sed -i 's/rainfallb1/rainfallb1_rep85/g' backtrackfollow_rep85.m
sed -i 's/rainfallb1/rainfallb1_rep85/g' rainupv2ori_rep85.m
sed -i 's/rainfallb1/rainfallb1_rep85/g' transitionb1_rep85.m
sed -i 's/rainfallb1/rainfallb1_rep85/g' Slow_recession_rep85.m
sed -i 's/rainfallb1/rainfallb1_rep85/g' Pre_rising_rep85.m
sed -i 's/current_basin/current_basin_rep85/g' backtrack_rep85.m
sed -i 's/current_basin/current_basin_rep85/g' backtrackfollow_rep85.m
sed -i 's/current_basin/current_basin_rep85/g' rainupv2ori_rep85.m
sed -i 's/current_basin/current_basin_rep85/g' transitionb1_rep85.m
sed -i 's/current_basin/current_basin_rep85/g' compare_rep85.m
sed -i 's/current_basin/current_basin_rep85/g' more_itra_rep85.m
sed -i 's/current_basin/current_basin_rep85/g' Pre_rising_rep85.m
sed -i 's/current_basin/current_basin_rep85/g' Slow_recession_rep85.m
sed -i 's/utchh/utchh_rep85/g' backtrack_rep85.m
sed -i 's/utchh/utchh_rep85/g' backtrackfollow_rep85.m
sed -i 's/utchh/utchh_rep85/g' rainupv2ori_rep85.m
sed -i 's/utchh/utchh_rep85/g' transitionb1_rep85.m
sed -i 's/utchh/utchh_rep85/g' Pre_rising_rep85.m
sed -i 's/utchh/utchh_rep85/g' Slow_recession_rep85.m


sed -i 's/filu/filu_rep85/g' backtrack_rep85.m
sed -i 's/filu/filu_rep85/g' backtrackfollow_rep85.m
sed -i 's/filu/filu_rep85/g' rainupv2ori_rep85.m
sed -i 's/bd.txt/bd_rep85.txt/g' backtrack_rep85.m
sed -i 's/bd1.txt/bd1_rep85.txt/g' backtrack_rep85.m
sed -i 's/bd.txt/bd_rep85.txt/g' backtrackfollow_rep85.m
sed -i 's/bd1.txt/bd1_rep85.txt/g' backtrackfollow_rep85.m
sed -i 's/loop.txt/loop_rep85.txt/g' rainupv2ori_rep85.m
sed -i 's/loop2.txt/loop2_rep85.txt/g' rainupv2ori_rep85.m


sed -i 's/qualflag.txt/qualflag_rep85.txt/g' rainupv2ori_rep85.m
sed -i 's/pasn.txt/pasn_rep85.txt/g' rainupv2ori_rep85.m
sed -i 's/stopit.txt/stopit_rep85.txt/g' rainupv2ori_rep85.m
sed -i 's/pernn/pernn_rep85/g' rainupv2ori_rep85.m
sed -i 's/intervalmulti/intervalmulti_rep85/g' rainupv2ori_rep85.m


sed -i 's/loop.txt/loop_rep85.txt/g' Slow_recession_rep85.m
sed -i 's/loop2.txt/loop2_rep85.txt/g' Slow_recession_rep85.m
sed -i 's/qualflag.txt/qualflag_rep85.txt/g' Slow_recession_rep85.m
sed -i 's/pasn.txt/pasn_rep85.txt/g' Slow_recession_rep85.m
sed -i 's/stopit.txt/stopit_rep85.txt/g' Slow_recession_rep85.m
sed -i 's/inflection_pt.txt/inflection_pt_rep85.txt/g' Slow_recession_rep85.m
sed -i 's/pre_rise_pt.txt/pre_rise_pt_rep85.txt/g' Slow_recession_rep85.m
sed -i 's/loop.txt/loop_rep85.txt/g' Pre_rising_rep85.m
sed -i 's/loop2.txt/loop2_rep85.txt/g' Pre_rising_rep85.m
sed -i 's/qualflag.txt/qualflag_rep85.txt/g' Pre_rising_rep85.m
sed -i 's/pasn.txt/pasn_rep85.txt/g' Pre_rising_rep85.m
sed -i 's/stopit.txt/stopit_rep85.txt/g' Pre_rising_rep85.m
sed -i 's/rising_pt.txt/rising_pt_rep85.txt/g' Pre_rising_rep85.m
sed -i 's/inflection_pt.txt/inflection_pt_rep85.txt/g' backtrack_rep85.m
sed -i 's/rising_pt.txt/rising_pt_rep85.txt/g' backtrack_rep85.m
sed -i 's/pre_rise_pt.txt/pre_rise_pt_rep85.txt/g' backtrack_rep85.m
sed -i 's/inflection_pt.txt/inflection_pt_rep85.txt/g' backtrackfollow_rep85.m
sed -i 's/rising_pt.txt/rising_pt_rep85.txt/g' backtrackfollow_rep85.m
sed -i 's/pre_rise_pt.txt/pre_rise_pt_rep85.txt/g' backtrackfollow_rep85.m

sed -i 's/pre_rise_pt.txt/pre_rise_pt_rep85.txt/g' rainupv2ori_rep85.m
sed -i 's/pre_rise_pt.txt/pre_rise_pt_rep85.txt/g' Pre_rising_rep85.m
sed -i 's/first_sttime.txt/first_sttime_rep85.txt/g' backtrack_rep85.m
sed -i 's/first_sttime.txt/first_sttime_rep85.txt/g' backtrackfollow_rep85.m
sed -i 's/moreitr.txt/moreitr_rep85.txt/g' more_itra_rep85.m

sed -i 's/loop.txt/loop_rep85.txt/g' backtrack_rep85.m
sed -i 's/loop2.txt/loop2_rep85.txt/g' backtrack_rep85.m
sed -i 's/loop.txt/loop_rep85.txt/g' backtrackfollow_rep85.m
sed -i 's/loop2.txt/loop2_rep85.txt/g' backtrackfollow_rep85.m

#sed -i 's/16)/12)/g' backtrack_rep85.m

printf "%d" "$datenm" >datenm_rep85.txt
filenm_num=4140
filenm=$((datenm * 100 - filenm_num * 100000))
middlef=$((datenm - filenm_num * 1000))
spinupsuffix='c'

# set a folder for temporarily storing Initial conditions 
mkdir /shared/dondo/home/ml423/world_mts_outputs/tmp_Basin"$gid"IC/
mkdir /shared/dondo/home/ml423/world_mts_outputs/tmp_Basin"$gid"IC/"$datenm"_IC/
# copy spin up results And make dir for the first day of IC in the 30-day window and IC right before the rainfall event
mkdir /shared/dondo/home/ml423/world_mts_outputs/Basin"$gid"_output_900m1hr_IC/
mkdir /shared/dondo/home/ml423/world_mts_outputs/Basin"$gid"_output_900m1hr_IC/"$datenm""$spinupsuffix"/
mkdir /shared/dondo/home/ml423/world_mts_outputs/Basin"$gid"_output_900m1hr_IC/"$datenm""$spinupsuffix"/TrueIC/
mkdir /shared/dondo/home/ml423/world_mts_runs/IRC1/ICall/
mkdir /shared/dondo/home/ml423/world_mts_runs/IRC1/ICall/Basin"$gid"/
mkdir /shared/dondo/home/ml423/world_mts_runs/IRC1/ICall/Basin"$gid"/"$filenm_num"
mkdir /shared/dondo/home/ml423/world_mts_runs/IRC1/ICall/Basin"$gid"/"$filenm_num"/"$datenm"

# check if spinup are done for this event, otherwise skip to the next event, this is to ignore events that have NaN values
if [ -f "/shared/dondo/home/ml423/world_mts_outputs/Basin"$gid"_output_900m1hr_"$spininfo"_spinup/"$lastday"/laststep_ws" ]; then
  echo "IC File exists."
else
  echo "IC File does not exist."
  continue
fi

cp /shared/dondo/home/ml423/world_mts_outputs/Basin"$gid"_output_900m1hr_"$spininfo"_spinup/"$lastday"/laststep* -r /shared/dondo/home/ml423/world_mts_outputs/Basin"$gid"_output_900m1hr_IC/"$datenm""$spinupsuffix"/

echo "$datenm is the target event for $gid"
# copy target date rainfall STIVDBKC, or etc
mkdir /shared/dondo/home/ml423/world_mts_outputs/Basin"$gid"_input_ERA5_900m_tmp/
mkdir /shared/dondo/home/ml423/world_mts_outputs/Basin"$gid"_input_ERA5_900m_tmp/"$datenm"/
cp /shared/dondo/home/ml423/world_mts/Precip/Event/Basin"$gid"_input_ERA5_900m_FIX/"$datenm"/* -r /shared/dondo/home/ml423/world_mts/Basin"$gid"_input_ERA5_900m_tmp/"$datenm"/

# load basin's characteristics, number of rows, cols, and npoin, and outlet row and col, so the hydrological model is suited to this individual basin

row=`sed -n '1p' /shared/dondo/home/ml423/world_mts/Basins_events_chars/Basin"$gid".txt`
col=`sed -n '2p' /shared/dondo/home/ml423/world_mts/Basins_events_chars/Basin"$gid".txt`
npoin=`sed -n '3p' /shared/dondo/home/ml423/world_mts/Basins_events_chars/Basin"$gid".txt`
rg=`sed -n '4p' /shared/dondo/home/ml423/world_mts/Basins_events_chars/Basin"$gid".txt`
cg=`sed -n '5p' /shared/dondo/home/ml423/world_mts/Basins_events_chars/Basin"$gid".txt`

sed -i "s/nrow=di1,ncol=di2,npoin=di3,nrg=di4,ncg=di5/nrow=$row,ncol=$col,npoin=$npoin,nrg=$rg,ncg=$cg/g" DCHM_900m1hr.f
sed -i "s/Basinxx/Basin$gid/g" DCHM_900m1hr.f
sed -i "s/XXX/11397/g" DCHM_900m1hr.f
mpif77 -mcmodel=medium -fno-automatic -o DCHM_900m1hr_Basin"$gid".exe DCHM_900m1hr.f
sed -i "s/11397/XXX/g" DCHM_900m1hr.f

mpirun -np 2 ./DCHM_900m1hr_Basin"$gid".exe "$datenm" '20140231' '5' '2' '95' 'SW' '721'


mkdir /shared/dondo/home/ml423/world_mts_runs/IRC1/Basin"$gid"outputs/
mkdir /shared/dondo/home/ml423/world_mts_runs/IRC1/Basin"$gid"outputs/refxtmp"$filenm"fullresults
cp /shared/dondo/home/ml423/world_mts_outputs/Basin"$gid"_output_900m1hr/"$datenm"_SW -r /shared/dondo/home/ml423/world_mts_runs/IRC1/Basin"$gid"outputs/refxtmp"$filenm"fullresults

# starts the classic 3 iterations run of the IRC.

for bl in $(seq 0 1 8)
do
printf "%d" "$bl" >bl_rep85.txt
printf "%d" "$filenm" >filenum_rep85.txt
printf "%d" "$filenm" >filenumbeforechange_rep85.txt

matlab -nodisplay -nosplash -nodesktop -r "run /shared/dondo/home/ml423/world_mts_runs/backtrack_rep85.m"
realICloc=`cat pre_rise_pt_rep85.txt`


filenm=`cat filu_rep85.txt`


mkdir /shared/dondo/home/ml423/world_mts_runs/IRC1/Basin"$gid"outputs/refxtmp"$filenm"fullresults
mkdir /shared/dondo/home/ml423/world_mts_runs/IRC1/Basin"$gid"outputs/refxtmp$filenm


matlab -nodisplay -nodesktop -r "run /shared/dondo/home/ml423/world_mts_runs/transitionb1_rep85.m"

echo "$bl is this value bl"
echo "$realICloc is this value realICloc"
echo "$evet is this value evet"

if [[ "$bl" -eq 0 ]]; then
echo "critical important here"
mpirun -np 2 ./DCHM_900m1hr_Basin"$gid".exe "$datenm" '20140231' '5' '2' '95' 'SW' "$realICloc"
#The script above should only be executed once, and this is to get the intermediate outputs regarding IC, and will be updated later on 
#first transfer the intermediate output (the IC conditions right before the event) to another folder
cp /shared/dondo/home/ml423/world_mts_outputs/Basin"$gid"_output_900m1hr/"$datenm"_SW/midout/laststep* -r /shared/dondo/home/ml423/world_mts_outputs/Basin"$gid"_output_900m1hr_IC/"$datenm""$spinupsuffix"/TrueIC/

fi

# 0 3 6 correspond to Window 2 in the classic 4-window IRC setup.
if [ "$bl" -eq 0 -o "$bl" -eq 3 -o "$bl" -eq 6 ]; then

matlab -nodisplay -nodesktop -r "run /shared/dondo/home/ml423/world_mts_runs/Pre_rising_rep85.m"
echo "$bl is this ICCpre"
# rerun with updated IC
mpirun -np 2 ./DCHM_900m1hr_Basin"$gid".exe "$datenm" '20140231' '58' '2' '95' 'SW' "$realICloc"
# 58 and 5 as istart is the same, just to distinguish whether or not read outputs from middleoutputs, which are now True IC right before the rainall events
fi


# if a window size is too big, requiring too much hourly IRC execution, coarsen the IRC technique to 2-hourly to improve speed, this strategy does not show a 
# significant change on the outcomes, but dramatically reduces computation time.
bd=`cat bd_rep85.txt`
bd1=`cat bd1_rep85.txt`

bdiff=$((bd-bd1))
echo "$bdiff"
modull=$(( bdiff % 2 ))

if [ "$bdiff" -gt 6 ]; then
squ=$(seq $bd -2 $bd1)
rainupmulti=2
else
squ=$(seq $bd -1 $bd1)
rainupmulti=1
fi
printf "%d" "$rainupmulti" >intervalmulti_rep85.txt

if [ "$bdiff" -gt 6 ] && [ "$modull" -eq 1 ] ; then
# the starting and the ending has to be bd and bd1, so the command below is to impose it
squ+=" $bd1"
fi

# loop through each step in the window that needs to do IRC

for loo in $squ
# {2..96}
do
printf "%d" "$loo" >loop_rep85.txt
loo2=0
printf "%d" "$loo2" >loop2_rep85.txt
qualflag=99
printf "%d" "$qualflag" >qualflag_rep85.txt


echo "$evet is this value evet"
echo "$bl is this IRCori"



matlab -nodisplay -nodesktop -r "run /shared/dondo/home/ml423/world_mts_runs/rainupv2ori_rep85.m"


value=`cat pasn_rep85.txt`

if [[ "$bl" -ge 7 ]]; then
# added 03/06/2025 to reduce data output sizes, only save one slice of intermediate results rather than full iterations of all windows results
mkdir /shared/dondo/home/ml423/world_mts_runs/IRC1/Basin"$gid"outputs/refxtmp"$filenm"fullresults/iteration$loo-$loo2
cp /shared/dondo/home/ml423/world_mts_outputs/Basin"$gid"_output_900m1hr/"$datenm"_SW -r /shared/dondo/home/ml423/world_mts_runs/IRC1/Basin"$gid"outputs/refxtmp"$filenm"fullresults/iteration$loo-$loo2
fi

for loo2 in $(seq 1 1 1)
do
printf "%d" "$loo2" >loop2_rep85.txt
value=`cat pasn_rep85.txt`

if [ "$loo2" -gt 1 ]; then
stov=`cat stopit_rep85.txt`
if [ "$stov" -eq 0 ]; then
echo "Continue iterations"
else
echo "Need to stop: minimum corrections not satisfied"
break
fi
fi

if [[ "$value" -eq 0 ]]; then
# before doing anything, every new run need to reset the xtmp in matlab and also rainfall in the input


mpirun -np 2 ./DCHM_900m1hr_Basin"$gid".exe "$datenm" '20140231' '58' '2' '95' 'SW' "$realICloc"
echo "This rep script rep85"
printf "%d" "$loo" >loop_rep85.txt



echo "loo value is $loo"
echo "iteration value is $loo2"
qualflag=100
printf "%d" "$qualflag" >qualflag_rep85.txt

echo "$bl is this IRCori"
matlab -nodisplay -nodesktop -r "run /shared/dondo/home/ml423/world_mts_runs/rainupv2ori_rep85.m"


else
printf "%d" "$loo" >loop_rep85.txt
echo "loo value is $loo"
qualflag=99
printf "%d" "$qualflag" >qualflag_rep85.txt
break
fi

if [[ "$bl" -ge 7 ]]; then
# added 03/06/2025 to reduce data output sizes, these are intermediate results
mkdir /shared/dondo/home/ml423/world_mts_runs/IRC1/Basin"$gid"outputs/refxtmp"$filenm"fullresults/iteration$loo-$loo2
cp /shared/dondo/home/ml423/world_mts_outputs/Basin"$gid"_output_900m1hr/"$datenm"_SW -r /shared/dondo/home/ml423/world_mts_runs/IRC1/Basin"$gid"outputs/refxtmp"$filenm"fullresults/iteration$loo-$loo2
fi

# after each time step, re-calculate travel time distributions and do a new tracking for next time step IRC. this turns out to be not sensitive and not needed.
matlab -nodisplay -nosplash -nodesktop -r "run /shared/dondo/home/ml423/world_mts_runs/backtrackfollow_rep85.m"

done
done


# 2 5 8 correspond to the recession window in the classic IRC scheme.
# this section below is to consider the numerical errors accumulated when tracking time is long, from begining of rainfall to the start of baseflow domination point
# this depends on each event, but this step turns out to not have any impacts on the results. 

if [ "$bl" -eq 2 -o "$bl" -eq 5 -o "$bl" -eq 8 ]; then

for icitr in 0 1 
#icitr is essentially pixel distance to a stream pixel
do
cp /shared/dondo/home/ml423/world_mts_outputs/Basin"$gid"_output_900m1hr_IC/"$datenm""$spinupsuffix"/TrueIC/laststep* -r /shared/dondo/home/ml423/world_mts_outputs/tmp_Basin"$gid"IC/"$datenm"_IC/

printf "%d" "$icitr" >flow_dis_str_rep85.txt

matlab -nodisplay -nodesktop -r "run /shared/dondo/home/ml423/world_mts_runs/Slow_recession_rep85.m"
echo "$bl is this ICC"
# rerun with updated IC
mpirun -np 2 ./DCHM_900m1hr_Basin"$gid".exe "$datenm" '20140231' '58' '2' '95' 'SV' "$realICloc"
# compare SW and SV if SW is better, then abandon the changes and move to next icitr, if SV is better, then adopt the changes and move to next icitr.
matlab -nodisplay -nodesktop -r "run /shared/dondo/home/ml423/world_mts_runs/compare_rep85.m"

kgev=`cat smkge_rep85.txt`
if [ "$kgev" -eq 1 ]; then
#sm change makes KGE get better, adopt this change, and move to next icitr
rm -r /shared/dondo/home/ml423/world_mts_outputs/Basin"$gid"_output_900m1hr/"$datenm"_SW
echo "remove worse IC corrections"
mv /shared/dondo/home/ml423/world_mts_outputs/Basin"$gid"_output_900m1hr/"$datenm"_SV /shared/dondo/home/ml423/world_mts_outputs/Basin"$gid"_output_900m1hr/"$datenm"_SW
else
# reverse IC changes, and then move to next icitr
cp /shared/dondo/home/ml423/world_mts_outputs/tmp_Basin"$gid"IC/"$datenm"_IC/laststep* -r /shared/dondo/home/ml423/world_mts_outputs/Basin"$gid"_output_900m1hr_IC/"$datenm""$spinupsuffix"/TrueIC/
fi

done
# every whole iteration save the IC just to avoid accidental interruption
cp /shared/dondo/home/ml423/world_mts_outputs/Basin"$gid"_output_900m1hr_IC/"$datenm""$spinupsuffix"/TrueIC/ -r /shared/dondo/home/ml423/world_mts_runs/IRC1/ICall/Basin"$gid"/"$filenm_num"/"$datenm"

fi


if [[ "$bl" -ge 7 ]]; then
# added 03/06/2025 to reduce data output sizes, these are intermediate results
cp /shared/dondo/home/ml423/world_mts_outputs/Basin"$gid"_output_900m1hr/"$datenm"_SW -r /shared/dondo/home/ml423/world_mts_runs/IRC1/Basin"$gid"outputs/refxtmp"$filenm"fullresults
fi

done


rm -rf /shared/dondo/home/ml423/world_mts_runs/IRC1/Basin"$gid"outputs/refxtmp"$middlef"03fullresults
rm -rf /shared/dondo/home/ml423/world_mts_runs/IRC1/Basin"$gid"outputs/refxtmp"$middlef"04fullresults
rm -rf /shared/dondo/home/ml423/world_mts_runs/IRC1/Basin"$gid"outputs/refxtmp"$middlef"12fullresults
rm -rf /shared/dondo/home/ml423/world_mts_runs/IRC1/Basin"$gid"outputs/refxtmp"$middlef"13fullresults
rm -rf /shared/dondo/home/ml423/world_mts_runs/IRC1/Basin"$gid"outputs/refxtmp"$middlef"14fullresults
rm -rf /shared/dondo/home/ml423/world_mts_runs/IRC1/Basin"$gid"outputs/refxtmp"$middlef"22fullresults
rm -rf /shared/dondo/home/ml423/world_mts_runs/IRC1/Basin"$gid"outputs/refxtmp"$middlef"23fullresults

matlab -nodisplay -nodesktop -r "run /shared/dondo/home/ml423/world_mts_runs/more_itra_rep85.m"

moreitr=`cat moreitr_rep85.txt`

moreitr=0
# more iterations are not needed for the paper because of good performance by only 3 iterations.  

if [ "$moreitr" -eq 0 ]; then

# after all IRC iterations done, copy the initial conditions over
mkdir /shared/dondo/home/ml423/world_mts_runs/IRC1/ICall/
mkdir /shared/dondo/home/ml423/world_mts_runs/IRC1/ICall/Basin"$gid"/
mkdir /shared/dondo/home/ml423/world_mts_runs/IRC1/ICall/Basin"$gid"/"$filenm_num"
mkdir /shared/dondo/home/ml423/world_mts_runs/IRC1/ICall/Basin"$gid"/"$filenm_num"/"$datenm"
cp /shared/dondo/home/ml423/world_mts_outputs/Basin"$gid"_output_900m1hr_IC/"$datenm""$spinupsuffix"/TrueIC/ -r /shared/dondo/home/ml423/world_mts_runs/IRC1/ICall/Basin"$gid"/"$filenm_num"/"$datenm"

# run the model for one last time after iteration 3 final window, this is to get results after the third IC change.
filenm=$((datenm * 100 - filenm_num * 100000 + 32))
# artificially set it as next iteration window 2.! careful!!!
mkdir /shared/dondo/home/ml423/world_mts_runs/IRC1/Basin"$gid"outputs/refxtmp"$filenm"fullresults
mpirun -np 2 ./DCHM_900m1hr_Basin"$gid".exe "$datenm" '20140231' '58' '2' '95' 'SW' "$realICloc"
cp /shared/dondo/home/ml423/world_mts_outputs/Basin"$gid"_output_900m1hr/"$datenm"_SW -r /shared/dondo/home/ml423/world_mts_runs/IRC1/Basin"$gid"outputs/refxtmp"$filenm"fullresults

else

for bl in $(seq 9 1 14)
do
printf "%d" "$bl" >bl_rep85.txt
printf "%d" "$filenm" >filenum_rep85.txt

matlab -nodisplay -nosplash -nodesktop -r "run /shared/dondo/home/ml423/world_mts_runs/backtrack_rep85.m"
#call execute_command_line(" matlab -nodesktop -nojvm -r 'backtrack_IRCv2.m; exit' ")


filenm=`cat filu_rep85.txt`


mkdir /shared/dondo/home/ml423/world_mts_runs/IRC1/Basin"$gid"outputs/refxtmp"$filenm"fullresults
mkdir /shared/dondo/home/ml423/world_mts_runs/IRC1/Basin"$gid"outputs/refxtmp$filenm


matlab -nodisplay -nodesktop -r "run /shared/dondo/home/ml423/world_mts_runs/transitionb1_rep85.m"

mpirun -np 2 ./DCHM_900m1hr_Basin"$gid".exe "$datenm" '20140231' '58' '2' '95' 'SW' "$realICloc"

# the below is to run ICC first before any IRC
if [ "$bl" -eq 9 -o "$bl" -eq 12 ]; then

matlab -nodisplay -nodesktop -r "run /shared/dondo/home/ml423/world_mts_outputs/Pre_rising_rep85.m"
echo "$bl is this ICCpre"
# rerun with updated IC
mpirun -np 2 ./DCHM_900m1hr_Basin"$gid".exe "$datenm" '20140231' '58' '2' '95' 'SW' "$realICloc"
fi



bd=`cat bd_rep85.txt`
bd1=`cat bd1_rep85.txt`
for loo in $(seq $bd -1 $bd1)
# {2..96}
do
printf "%d" "$loo" >loop_rep85.txt
loo2=0
printf "%d" "$loo2" >loop2_rep85.txt
qualflag=99
printf "%d" "$qualflag" >qualflag_rep85.txt

# ML add 2024-05-27
#matlab -nodisplay -nodesktop -r "run /shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/ICup_b1_250_rep85.m"
# ended add

echo "$bl is this IRCori"
matlab -nodisplay -nodesktop -r "run /shared/dondo/home/ml423/world_mts_runs/rainupv2ori_rep85.m"




value=`cat pasn_rep85.txt`

mkdir /shared/dondo/home/ml423/world_mts_runs/IRC1/Basin"$gid"outputs/refxtmp"$filenm"fullresults/iteration$loo-$loo2
cp /shared/dondo/home/ml423/world_mts_outputs/Basin"$gid"_output_900m1hr/"$datenm"_SW -r /shared/dondo/home/ml423/world_mts_runs/IRC1/Basin"$gid"outputs/refxtmp"$filenm"fullresults/iteration$loo-$loo2

for loo2 in $(seq 1 1 1)
do
printf "%d" "$loo2" >loop2_rep85.txt
value=`cat pasn_rep85.txt`

if [ "$loo2" -gt 1 ]; then
stov=`cat stopit_rep85.txt`
if [ "$stov" -eq 0 ]; then
echo "Continue iterations"
else
echo "Need to stop: minimum corrections not satisfied"
break
fi
fi

if [ "$value" -eq 0 ]; then
# before doing anything, every new run need to reset the xtmp in matlab and also S4xxxx rainfall in the input
#mpif77 -mcmodel=medium -fno-automatic -o LSHM_IPHExHindcast_250m5min_Basin02b2fixZB2.exe LSHM_IPHExHindcast_250m5min_Basin02b2fixZB2.f
#mpirun -np 2 ./LSHM_IPHExHindcast_250m5min_Basin01c2.exe '2014051400' '2014051400' '1' '1' '1' 'SW'
#mpif77 -mcmodel=medium -fno-automatic -o LSHM_IPHExForecast_250m5min_Basin01.exe LSHM_IPHExForecast_250m5min_Basin01.f
#mpirun -np 2 ./LSHM_IPHExForecast_250m5min_Basin01.exe '2014050300' '2014050200' '2' '1' 'SW'

mpirun -np 2 ./DCHM_900m1hr_Basin"$gid".exe "$datenm" '20140231' '58' '2' '95' 'SW' "$realICloc"
#mpirun -np 2 ./LSHM_IPHExForecast_250m5min_Basin01.exe '2014051600' '2014051500' '2' '2' 'SW'

printf "%d" "$loo" >loop_rep85.txt
echo "loo value is $loo"
echo "iteration value is $loo2"
qualflag=100
printf "%d" "$qualflag" >qualflag_rep85.txt
# below commented out 2014-01-27 No
#matlab -nodisplay -nodesktop -r "run /shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/ICup_b1_250_rep85.m"


echo "$bl is this IRCori"
matlab -nodisplay -nodesktop -r "run /shared/dondo/home/ml423/world_mts_runs/rainupv2ori_rep85.m"



else
printf "%d" "$loo" >loop_rep85.txt
echo "loo value is $loo"
qualflag=99
printf "%d" "$qualflag" >qualflag_rep85.txt
break
fi


mkdir /shared/dondo/home/ml423/world_mts_runs/IRC1/Basin"$gid"outputs/refxtmp"$filenm"fullresults/iteration$loo-$loo2
cp /shared/dondo/home/ml423/world_mts_outputs/Basin"$gid"_output_900m1hr/"$datenm"_SW -r /shared/dondo/home/ml423/world_mts_runs/IRC1/Basin"$gid"outputs/refxtmp"$filenm"fullresults/iteration$loo-$loo2

done
done

#
if [ "$bl" -eq 11 -o "$bl" -eq 14 ]; then

for icitr in 0 1 2 3 4 5
#icitr is essentially pixel distance to a stream pixel
do

cp /shared/dondo/home/ml423/world_mts_outputs/Basin"$gid"_output_900m1hr_IC/"$datenm""$spinupsuffix"/TrueIC/laststep* -r /shared/dondo/home/ml423/world_mts_outputs/tmp_Basin"$gid"IC/"$datenm"_IC/

printf "%d" "$icitr" >flow_dis_str_rep85.txt

matlab -nodisplay -nodesktop -r "run /shared/dondo/home/ml423/world_mts_outputs/Slow_recession_rep85.m"
echo "$bl is this ICC"
# rerun with updated IC
mpirun -np 2 ./DCHM_900m1hr_Basin"$gid".exe "$datenm" '20140231' '58' '2' '95' 'SW' "$realICloc"
# compare SW and SV if SW is better, then abandon the changes and move to next icitr, if SV is better, then adopt the changes and move to next icitr.
matlab -nodisplay -nodesktop -r "run /shared/dondo/home/ml423/world_mts_outputs/compare_rep85.m"
kgev=`cat smkge_rep85.txt`


if [ "$kgev" -eq 1 ]; then
#sm change makes KGE get better, adopt this change, and move to next icitr
rm -r /shared/dondo/home/ml423/world_mts_outputs/Basin"$gid"_output_900m1hr/"$datenm"_SW
echo "remove worse IC corrections"
mv /shared/dondo/home/ml423/world_mts_outputs/Basin"$gid"_output_900m1hr/"$datenm"_SV /shared/dondo/home/ml423/world_mts_outputs/Basin"$gid"_output_900m1hr/"$datenm"_SW
else
# reverse IC changes, and then move to next icitr
cp /shared/dondo/home/ml423/world_mts_outputs/tmp_Basin"$gid"IC/"$datenm"_IC/laststep* -r /shared/dondo/home/ml423/world_mts_outputs/Basin"$gid"_output_900m1hr_IC/"$datenm""$spinupsuffix"/TrueIC/
fi

done

# every whole iteration save the IC just to avoid accidental interruption
cp /shared/dondo/home/ml423/world_mts_outputs/Basin"$gid"_output_900m1hr_IC/"$datenm""$spinupsuffix"/TrueIC/ -r /shared/dondo/home/ml423/world_mts_runs/IRC1/ICall/Basin"$gid"/"$filenm_num"/"$datenm"
fi

cp /shared/dondo/home/ml423/world_mts_outputs/Basin"$gid"_output_900m1hr/"$datenm"_SW -r /shared/dondo/home/ml423/world_mts_runs/IRC1/Basin"$gid"outputs/refxtmp"$filenm"fullresults


done

# after all IRC iterations done, copy the initial conditions over

mkdir /shared/dondo/home/ml423/world_mts_runs/IRC1/ICall/
mkdir /shared/dondo/home/ml423/world_mts_runs/IRC1/ICall/Basin"$gid"/
mkdir /shared/dondo/home/ml423/world_mts_runs/IRC1/ICall/Basin"$gid"/"$filenm_num"
mkdir /shared/dondo/home/ml423/world_mts_runs/IRC1/ICall/Basin"$gid"/"$filenm_num"/"$datenm"
cp /shared/dondo/home/ml423/world_mts_outputs/Basin"$gid"_output_900m1hr_IC/"$datenm""$spinupsuffix"/TrueIC/ -r /shared/dondo/home/ml423/world_mts_runs/IRC1/ICall/Basin"$gid"/"$filenm_num"/"$datenm"

## run the model for one last time after iteration 3 final window, this is to get results after the third IC change.
#filenm=$((datenm * 100 - filenm_num * 100000 + 52))
## artificially set it as next iteration window 2.! careful!!!
#mkdir /shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/refxtmp"$filenm"fullresults
#mpirun -np 2 ./LSHM_IPHExHindcast_250m5min_Basin01"$datenm".exe "$datenm" '20140231' '5' '1' '5' 'SW'
#cp /shared/dondo/home/ml423/HMIOP/NoDAre/LB/Basin01_output_250m5min_"$datenm"_SW -r /shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/refxtmp"$filenm"fullresults
#

# fi below is for wethear need more iterations for IRC
fi

fi

done