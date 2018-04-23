#!/bin/bash
#
# This script runs calc_image for different times. You can define 
# parameters in this script, the script then creates input files for
# each timeframe with your parameters in them, and submits an array
# jobscript to cheops.
# Needs correct magspot folder to work. Should be in the same folder as
# this script!
########################################################################
#  DEFINITIONS                                                         #
########################################################################

# Black Hole parameters
gM=20.0                                  # BH mass (in solar masses)
a=.5                                   # BH spin parameter
incl=79.0                                 # observer's inclination
gDist=1e3                               # BH distance (in parsec)
# blob/spot properties
Rg=1                # Rg=1: gRadialPos given in units of r_g
                    # Rg=0: in units of ISCO
gRadialPos=15.0      # initial radial position of the blob/spot (in r_g)
gBlobSize=5.0       # size of the blob/spot
startAngle=90.0      # where the blob starts orbit (in deg)
gB=1                                      # 0 = vertical B-field
                                          # 1 = azimuthal B-field
bExp=0.0             # expansion,
zVel=0.05             # z- and radial velocity of the blob
rVel=0.02             # (only works with gModel=5)
# time parameters
startTime=0                            # time for-loop parameters
stopTime=100
timeStep=20
# model
gModel=5                               # 0 = Spot (supposed 2D)
                                       # 1 = Disk
                                       # 2 = Disk + Spot
                                       # 3 = Disk + Jet
                                       # 4 = Blob (3D)
                                       # 5 = model-diskblob2.c 
gThickdisk=0
# code properties
gNx=200                    # = param_nx (resolution)
gStep=0.5                 # = param_step; step along geodesics (GM/c²)
gFrames=100                # = param_frames (number of time frames)
gImg=1                    # 1 => calc_image(), 0 => calc_time_evol()
gRmax=30                # = max domain radius (GM/c2); (param_rmax)
gVis=10
# folder where this script is run and the magspot folder should be                     
ScriptFolder=/afs/ph1.uni-koeln.de/home/khandu/Documents/Munawwar/Raytrace_Dec17/

##### SCRIPT STARTS HERE ###############################################
if [ -f $ScriptFolder/input_img.x ]; then 
gNum=`wc -l "input_img.x" | awk '{print $1'}` 
numStart=`echo $gNum - 1 | bc` 
gNum=`echo $gNum - 2 | bc` 
else
   echo "# gImg gNx gFrames gB gModel gThickdisk a incl Rg gRadialPos gBlobSize gStep gM gDist gRmax gTime startAngle startFrame zVel rVel bExp gVis" > \
$ScriptFolder/input_img.x
gNum=-1;
numStart=0 
fi

                                    # go through these time frames
for (( gTimeSecond = $startTime ; gTimeSecond <= $stopTime ; gTimeSecond += $timeStep ))
do 

gNum=`echo $gNum + 1 | bc`;

#gTime=`echo $gTimeSecond + 0.5000 | bc`
gTime=`echo $gTimeSecond \\* 1 | bc`
#gTime=`echo $gTimeSecond \\* 0.5 | bc`

                                    # create input.x, which is needed
                                    # by magspot to run properly
   echo "$gImg $gNx $gFrames $gB $gModel $gThickdisk $a $incl $Rg $gRadialPos $gBlobSize $gStep $gM $gDist $gRmax $gTime $startAngle 90 $zVel $rVel $bExp $gVis" >> \
$ScriptFolder/input_img.x

   echo "$gImg $gNx $gFrames $gB $gModel $gThickdisk $a $incl $Rg $gRadialPos $gBlobSize $gStep $gM $gDist $gRmax $gTime $startAngle 90 $zVel $rVel $bExp $gVis" > \
$ScriptFolder/input_img.$gNum

done # time for


# run magspot for all new input_img files
echo "Input files created. Starting with magspot."
for (( number = $numStart ; number <= $gNum ; number += 1 ))
do
echo "Running .$ScriptFolder magspot with input_img.$number..."
   ./sphere < input_img.$number > output_img.$number 
echo "Done. Result is in output_img.$number"
done
echo "Finished." 

                                    # submit the jobarray to cheops
#sbatch --array=$numStart-$gNum imgjob.sh

##### EOF ##############################################################
