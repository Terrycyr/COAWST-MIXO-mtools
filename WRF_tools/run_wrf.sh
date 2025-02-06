#!/bin/bash
export ROOT_DIR="/home/ychen/WRF-master/"
export RUN_DIR="$ROOT_DIR/WRF/run"
export Project_DIR="$PWD"
export GFS_DIR="/home/ychen/WRF-GFS-Hurricane/GFS_data"
# !!! namelist files for following two steps should be prepared first!!!
echo "------------------------------------------------------------------------"
echo "Please select from among the following options:"
echo "1. Complete Run"
echo "2. Run Preprocessing"
echo "3. Clean ALL + Run WRF"
echo "4. Clean Last Run"
echo "5. Clean Preprocessing"

read -p "Enter your selection: " select </dev/tty


############################################################
#                                                          #
#    STEP I - Generate Input Files (Boundaries ...)        #
#                                                          #
############################################################
if [ $select -eq 1 ] || [ $select -eq 2 ]
then

cd $Project_DIR/preprocessing/
rm ./met_em*
rm ./FILE*
rm ./PFILE*
rm ./geo*.nc
rm ./GRIBFILE*
rm *.log
rm log.*
ln -sf $ROOT_DIR/WPS/geogrid.exe ./
ln -sf $ROOT_DIR/WPS/geogrid ./
ln -sf $ROOT_DIR/WPS/ungrib.exe ./
ln -sf $ROOT_DIR/WPS/ungrib ./
ln -sf $ROOT_DIR/WPS/metgrid.exe ./
ln -sf $ROOT_DIR/WPS/metgrid ./
ln -sf $ROOT_DIR/WPS/link_grib.csh

./geogrid.exe | tee -a log.geogrid
./link_grib.csh $GFS_DIR/gfs
ln -sf ungrib/Variable_Tables/Vtable.GFS Vtable
./ungrib.exe
./metgrid.exe | tee -a log.metgrid

fi

if [ $select -eq 5 ]
then
cd $Project_DIR/preprocessing/
rm ./met_em*
rm ./FILE*
rm ./PFILE*
rm ./geo*.nc
rm ./GRIBFILE*
rm *.log
rm log.*
fi

############################################################
#                                                          #
#    STEP II - Run WRF (REAL-TIME CASES)                   #
#                                                          #
############################################################
cd $Project_DIR

if [ $select -eq 1 ] || [ $select -eq 3 ]
then
rm ./wrfinput*
rm ./wrfbdy*
rm ./wrfout*
rm ./wrfrst*
rm ./rsl.*
rm ./*.dat
rm ./*.output
rm ./slurm*

##ln -sf $ROOT_DIR/WRF/run/* ./
ln -sf ./preprocessing/met_em* ./
ln -sf $ROOT_DIR/WRF/run/real.exe ./
ln -sf $ROOT_DIR/WRF/run/wrf.exe ./

##############   Required Input Files  #####################
cp  $RUN_DIR/CAMtr_volume_mixing_ratio ./
cp  $RUN_DIR/LANDUSE.TBL ./
cp  $RUN_DIR/ozone_plev.formatted ./
cp  $RUN_DIR/ozone_lat.formatted ./
cp  $RUN_DIR/ozone.formatted ./
cp  $RUN_DIR/RRTMG_LW_DATA ./
cp  $RUN_DIR/RRTMG_SW_DATA ./
cp  $RUN_DIR/RRTM_DATA ./
cp  $RUN_DIR/VEGPARM.TBL ./
cp  $RUN_DIR/SOILPARM.TBL ./
cp  $RUN_DIR/GENPARM.TBL ./
cp  $RUN_DIR/URBPARM.TBL ./
############################################################

./real.exe

rm ./met_em*

fi

if [ $select -eq 4 ]
then
rm ./wrfout*
rm ./wrfrst*
rm ./rsl.*
rm ./*.dat
rm ./*.output
rm ./slurm*
fi

#sbatch my_job_wrf.sh
