year=2022
month=9
cycle=000

for ((day=1; day<=26; day+=1))
do

fday=`printf "%02d\n" "${day}"`
fmonth=`printf "%02d\n" "${month}"`

for ((i=0; i<=18; i+=6))
do

ftime=`printf "%02d\n" "${i}"`
server=https://www.ncei.noaa.gov/thredds/fileServer/model-gfs-g4-anl-files
directory=${year}${fmonth}/${year}${fmonth}${fday}
file=gfs_4_${year}${fmonth}${fday}_${ftime}00_${cycle}.grb2


url=${server}/${directory}/${file}

echo $url

wget ${url}

done
done

