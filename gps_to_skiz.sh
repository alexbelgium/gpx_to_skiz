#!/bin/bash

########################
# Install dependencies #
########################
apt-get update
apt-get install -yqq gpsbabel

#############################################
# Convert gpx to individual skiz components #
#############################################

# For all files
for input in /share/*.gpx, do
  # Extract filename
  filename="${input%.*}"
  # Create directory
  mkdir "$filename"

# Create nodes.csv
##################
  # duplicate template
  cp "$input" temp.gpx
  # Remove until first data point
  sed '1,/^<name>$/d' temp.gpx
  # Convert all times to epoch
#TODO
  # Create nodes.csv 
  gpsbabel -t -i gpx -f temp.gpx -x track,merge,speed,rptdigits=6 -o unicsv,fields=time+lat+lon+alt+speed -F nodes.csv


  
  
  echo "$(sed -n "s|.*lat=\"(.*)" lon.*|\1|p" temp.gpx)"
sed '/auir/!d;q'
time="$(date -d "$time" +"%s")"

gpsbabel -i unicsv,fields=time+lat+lon+alt+speed -f file.csv -o gpx -F file.gpx


rm temp.gpx
done

1. Rename .skiz to .zip

2. Create Nodes.csv

- Delete file until first "<trkpt lat="

- - Extract first lat, round to 6 after decimal
- Extract first lon, round to 6 after decimal
- Extract altitude, round to 2 after decimal
- 0.0
- 0.00
- 

"$time,$lat,$lon,$alt,$speed,$??;$??,$??"
1681210036.999,46.007251,6.691240,1723.58,0.0,0.00,24.000000,0.250000


Segment.csv
Track.xml
