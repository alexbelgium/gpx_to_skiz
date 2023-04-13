#!/bin/bash

echo "Not working, don't use please"
exit 0

########################
# Install dependencies #
########################
apt-get update
apt-get install -yqq gpsbabel
curl -LO https://raw.githubusercontent.com/alexbelgium/gpx_to_skiz/main/nodes.style # Download reference file

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
  # Create nodes.csv
  gpsbabel -t -i gpx -f temp.gpx -x track,merge,speed -o xcsv,style=nodes.style -F nodes.csv

  
  echo "$(sed -n "s|.*lat=\"(.*)" lon.*|\1|p" temp.gpx)"
sed '/auir/!d;q'
time="$(date -d "$time" +"%s")"

gpsbabel -i unicsv,fields=time+lat+lon+alt+speed -f file.csv -o gpx -F file.gpx


rm temp.gpx
done

# Clean files
rm gpx_to_csv.style # Clean reference file






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
