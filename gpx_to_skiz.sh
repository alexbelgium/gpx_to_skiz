#!/bin/bash

clear -x

echo "###########################################"
echo "               .GPX TO .SKIZ               "
echo "               by Alexbelgium              "
echo "###########################################"
echo ""
echo "This app will take all .gpx from the folder it is executed, and convert it to .skiz format compatible with Ski Tracks"
echo "The converted files will be stored in a new folder named gps_to_skiz"
echo ""
echo "https://github.com/alexbelgium/gpx_to_skiz"
echo ""
echo "###########################################"

########################
# Install dependencies #
########################
echo "Installing dependencies, please wait"
if command -v "apk" &>/dev/null; then
    # If apk based
    apk add zip --no-cache >/dev/null
    apk add gpsbabel --no-cache --repository=https://dl-cdn.alpinelinux.org/alpine/edge/community >/dev/null
elif command -v "apt" &>/dev/null; then
    # If apt-get based
    apt-get update >/dev/null
    apt-get install -yqq zip >/dev/null
    apt-get install -yqq gpsbabel >/dev/null
else
    echo "Filesystem not supported, please use alpine or ubuntu"
    exit 1
fi

# Prepare directory
echo "Making gpx_to_skiz directory"
mkdir -p gpx_to_skiz

# Download reference file
echo "Downloading nodes.style"
curl -f -L -s -S https://raw.githubusercontent.com/alexbelgium/gpx_to_skiz/main/helper/nodes.style --output gpx_to_skiz/nodes.style
echo "Downloading template"
curl -f -L -s -S https://raw.githubusercontent.com/alexbelgium/gpx_to_skiz/main/helper/Track.xml --output gpx_to_skiz/Track.xml

#############################################
# Convert gpx to individual skiz components #
#############################################

# For all files
for input in *.gpx; do
  # Text
  echo "Converting $input to skiz format"
  # Extract filename
  filename="${input%.*}"
  # Create directory
  mkdir -p gpx_to_skiz/"$filename"
  # Create Nodes.csv
  gpsbabel -t -i gpx -f "$input" -x track,merge,speed -o xcsv,style=gpx_to_skiz/nodes.style -F gpx_to_skiz/"$filename"/Nodes.csv
  # Create Photos.csv
  touch gpx_to_skiz/"$filename"/Photos.csv
  # Create Segment.csv
  touch gpx_to_skiz/"$filename"/Segment.csv
  # Create Track.xml
  cp gpx_to_skiz/Track.xml gpx_to_skiz/"$filename"/Track.xml
  # Create Track.xml - UID
  sed -i "s=VAR_parseObjectId=$(echo $RANDOM | md5sum | head -c 20; echo;)=g" gpx_to_skiz/"$filename"/Track.xml
  # Create Track.xml - Times
  START_TIME="$(cat gpx_to_skiz/"$filename"/Nodes.csv | awk -F "," '{ print $1 }' | head -1 | sed 's/^/@/' | xargs date +"%Y-%m-%dT%H:%M:%S%:z" -d)"
  END_TIME="$(cat gpx_to_skiz/"$filename"/Nodes.csv | awk -F "," '{ print $1 }' | tail -1 | sed 's/^/@/' | xargs date +"%Y-%m-%dT%H:%M:%S%:z" -d)"
  sed -i "s=VAR_start=$START_TIME=g" gpx_to_skiz/"$filename"/Track.xml
  sed -i "s=VAR_end=$END_TIME=g" gpx_to_skiz/"$filename"/Track.xml
  sed -i "s=VAR_duration=$(($(date -d ${END_TIME} +%s)-$(date -d ${START_TIME} +%s))).000=g" gpx_to_skiz/"$filename"/Track.xml
  sed -i "s=VAR_tz=${START_TIME: -3}=g" gpx_to_skiz/"$filename"/Track.xml
  # Create Track.xml - Speeds
  sed -i "s=VAR_maxspeed=$(cut -d, -f4,4 < gpx_to_skiz/"$filename"/Nodes.csv | sort -nr | head -1)=g" gpx_to_skiz/"$filename"/Track.xml
  # Create skiz
  zip -j -r gpx_to_skiz/"$filename".skiz gpx_to_skiz/"$filename"/*
  # Remove temporary folder
  rm -r gpx_to_skiz/"$filename"
done

echo "All file converted"
