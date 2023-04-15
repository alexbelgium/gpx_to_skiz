# gpx_to_skiz
This app will take all .gpx from the folder it is executed, and convert it to .skiz format compatible with Ski Tracks. 

## Usage
To run, download and execute the .sh script in a linux Ubuntu environement. Here is a command line script for that :
```bash
curl -fLO https://raw.githubusercontent.com/alexbelgium/gpx_to_skiz/main/gpx_to_skiz.sh
chmod +x gpx_to_skiz.sh
./gpx_to_skiz.sh
rm gpx_to_skiz.sh
```

It will convert all gpx files in the run folder.

The converted files will be stored in a new folder named gps_to_skiz. 

Import them into ski tracks, go in the segments section, and scroll down to trigger a refresh. This will regenerate all segments data and statistics

## Limitations
- GPX files don't seem to contain gps precision data. This will be set to 0.
- Bearing not supported by gpsbabel. This will be set to 0.
- I do not know what the last field is but it was 0.25 in all files I checked.  This will be set to 0.25.

