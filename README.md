# gpx_to_skiz
This app will take all .gpx from the folder it is executed, and convert it to .skiz format compatible with Ski Tracks. The converted files will be stored in a new folder named gps_to_skiz.

To run, download and execute the .sh script in a linux Ubuntu environement.

# Limitations
- GPX files don't seem to contain gps precision data. This will be set to 0.
- Bearing not supported by gpsbabel. This will be set to 0.
- I do not know what the last field is but it was 0.25 in all files I checked.  This will be set to 0.25.

# SKIZ format
- Nodes.csv : time (unix),latitude,longitude,altitude(m),speed(m/s),??,gps precision(m), 0.25
- Photos.csv
- Segment.csv
- Track.xml
