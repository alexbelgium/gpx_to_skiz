# gpx_to_skiz
This app will take all .gpx from the folder it is executed, and convert it to .skiz format compatible with Ski Tracks. The converted files will be stored in a new folder named gps_to_skiz.

To run, download and execute the .sh script in a linux Ubuntu environement.

# Limitations
- GPX files do not contain gps precision data. They will appear as 0 in all converted skiz files.
- I do not know what the last field is. It will therefore be fixed at 0.25 for all files

# SKIZ format
- Nodes.csv : time (unix),latitude,longitude,altitude(m),speed(m/s),??,gps precision(m), 0.25
- Photos.csv
- Segment.csv
- Track.xml
