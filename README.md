# gpx_to_skiz converter
This app will take all .gpx from the folder it is executed, and convert it to .skiz format compatible with Ski Tracks.
I have no linked with Ski Tracks devs, and am only doing this as I am a huge fan of this super app and want to import gpx generated through other apps (the lead dev, Steve, also offers to do it himself if you send him an email containing gpx files ; but I didn't want to bother him every time).

## Usage
To run, download and execute the .sh script in a linux Ubuntu environement. Here is a command line script for that :
```bash
curl --progress-bar -fL https://raw.githubusercontent.com/alexbelgium/gpx_to_skiz/main/gpx_to_skiz.sh > gpx_to_skiz.sh
chmod +x gpx_to_skiz.sh
./gpx_to_skiz.sh
rm gpx_to_skiz.sh
```

It will convert all gpx files in the run folder.
The converted files will be stored in a new folder named gps_to_skiz. 
Import them into ski tracks. All segment data and statistics should be automatically filled. Otherwise, go in the segments section, and scroll down to trigger a refresh.

## Limitations
- Only Ubuntu supported for the moment
- GPX files don't all have GPS precision embedded. This will be set to 0 if not found.
- I do not know what the last field is but it was 0.25 in all files I checked.  This will be set to 0.25.
