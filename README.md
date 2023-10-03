# RPi HQ camera for particle detection

Set of scripts and codes to process images taken with the RPi HQ camera for particle detection

```sequence {theme=hand}
Title: Basic workflow
Start!           ->  raw DNGs: takeImages.sh
raw DNGs         ->  FITS images: raw2fits.py
raw DNGs         --> FITS images: (raw2fits.py)
FITS images      ->  median (master bias): checkConsistencyAndComputeMedian.exe
master bias      ->  equalized images: subtract.exe
FITS images      ->  equalized images: 
equalized images ->  catalog of tracks: extract.exe
```

+ *takeImages.sh*: script to take raw images in DNG format in the RPi
+ *raw2fits.py*: reads raw data in DNG format and produces a FITS image
