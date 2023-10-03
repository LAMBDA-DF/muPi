# RPi HQ camera for particle detection

Set of scripts, software and documentation to process images taken with the RPi HQ camera for particle detection.

```mermaid
sequenceDiagram
Title: Data taking and processing chain
Start!           ->  raw DNGs: takeImages.sh
raw DNGs         ->  FITS images: raw2fits.py
raw DNGs         --> FITS images: (procImages.sh)
critical subtract master bias
  FITS images      ->  master bias: computeMedian.exe
  master bias      ->  equalized images: subFits.py
  FITS images      --> equalized images: (subImages.sh)
end
equalized images ->  catalog of tracks: extract.exe
```

+ **takeImages.sh**: script that runs in the RPi to take raw images in DNG format
+ **raw2fits.py**: reads raw data in DNG format and produces a FITS image
+ **subFits.py**: subtracts two FITS images
+ **fits2root.py**: converts a FITS images to a ROOT TTree