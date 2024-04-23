#!python
import sys,os
import modulib
import essio

if len(sys.argv) > 1:
    infile = sys.argv[1]
else:
    print("Input data file path is not provided")

filenam=infile.split("/")[-1]
inpath=infile.split("/"+filenam)[0]
obstype=filenam.split(".")[0].lower()
subtyplst=[50100,50500]
print(obstype)
columns=["obsgroup",  "subtype",  "TCWV", "Latitude", "Longitude", "SatView", "SatID" ]

dataset=modulib.obstore_read_file(inpath,obstype,subtyplst=subtyplst)
datframe=dataset["data"]
datset=essio.xar_framegrid(datframe,varlst=columns)

print(datset)

exit()

