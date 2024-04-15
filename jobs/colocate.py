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
print(obstype)
columns=["obsgroup",  "subtype",  "Year" "TCWV"]

dataset=modulib.obstore_read_file(inpath,obstype)
datframe=dataset["data"]
datset=essio.xar_framegrid(datframe)

print(datset)

exit()
