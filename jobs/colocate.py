#!python
import sys,os
import modulib
import essio
import obslib

if len(sys.argv) > 1:
    infile = sys.argv[1]
else:
    print("Input data file path is not provided")

if len(sys.argv) > 2:
    PDY = sys.argv[2]
    print("pdy",sys.argv[2])
else:
    print("Input date is not provided")

#print(PDY)

filenam=infile.split("/")[-1]
inpath=infile.split("/"+filenam)[0]
obstype=filenam.split(".")[0].lower()
#subtyplst=[50100,50500]
print(obstype)
columns=["obsgroup",  "subtype",  "TCWV", "Latitude", "Longitude", "SatView", "SatID" ]

dataset=modulib.obstore_read_file(inpath,obstype,subtyplst=subtyplst)
datframe=dataset["data"]
datset=essio.xar_framegrid(datframe,varlst=columns)

plotdir="/scratch/meena/validation"
obslib.mkdir(plotdir)
#outfile=plotdir+'/sattcwv_{}.nc'.format(PDY)
#ncfile = essio.nix_write(datset,outfile)
print(datset)

exit()

