#!python
import sys,os
import modulib

if len(sys.argv) > 1:
    inpath = sys.argv[1]
else:
    print("Input data file path is not provided")

filenam=inpath.split("/")[-1]
print(filenam)
obstype=filenam.split(".")[0]
print(obstype)

dataset=modulib.obstore_read_file(inpath,obstype)
print(dataset)

exit()
