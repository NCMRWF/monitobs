#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 12 15:45:39 2021

@author: gibies
"""
#  Some of the basic elements of this program was derived from 
#  automatically generated code with bufr_dump -Dpython
#  Using ecCodes version: 2.16.0

from __future__ import print_function
import traceback
import sys,os
#from eccodes import *
import pandas
import numpy
import collections
CURR_PATH=os.path.dirname(os.path.abspath(__file__))
PKGHOME=os.path.dirname(CURR_PATH)
OBSLIB=os.environ.get('OBSLIB',PKGHOME+"/pylib")
sys.path.append(OBSLIB)
OBSDIC=os.environ.get('OBSDIC',PKGHOME+"/pydic")
sys.path.append(OBSDIC)
OBSNML=os.environ.get('OBSNML',PKGHOME+"/nml")
sys.path.append(OBSNML)
ECBUFRNML=OBSNML+"/ecbufr_fieldname.nml"
AAPP_NML=OBSNML+"/aapp_fieldname.nml"
import obslib

def read_elist(ibufr):
    iterid = codes_keys_iterator_new(ibufr)
    while iterid :
    	iterid = codes_keys_iterator_next(iterid)
	#print("ID:"+str(iterid))
        if not iterid:
        	break
   	key_names = codes_keys_iterator_get_name(iterid)
    	#print("Name:"+key_names)
	#val = codes_get(ibufr, key_names)
	#print("Value:"+str(val))
	if key_names is 'md5Data':
		break
    codes_keys_iterator_delete(ibufr)

def get_dtype(elename,eleindxmaptbl=ECBUFRNML):
	elename=elename.split('#')[-1]
	#print(elename)
	dtype=pandas.read_table(eleindxmaptbl).query("fieldname == @elename").datatype.values[0]
	#print(dtype)
	return(dtype)

def rename_field(data):
    fieldlist=list(data.columns[:])
    print(fieldlist)
    return(data)
     
def get_ecbufr_fieldlist(elemlist=None,nmlfile=None,eleindxmaptbl=ECBUFRNML,subtyp=None):
	if elemlist is None : 
		if subtyp is None: 
			elemlist=obslib.get_key_info(nmlfile,key="elemlist")
		else:
			elemlist=obslib.get_key_info(nmlfile,key="elemlist_"+str(subtyp))
	if isinstance(elemlist, str): elemlist=numpy.fromstring(elemlist[1:-1],sep=',',dtype=int)
	#print(elemlist)
	elist=pandas.DataFrame(collections.Counter(elemlist).items(),columns=["elem","elcnt"])
	elist=elist.sort_values(by=['elem'])
	obs_fieldlist=[]
	ecb_fieldlist=[]
	obsindxlist=[]
	levlist=[]
        #print(elist)
	for elem in elist.elem:
		#print(elem)
		fieldname=pandas.read_table(eleindxmaptbl).query("indx == @elem").fieldname.values[0]
		elename=pandas.read_table(eleindxmaptbl).query("indx == @elem").elename.values[0]
		count=elist.query("elem == @elem").elcnt.values[0]
		#print(count)
		if count == 1 :
		    if fieldname is not numpy.nan:
			#print(fieldname)
			if elem > 261 :
				ecb_fieldlist=ecb_fieldlist+["#"+str(count)+"#"+fieldname]
			else :
				ecb_fieldlist=ecb_fieldlist+[fieldname]
			obs_fieldlist=obs_fieldlist+[elename]
			obsindxlist=obsindxlist+[elem]
			levlist=levlist+[1]
		else :
			chnlmap=obslib.get_key_info(nmlfile,key=fieldname)
			chnlmap=numpy.fromstring(chnlmap[1:-1],sep=',',dtype=int)
			#print(chnlmap)
			for i,indx in enumerate(chnlmap,start=1):
			    if indx > 0:
				#print(i)
				if fieldname is not numpy.nan: 
					ecb_fieldlist=ecb_fieldlist+["#"+str(indx)+"#"+fieldname]
					obs_fieldlist=obs_fieldlist+[elename+"_"+str(i)]
					obsindxlist=obsindxlist+[elem]
					levlist=levlist+[i]
	fieldlist=pandas.DataFrame({"obsindx":obsindxlist,"chnlev":levlist,"elename":obs_fieldlist,"fieldname":ecb_fieldlist})
	fieldlist=fieldlist.set_index(["obsindx","chnlev"]).sort_index(ascending=True)
	#print(fieldlist)
	return(fieldlist)


def read_element(ibufr,field,count,data,eleindxmaptbl=ECBUFRNML):
	elename=field.elename
	fieldname=field.fieldname
	dtype=get_dtype(fieldname,eleindxmaptbl)
	#print(elename,fieldname)
	if dtype in [ "iVal", "rVal" ] : 
		data1=codes_get(ibufr,fieldname)
		data1=[data1]*count
	else :
		data1=codes_get_array(ibufr,fieldname)
	if len(data1) == 1:
		data1=[data1[0]]*count
	#print(len(data))
	#print(len(data1))
	data[elename]=data1
	return(data)

def message_read(ibufr,fieldlist=None,eleindxmaptbl=ECBUFRNML):
    codes_set(ibufr, 'unpack', 1)
    count=codes_get(ibufr, 'numberOfSubsets')
    data=pandas.DataFrame(index=range(1,(count+1),1))
    #print("Message "+str(ibufr)+" contain "+str(count)+" observations")
    #read_elist(ibufr)
    if fieldlist is None:
	fieldlist = ['year', 'month', 'latitude', 'longitude', ]
   
    for indx,field in fieldlist.iterrows() :
	data=read_element(ibufr,field,count,data,eleindxmaptbl=eleindxmaptbl)
    return(data)

def bufr_decode(input_file,nmlfile,eleindxmaptbl=ECBUFRNML,elemlist=None,subtype=None):
    fieldlist=get_ecbufr_fieldlist(nmlfile=nmlfile,eleindxmaptbl=eleindxmaptbl,elemlist=elemlist,subtyp=subtype)
    f = open(input_file, 'rb')
    data=pandas.DataFrame()
    while 1 :
    	ibufr = codes_bufr_new_from_file(f)
    	#print(ibufr)
    	if ibufr is None: break
    	codes_set(ibufr, 'unpack', 1)
    	#print ('Decoding message number '+str(ibufr))
        data1=message_read(ibufr,fieldlist,eleindxmaptbl=eleindxmaptbl)
    	data=data.append(data1, ignore_index=True)
	data=obslib.reset_index(data)
    	codes_release(ibufr)
    if subtype is None: subtype=obslib.get_key_info(nmlfile,"obsubtyp")
    data=data.assign(subtype=[int(subtype)]*len(data))
    f.close()
    return(data)

def bufr_decode_files(inpath,Tnode,slctstr,nmlfile,eleindxmaptbl=None,elemlist=None,subtype=None,keyfieldlst=[],minval=-99999.99,maxval=99999.99):
	if eleindxmaptbl is None : eleindxmaptbl=ECBUFRNML
	searchstring=inpath+"/"+slctstr
	infiles=obslib.globlist(searchstring)
	if len(infiles) == 0: print("File not found: "+searchstring)
	print("Received "+str(len(infiles))+" bufr files")
	data=obslib.DataFrame()
	for infile in infiles[:] :
		print(infile)
		data1=bufr_decode(infile,nmlfile,eleindxmaptbl=eleindxmaptbl,elemlist=elemlist,subtype=subtype)
		for field in keyfieldlst:
			data1=obslib.frame_window_filter(data1,item=field,minval=minval,maxval=maxval)
		print(data1)
		data=data.append(data1.copy())
	data=obslib.reset_index(data)
	#print(data)
	return(data)
