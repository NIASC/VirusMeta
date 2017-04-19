#! /usr/bin/env python
import os
import sys

path = sys.argv[1]
listing = os.listdir(path)
xml_files = 0
success_files = 0
for infile in listing:
   if  infile[-8:] ==".segment":
       xml_files += 1
   if  infile[-11:] ==".xmlSuccess":
       success_files += 1

print "success_files: ",success_files
print "xml_files:",xml_files

print "progress ", (float(success_files)/float(xml_files))*100,"% is done"



