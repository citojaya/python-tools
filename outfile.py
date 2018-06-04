#!/usr/bin/env python
#import vtk
#reader = vtk.vtkDataSetReader()
#reader.SetFileName("frame1.vtu")

import os
import sys
import re
from optparse import OptionParser, OptionValueError

# parse command line
p = OptionParser(usage="""usage: %prog [options] <data> <mesh> <output>
Interpolates everything in <data> vtu file onto the new mesh given in
<mesh> vtu file.  Writes out into <output> vtu file.
""")
#p.add_option("-v", action="store_true",        dest="verbose",  help="Verbose")
#p.add_option("-a", action="store_true", dest="ascii",  help="ASCII, instead of binary, VTK output")
(opts, args) = p.parse_args()

# get the com filename
if len(args) < 2:
   print 'check input and output files'
   sys.exit(2)
(inputfn1, inputfn2) = args

#fin1 = open(inputfn1, 'r')
#fin2 = open(inputfn2, 'r')
#fin2 = open(inputfn2, 'w')

string1 ='vtk_mix.py -p "1,-1" '+inputfn1+' '+inputfn2+' xxx.vtp'
print string1
os.system(string1)

string4 ='vtk_extend.py xxx.vtp 10000 ext_xxx.vtp'
print string4
os.system(string4)

string2 ='vtk_constrainz.py -g 2 ext_xxx.vtp -o yyy.vtp'
print string2
os.system(string2)

string3 ='vtk_mix.py -p "1,-1" '+inputfn1+' yyy.vtp '+'new_'+inputfn2
print string3
os.system(string3)

name1 = tuple(inputfn1)
name1 = name1[:-4]

name2 = tuple(inputfn2)
name2 = name2[:-4]

string11=''

for i in range(0, len(name1)):
  string11 = string11+name1[i]

string22=''
for i in range(0, len(name2)):
  string22 = string22+name2[i]

string4 ='vtk_mix.py -p "1,-1" '+inputfn1+' '+'new_'+inputfn2+' thick_'+string11+'_'+string22+'.vtp'
#string4 ='vtk_mix.py -p "1,-1" '+inputfn1+' '+inputfn2+' thick_'+string11+'_'+string22+'.vtp'
print string4
os.system(string4)


#name = tuple(inputfn2)
#name = name[:-4]

#string4 = ''
#for i in range(0, len(inputfn2)-4):
#  string4 = string4+name[i]
#string4 = string4+'_1.vtp'

#print string4

sys.exit(0)







