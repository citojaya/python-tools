#!/tools/python/2.6-x86_64/bin/python
import os
import sys

import re
from optparse import OptionParser, OptionValueError

import numpy as np

import glob
import vtk
from vtk.util import numpy_support

# so we can find our ../lib no matter how we are called
findbin = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(findbin + '/../lib')


# parse command line
p = OptionParser(usage="""

""")
p.add_option("-v", action="store_true", dest="verbose",  help="Verbose")
p.add_option("-a", action="store_true", dest="ascii",  help="ASCII, instead of binary, .vtu output")
(opts, args) = p.parse_args()


# Get the arguments
if len(args) != 2:
  print "Program must be called with 2 arguments"
  p.print_help()
  sys.exit(1)
(casename, node) = args



list_of_files = glob.glob(casename+'*.history')

s = ''

for file_name in list_of_files:
  f = open(file_name, 'r')
  data = f.readlines()
  f.close()

  outdata = []

  for i in range(data.__len__()):
    line = data[i].strip()
	tuple = line.split()
    if(int(tuple[3]) != int(node)):
	  outdata.append(line+'\n')
  
  fout = open("new_"+file_name, 'w')
  fout.writelines(outdata)
  fout.close()

print "DONE"


