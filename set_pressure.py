#!/tools/python/2.6-x86_64/bin/python
import numpy as np
import os
import random
import math
import sys
import glob
from optparse import OptionParser, OptionValueError



# so we can find our ../lib no matter how we are called
findbin = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(findbin + '/../lib')


# parse command line
p = OptionParser(usage="""usage: %prog [options] <inputfile>


""")
p.add_option("-v", action="store_true",        dest="verbose",  help="Verbose")

(opts, args) = p.parse_args()
# get the com filename
if len(args) != 2:
   p.print_help()
   sys.exit(1)
(casename, p) = args


list_of_files = glob.glob(casename+'*.tar')

#for file_name in list_of_files:
#  phrase = file_name[0:len(file_name)-4]
#  print phrase


for file_name in list_of_files:
  os.system("tar -xvf "+file_name)
  phrase = file_name[0:len(file_name)-10]
  f = open(phrase+"iniporepressure.bdy")

  data = f.readlines()
  f.close()

  for i in range(data.__len__()):
    line = data[i].strip()
    tuple = line.split()
    if(float(tuple[4]) < (-1.0*float(p))):
      tuple[4] = str(-float(p))
      line = ' '.join(tuple)
      data[i] = line+'\n'

  os.system("rm "+phrase+"iniporepressure.bdy")

  file = phrase+"iniporepressure.bdy"
  f = open(file, 'w')
  f.writelines(data)
  f.close()

  os.system("tar -cvf "+file_name+" *.bdy *.history")

  os.system("rm *.bdy *.history")




print "DONE"
