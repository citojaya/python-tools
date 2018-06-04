#!/tools/python/2.6-x86_64/bin/python
import numpy as np
import os
import random
import math
import sys
from optparse import OptionParser, OptionValueError



# so we can find our ../lib no matter how we are called
findbin = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(findbin + '/../lib')


# parse command line
p = OptionParser(usage="""usage: %prog [options] <outputfile>
Writes sample dump file


""")
p.add_option("-v", action="store_true",        dest="verbose",  help="Verbose")
p.add_option("-a", action="store_true", dest="ascii",  help="ASCII, instead of binary, .vtu output")

(opts, args) = p.parse_args()
# get the com filename
if len(args) != 3:
   p.print_help()
   sys.exit(1)
(in_file, min_frame, max_frame) = args

for i in range(int(min_frame), int(max_frame)):
  infile = in_file+str(i)+'.dat'
  outfile = in_file+str(i)+'.vtp'
  print infile
  print outfile
  os.system('vtk_from_xyz.py '+infile+' '+outfile)


print "DONE"
