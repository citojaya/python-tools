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
p = OptionParser(usage="""usage: %prog [options] <particle_smallest_diameter> <total_np> <outputfile>
Writes sample dump file


""")
p.add_option("-v", action="store_true",        dest="verbose",  help="Verbose")
p.add_option("-a", action="store_true", dest="ascii",  help="ASCII, instead of binary, .vtu output")
p.add_option("-s", action="store", dest="input_dump", type="str", help="read sample dump file")

(opts, args) = p.parse_args()
# get the com filename
if len(args) != 3:
   p.print_help()
   sys.exit(1)
(s_bdy,p_bdy,curve) = args

nodes_to_be_changed = {}

f = open(s_bdy, 'r')
s_ini = f.readlines()
f.close()

print s_ini.__len__()

for i in range(s_ini.__len__()):
  line = s_ini[i].strip()
  tuple = line.split()

  if float(tuple[1]) < 0.1:
    tuple[1] = str(0.1)
    line = ' '.join(tuple)
    s_ini[i] = line+'\n'
    nodes_to_be_changed.append(i)
f.close()

f = open("new"+s_bdy, 'w')
f.writelines(s_ini)
f.close()

exit(0)

f = open(curve, 'r')
curve = f.readlines()
f.close()

new_p = 0
for i in range(curve.__len__()):
  new_p = 10.0



f = open(p_bdy, 'r')
p_ini = f.readlines()
f.close()

for i in range(p_ini.__len__()):
  if i in nodes_to_be_changed:
    print i
    line = p_ini[i].strip()
    tuple = line.split()
    tuple[1] = str(new_p)



print "DONE"
