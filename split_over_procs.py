#!/usr/bin/env python
import math
import os
import sys
from optparse import OptionParser, OptionValueError

# so we can find our ../lib no matter how we are called
findbin = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(findbin + '/../lib')


# parse command line
p = OptionParser(usage="""usage: %prog [options] <casename> <multikey> <numprocs> <extension> <colnumber>
Split the file <casename>.<extension> over <numprocs> in preparation
for reading by cosflow.

The <colnumber> in <casename>.<extension> is assumed to contain integers
which are either element numbers (in which case %prog should be
called with --type=elemental), or node numbers (%prog should be
called with --type=nodal).  Then <numprocs> files are created,
called <casename><multikey><processor_number>.<extension>

""")
p.add_option("-v", action="store_true", dest="verbose",  help="Verbose")
p.add_option("--type", action="store", dest="split_type", type="str", default="elemental", help="Treat the numbers contained in <colnumber> of <casename>.<extension> as this type.  Default=%default.  Currently you can choose elemental or nodal for this option.")
(opts, args) = p.parse_args()

# get the arguments
if len(args) != 5:
   p.print_help()
   sys.exit(1)
(casename, multikey, numprocs, extension, colnumber) = args

# check the arguments and options are sensible:
numprocs = int(numprocs)
if numprocs <= 1:
   print "The number of processors must be greater than one"
   sys.exit(2)
padit = int(math.log10(numprocs-1)+1)

colnumber = int(colnumber)
if colnumber < 1:
   print "The colnumber must be greater than zero"
   sys.exit(4)

if opts.split_type == "elemental":
   split_ext = ".ele"
   key_posn = 1
elif opts.split_type == "nodal":
   split_ext = ".nod"
   key_posn = 0
else:
   print "The --type flag can only be elemental or nodal, not", opts.split_type
   sys.exit(3)




# Read all ele or nod files and assign their processsor number
if opts.verbose: print "Reading all", split_ext, "files"
cpu_ref_dict = {}
for proc in range(numprocs):
   filename = casename + multikey + str(proc).zfill(padit) + split_ext
   try:
      f = open(filename, 'r')
      for line in f:
         tuple = line.split()
         if not tuple or len(tuple) < key_posn + 1:
            continue
         cpu_ref_dict[tuple[key_posn]] = proc
   except:
      print "Error in processing", filename
      sys.exit(6)


# Read data in original file
if opts.verbose: print "Reading and splitting", casename+"."+extension
try:
   f = open(casename+"."+extension, 'r')
except:
   "Cannot open file", casename+"."+extension
   sys.exit(5)

split_array = [[] for x in range(numprocs)]
for line in f:
   if not line.strip():
      continue
   if len(line.split()) < colnumber:
      print "Not enough columns in line", line
      continue
   num = line.split()[colnumber - 1]
   if num not in cpu_ref_dict:
      print "ERROR", num, "not found in any", split_ext," file.  Line=", line
      continue
   ncpu = cpu_ref_dict[num]
   split_array[ncpu].append(line)
f.close()


if opts.verbose: print "Outputting results"
for proc in range(numprocs):
   filename = casename + multikey + str(proc).zfill(padit) + "." + extension
   try:
      f = open(filename, "w")
      f.write("".join(split_array[proc]))
      f.close()
   except:
      print "ERROR: cannot open or write to file", filename

sys.exit(0)



