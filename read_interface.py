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
p = OptionParser(usage="""usage: %prog [options] <dummy_sloid.vtu> <casename> <in_file.vtu> <solids.dat> <joint_ele.dat> <out_file.vtu> 
Assign material properties according to the surfaces given in "solids.dat" file
Write new_*.ele and new_*.ine files

""")
p.add_option("-v", action="store_true", dest="verbose",  help="Verbose")
p.add_option("-a", action="store_true", dest="ascii",  help="ASCII, instead of binary, .vtu output")
p.add_option("-m", action="store", dest="multikey", type="str", help="Multikey of the mesh partition.  If given, then <out_vtu> will contain cosflow_multikey")
(opts, args) = p.parse_args()


# Get the arguments
if len(args) != 2:
  print "Program must be called with 5 arguments"
  p.print_help()
  sys.exit(1)
(in_file, out_file) = args

r = vtk.vtkXMLUnstructuredGridReader()
r.SetFileName(in_file)
r.Update()

mesh = r.GetOutput()

array = mesh.GetCellData().GetArray("cfelement_attached")
print mesh.GetCellData()

no_of_tuples = mesh.GetCells().GetData().GetNumberOfTuples()
#for i in range(no_of_tuples):
#  print mesh.GetCells().GetData().GetValue(i)


print "DONE"
