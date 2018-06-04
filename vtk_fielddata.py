#!/usr/bin/env python

import os
import sys
import math
from optparse import OptionParser, OptionValueError

import numpy as np
import scipy.linalg

import vtk
from vtk.util import numpy_support

# so we can find our ../lib no matter how we are called
findbin = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(findbin + '/../lib')
import misc.common
from mymath.geom import centroid, bestfit

# parse command line
p = OptionParser(usage="""usage: %prog <vtk_file>
Print field data to stdout in tabular form

""")
(opts, args) = p.parse_args()

# check that arguments are correct
if len(args) != 1:
   p.print_help()
   print "Incorrect arguments"
   sys.exit(1)

(vtk_file) = args[0]


if vtk_file.endswith(".vtu"):
   r = vtk.vtkXMLUnstructuredGridReader()
elif vtk_file.endswith(".vtp"):
   r = vtk.vtkXMLPolyDataReader()
elif vtk_file.endswith(".vtr"):
   r = vtk.vtkRectilinearGridReader()
else:
   print "Not currently configured to read file of type", vtk_file
r.SetFileName(vtk_file)
r.Update()
r = r.GetOutput()

fd = r.GetFieldData()
for fd_id in range(fd.GetNumberOfArrays()):
   the_fd = fd.GetAbstractArray(fd_id)
   sys.stdout.write("FieldDataName=" + the_fd.GetName() + "\n")
   vals = [str(pt) + " " + str(the_fd.GetValue(pt)) for pt in range(the_fd.GetNumberOfTuples())]   # note for future: not sure what happens with vector/tensors here
   sys.stdout.write("FieldDataValues=\n")
   sys.stdout.write("\n".join(vals) + "\n")

sys.exit(0)
