#!/usr/bin/env python

import os
import sys
import subprocess
import re
from optparse import OptionParser, OptionValueError

import numpy as np
import scipy.linalg

import vtk
from vtk.util import numpy_support

# so we can find our ../lib no matter how we are called
findbin = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(findbin + '/../lib')


# parse command line
p = OptionParser(usage="""usage: %prog [options] <vtkfile_in> <xyzfile_out> 
Converts the input vtk file into xyz format

For example:
%prog surface.vtp raw.xyz 

""")
p.add_option("-v", action="store_true",        dest="verbose",  help="Verbose")
p.add_option("--no_pointdata", action="store_true",        dest="no_pointdata",  help="By default, pointdata in vtkfile_in is stored as extra columns in xyzfile_out.  This option turns that off")
p.add_option("--no_header_info", action="store_true",  dest="no_header_info", help="Do not put header info into xyzfile_out")
(opts, args) = p.parse_args()

# get the com filename
if len(args) != 2:
   p.print_help()
   print "Number of arguments incorrect"
   sys.exit(1)
(vtk_file, xyz_file) = args


if opts.verbose: print "Reading", vtk_file
if vtk_file.endswith(".vtu"):
   r = vtk.vtkXMLUnstructuredGridReader()
elif vtk_file.endswith(".vtp"):
   r = vtk.vtkXMLPolyDataReader()
else:
   print "Not currently configured to read file of type", vtk_file
r.SetFileName(vtk_file)
r.Update()
r = r.GetOutput()


f = open(xyz_file, 'w')

if not opts.no_header_info:
   if opts.verbose: print "Writing header info into", xyz_file
   f.write("# This file was generated using the command\n")
   f.write("# " + ' '.join(sys.argv) + "\n")
   f.write("#\n")
   fd = r.GetFieldData()
   if fd:
      f.write("# Field data from " + vtk_file + "\n")
      for fd_id in range(fd.GetNumberOfArrays()):
         the_fd = fd.GetAbstractArray(fd_id)
         f.write("# name=\n" + "#" + the_fd.GetName() + "\n")
         vals = [str(the_fd.GetValue(pt)) for pt in range(the_fd.GetNumberOfTuples())]  # note for future: not sure what happens with vector/tensors here
         f.write("# values=\n#" + "\n#".join(vals) + "\n")
      f.write("# End field data\n")
      f.write("#\n")


if opts.verbose: print "Formatting lines in", xyz_file
if not opts.no_header_info: f.write("# Format for rest of file is\n")
form = "# x y z "
pd = r.GetPointData()
# note for future: not sure what happens with vector/tensors here
if pd and not opts.no_pointdata:
   pd_array = [pd.GetArray(pd_id) for pd_id in range(pd.GetNumberOfArrays())]
   form += " ".join([str(p_data.GetName()) for p_data in pd_array])
else:
   pd_array = []
if not opts.no_header_info: f.write(form + "\n#\n")


if opts.verbose: print "Writing data to", xyz_file
for pt_id in range(r.GetNumberOfPoints()):
   xyz = map(str, r.GetPoint(pt_id))
   xyz += [str(p_data.GetValue(pt_id)) for p_data in pd_array]
   f.write(' '.join(xyz) + "\n")

sys.exit(0)
