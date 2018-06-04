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
p = OptionParser(usage="""usage: %prog <name> <id> <val> <vtk_file>
Change field data with <name> and <id> to <val> in <vtk_file>.

Note: name cannot be provenance
Note: FieldData is used to store vital model info

""")
p.add_option("-v", action="store_true",        dest="verbose",  help="Verbose")
p.add_option("-a", action="store_true", dest="ascii",  help="ASCII, instead of binary vtk output")
p.add_option("-o", action="store", type="str", dest="output",  help="Output to this file instead of vtk_file.")
(opts, args) = p.parse_args()

# check that arguments are correct
if len(args) != 4:
   p.print_help()
   print "Incorrect arguments"
   sys.exit(1)

(fd_name, fd_id, fd_val, vtk_file) = args
fd_id = int(fd_id)

if fd_name == "provenance":
   print "You cannot change provenance with this program"
   sys.exit(2)

if fd_id < 0:
   print "VTK does not use negative indices"
   sys.exit(5)

output = vtk_file
if opts.output: output = opts.output

if os.path.splitext(vtk_file)[1] != os.path.splitext(output)[1]:
   print "Input and output vtk types must be identical"
   sys.exit(3)


if vtk_file.endswith(".vtu"):
   r = vtk.vtkXMLUnstructuredGridReader()
   writer = vtk.vtkXMLUnstructuredGridWriter()
elif vtk_file.endswith(".vtp"):
   r = vtk.vtkXMLPolyDataReader()
   writer = vtk.vtkXMLPolyDataWriter()
elif vtk_file.endswith(".vtr"):
   r = vtk.vtkRectilinearGridReader()
   writer = vtk.vtkRectilinearGridWriter()
else:
   print "Not currently configured to read file of type", vtk_file

if opts.verbose: print "Reading", vtk_file
r.SetFileName(vtk_file)
r.Update()
r = r.GetOutput()

fd = r.GetFieldData().GetAbstractArray(fd_name)
if not fd:
   print "Cannot find FieldData with name", fd_name, "in", vtk_file
   sys.exit(4)


if fd_id > fd.GetNumberOfTuples():
   print "You can only choose index between 0 and", fd.GetNumberOfTuples(), "for", fd_name
   sys.exit(6)


if fd.GetDataTypeAsString() == "string":
   pass # nothing required
elif fd.GetDataTypeAsString() == "double":
   fd_val = float(fd_val)
else:
   print "Andy has not yet coded for your data types - ask him to include datatype=", fd.GetDataTypeAsString()
   sys.exit(7)


if opts.verbose: print "Adding required value"
if fd_id == fd.GetNumberOfTuples():
   fd.InsertNextValue(fd_val)
else:
   fd.InsertValue(fd_id, fd_val)
r.GetFieldData().AddArray(fd)


if opts.verbose: print "Adding provenance"
if r.GetFieldData().HasArray("provenance"):
   command_used = r.GetFieldData().GetAbstractArray("provenance")
else:
   command_used = vtk.vtkStringArray()
   command_used.SetName("provenance")
command_used.InsertNextValue(" ".join(sys.argv))
r.GetFieldData().AddArray(command_used)


if opts.verbose: print "Outputting", output
writer.SetFileName(output)
if opts.ascii:
   writer.SetDataModeToAscii()
else:
   writer.SetDataModeToBinary()
writer.SetInputConnection(r.GetProducerPort())
writer.Write()


sys.exit(0)
