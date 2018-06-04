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
p = OptionParser(usage="""usage: %prog [options] <vtk_in> <vtk_out>
Translates vtk_in by given translation.

The types (.vtp, .vtu, etc) of the input and output must be the same.


For example:
%prog -x 1000 -y -1000 -z -123 -r 10 orig.vtp final.vtp
will first translated orig.vtp by (1000, -1000, -123), and then rotate anticlockwise by 10 degrees about the origin.

""")
p.add_option("-v", action="store_true",        dest="verbose",  help="Verbose")
p.add_option("-a", action="store_true", dest="ascii",  help="ASCII, instead of binary, .vtp output")
p.add_option("-x", "--xtrans", action="store", type="float", default=0.0, dest="xtrans", help="X translation applied to data (this is done before any rotation)")
p.add_option("-y", "--ytrans", action="store", type="float", default=0.0, dest="ytrans", help="Y translation applied to data (this is done before any rotation)")
p.add_option("-z", "--ztrans", action="store", type="float", default=0.0, dest="ztrans", help="Z translation applied to data (this is done before any rotation)")
p.add_option("-r", "--rot", action="store", type="float", default=0.0, dest="rot", help="rotation of data in degrees anticlockwise about origin (this is done after any translation)")
(opts, args) = p.parse_args()

# get the com filename
if len(args) != 2:
   p.print_help()
   sys.exit(1)
(in_file, out_file) = args

if opts.verbose: print "Reading", in_file
if in_file.endswith(".vtp") and out_file.endswith(".vtp"):
   indata = vtk.vtkXMLPolyDataReader()
   writer = vtk.vtkXMLPolyDataWriter()
   outdata = vtk.vtkTransformPolyDataFilter()
elif in_file.endswith(".vtu") and out_file.endswith(".vtu"):
   indata = vtk.vtkXMLUnstructuredGridReader()
   writer = vtk.vtkXMLUnstructuredGridWriter()
   outdata = vtk.vtkTransformFilter()
elif in_file.endswith(".vtr") and out_file.endswith(".vtr"):
   indata = vtk.vtkXMLRectilinearGridReader()
   writer = vtk.vtkXMLRectilinearGridWriter()
   outdata = vtk.vtkTransformFilter()
else:
   print "This program has not yet been configured to read", in_file, "and write", out_file
   p.print_help()
   sys.exit(2)


indata.SetFileName(in_file)
indata.Update()
indata = indata.GetOutput()

if opts.verbose: print "Transforming"
trans = vtk.vtkTransform()
trans.PostMultiply()
trans.Translate(opts.xtrans, opts.ytrans, opts.ztrans)
trans.RotateZ(opts.rot) # this will come after the translation

outdata.SetInput(indata)
outdata.SetTransform(trans)
outdata.Update()

if outdata.GetOutput().GetPointData().GetArray("Elevation"):
   if opts.verbose: print "Creating new Elevation array"
   outdata.GetOutput().GetPointData().RemoveArray("Elevation")
   bounds = outdata.GetOutput().GetBounds()
   to_write = vtk.vtkElevationFilter()
   to_write.SetHighPoint(0, 0, bounds[5])
   to_write.SetLowPoint(0, 0, bounds[4])
   to_write.SetScalarRange(bounds[4], bounds[5])
   to_write.SetInput(outdata.GetOutput())
   to_write.Update()
else:
   to_write = outdata


to_write = to_write.GetOutput()   
if indata.GetFieldData().HasArray("provenance"):
   command_used = indata.GetFieldData().GetAbstractArray("provenance")
else:
   command_used = vtk.vtkStringArray()
   command_used.SetName("provenance")
command_used.InsertNextValue(" ".join(sys.argv))
to_write.GetFieldData().AddArray(command_used)


if opts.verbose: print "Writing", out_file
writer.SetFileName(out_file)
if opts.ascii:
   writer.SetDataModeToAscii()
else:
   writer.SetDataModeToBinary()
writer.SetInputConnection(to_write.GetProducerPort())
writer.Write()


sys.exit(0)
