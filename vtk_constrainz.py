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
p = OptionParser(usage="""usage: %prog [options] <vtp>
Constrains the z values and elevation values of <vtp> to be less and/or greater than specified amounts.

Overwrites <vtp> unless the -o option is given.
""")
p.add_option("-v", action="store_true",        dest="verbose",  help="Verbose")
p.add_option("-g", action="store", type="float", dest="greater_than",  help="z values and elevation values must be greater than or equal to this value")
p.add_option("-l", action="store", type="float", dest="less_than",  help="z values and elevation values must be less than or equal to this value")
p.add_option("-a", action="store_true", dest="ascii",  help="ASCII, instead of binary, .vtp output")
p.add_option("-o", action="store", type="str", dest="output", help="Output to this file instead of back to <vtp>")
(opts, args) = p.parse_args()


if len(args) != 1:
   p.print_help()
   sys.exit(1)
   print "Incorrect number of arguments"
(vtp_in) = args[0]

vtp_out = vtp_in
if opts.output: vtp_out = opts.output

if opts.greater_than != None and opts.less_than != None and opts.greater_than > opts.less_than:
   print "Cannot have greater_than amount larger than less_than amount"
   sys.exit(2)



if opts.verbose: print "Reading", vtp_in
orig = vtk.vtkXMLPolyDataReader()
orig.SetFileName(vtp_in)
orig.Update()
orig = orig.GetOutput()



if opts.verbose: print "Bounding"
bounded_points = orig.GetPoints()
if opts.greater_than != None:
   for ptid in range(bounded_points.GetNumberOfPoints()):
      xyz = bounded_points.GetPoint(ptid)
      if xyz[2] < opts.greater_than:
         bounded_points.InsertPoint(ptid, xyz[0], xyz[1], opts.greater_than)
if opts.less_than != None:
   for ptid in range(bounded_points.GetNumberOfPoints()):
      xyz = bounded_points.GetPoint(ptid)
      if xyz[2] > opts.less_than:
         bounded_points.InsertPoint(ptid, xyz[0], xyz[1], opts.less_than)
orig.SetPoints(bounded_points)




if opts.verbose: print "Putting Elevation PointData into the output"
bounds = orig.GetBounds()
elev = vtk.vtkElevationFilter()
elev.SetHighPoint(0, 0, bounds[5])
elev.SetLowPoint(0, 0, bounds[4])
elev.SetScalarRange(bounds[4], bounds[5])
elev.SetInputConnection(orig.GetProducerPort())
elev = elev.GetOutput()



if opts.verbose: print "Adding provenance"
command_used = orig.GetFieldData().GetAbstractArray("provenance")
if not command_used:
   command_used = vtk.vtkStringArray()
   command_used.SetName("provenance")
command_used.InsertNextValue(" ".join(sys.argv))
elev.GetFieldData().AddArray(command_used)




if opts.verbose: print "Writing", vtp_out
writer = vtk.vtkXMLPolyDataWriter()
writer.SetFileName(vtp_out)
if opts.ascii:
   writer.SetDataModeToAscii()
else:
   writer.SetDataModeToBinary()
writer.SetInputConnection(elev.GetProducerPort())
writer.Write()


sys.exit(0)
