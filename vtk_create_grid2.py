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
p = OptionParser(usage="""usage: %prog [options] <vtr_out>

For example:
%prog 

""")
p.add_option("-v", action="store_true", dest="verbose",  help="Verbose")
p.add_option("-a", action="store_true", dest="ascii",  help="ASCII, instead of binary, .vtr output")
(opts, args) = p.parse_args()

# get the com filename
if len(args) != 3:
   p.print_help()
   sys.exit(1)
(xd, yd, vtr_out) = args

#vtr_out = "out.vtr"
#if opts.verbose: print "Reading", bot_file
#bot = vtk.vtkXMLPolyDataReader()
#bot.SetFileName(bot_file)
#bot.Update()
#bot = bot.GetOutput()

xdiv = int(xd)
ydiv = int(yd)
grid = vtk.vtkRectilinearGrid()
grid.SetDimensions(xdiv+1,ydiv+1,1)

#xmin = 312281.0
#xmax = 319638.0

#ymin = 6402884.0
#ymax = 6413777.0

xmin = 313000.0
xmax = 318000.0

ymin = 6403884.0
ymax = 6410777.0


dx = 0.0
dy = 0.0
dx = (xmax-xmin)/xdiv
dy = (ymax-ymin)/ydiv

xArray = vtk.vtkDoubleArray()
for i in range(0,xdiv+1):
  xArray.InsertNextValue(xmin + dx*i)

yArray = vtk.vtkDoubleArray()
for i in range(0,ydiv+1):
  yArray.InsertNextValue(ymin + dy*i)
  
zArray = vtk.vtkDoubleArray()
zArray.InsertNextValue(400.0)

grid.SetXCoordinates(xArray)
grid.SetYCoordinates(yArray)
grid.SetZCoordinates(zArray)

#print grid.GetNumberOfPoints()
#to_write.GetPointData().RemoveArray("vtkValidPointMask")


#if opts.verbose: print "Adding provenance"
#command_used = vtk.vtkStringArray()
#command_used.SetName("provenance")
#command_used.InsertNextValue(" ".join(sys.argv))
#to_write.GetFieldData().AddArray(command_used)
#bounds = orig.GetBounds()
#to_write = vtk.vtkElevationFilter()
#elev.SetHighPoint(0, 0, bounds[5])
#elev.SetLowPoint(0, 0, bounds[4])
#elev.SetScalarRange(bounds[4], bounds[5])
#elev.SetInputConnection(orig.GetProducerPort())
#to_write = to_write.GetOutput()
#rectilinearGridToTetrahedra = vtk.vtkRectilinearGridToTetrahedra()
#rectilinearGridToTetrahedra.SetInputData(grid)
#rectilinearGridToTetrahedra.Update()

if opts.verbose: print "Writing", vtr_out
writer = vtk.vtkXMLRectilinearGridWriter()
writer.SetFileName(vtr_out)
#writer.SetInputData(grid)
if opts.ascii:
   writer.SetDataModeToAscii()
else:
   writer.SetDataModeToBinary()
writer.SetInputConnection(grid.GetProducerPort())
writer.Write()

sys.exit(0)
