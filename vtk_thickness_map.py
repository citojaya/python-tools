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
p = OptionParser(usage="""usage: %prog [options] <vtu_in> <layer_name> <vtp_out>
Forms a thickness map of elements of name layer_name in vtu_in: writes output as a vtp


""")
p.add_option("-v", action="store_true",        dest="verbose",  help="Verbose")
p.add_option("-a", action="store_true", dest="ascii",  help="ASCII, instead of binary, .vtp output")
(opts, args) = p.parse_args()


if len(args) != 3:
   p.print_help()
   sys.exit(1)
(vtu_in, layer_name, vtp_out) = args


# read vtu_in
if opts.verbose: print "Reading", vtu_in
r = vtk.vtkXMLUnstructuredGridReader()
r.SetFileName(vtu_in)
r.Update()
r = r.GetOutput()



# extract material_names from those stored in vtu_in
mat_names = r.GetFieldData().GetAbstractArray("material_names")
if not mat_names:
   print "Cannot find material_names fielddata in", vtu_in
   sys.exit(1)
for id in range(mat_names.GetNumberOfValues()):
   if mat_names.GetValue(id) == layer_name:
      break
else:
   print "Cannot find", layer_name, "listed in the material_names fielddata"
   sys.exit(2)
   


# use threshold to get the eles of interest
if not r.GetCellData().HasArray("material"):
   print "Cannot find celldata with name 'material' in", vtu_in
   sys.exit(3)
threshold = vtk.vtkThreshold()
threshold.SetInput(r)
threshold.SetInputArrayToProcess(0, 0, 0, 1, "material")
threshold.ThresholdBetween(id-0.5, id+0.5)
threshold.Update()
threshold = threshold.GetOutput()



# find min and max
if opts.verbose: print "Finding the min and max values of all (x,y) points"
min_max_z = {}
for ptid in range(threshold.GetNumberOfPoints()):
   xyz = threshold.GetPoint(ptid)
   if (xyz[0],xyz[1]) not in min_max_z:
      # have not yet encountered this (x,y) posn
      min_max_z[(xyz[0], xyz[1])] = (xyz[2], xyz[2])
   else:
      # have encountered this (x,y) posn
      (min_so_far, max_so_far) = min_max_z[(xyz[0], xyz[1])]
      min_max_z[(xyz[0], xyz[1])] = (min(min_so_far, xyz[2]), max(max_so_far, xyz[2]))



# make a vtk.vtkPoints array
if opts.verbose: print "Putting nodes into vtp", vtp_out
pts = vtk.vtkPoints()
pts.SetNumberOfPoints(len(min_max_z))
ptid = 0
for k, v in min_max_z.items():
   pts.InsertPoint(ptid, k[0], k[1], v[1]-v[0])
   ptid += 1
vtp = vtk.vtkPolyData()
vtp.SetPoints(pts)




# do Delaunay and put elevation
if opts.verbose: print "Performing Delaunay and elevation"
delaunay = vtk.vtkDelaunay2D()
delaunay.SetInput(vtp)
delaunay.Update()
no_orphaned_points = vtk.vtkCleanPolyData()
no_orphaned_points.SetInput(delaunay.GetOutput())
no_orphaned_points.Update()
bounds = no_orphaned_points.GetOutput().GetBounds()
elev = vtk.vtkElevationFilter()
elev.SetHighPoint(0, 0, bounds[5])
elev.SetLowPoint(0, 0, bounds[4])
elev.SetScalarRange(bounds[4], bounds[5])
elev.SetInput(no_orphaned_points.GetOutput())
elev = elev.GetOutput()


if r.GetFieldData().HasArray("provenance"):
   command_used = r.GetFieldData().GetAbstractArray("provenance")
else:
   command_used = vtk.vtkStringArray()
   command_used.SetName("provenance")
command_used.InsertNextValue(" ".join(sys.argv))
elev.GetFieldData().AddArray(command_used)




writer = vtk.vtkXMLPolyDataWriter()
if opts.verbose: print "Writing", vtp_out
writer.SetFileName(vtp_out)
if opts.ascii:
   writer.SetDataModeToAscii()
else:
   writer.SetDataModeToBinary()
writer.SetInputConnection(elev.GetProducerPort())
writer.Write()


sys.exit(0)
