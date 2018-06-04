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
p = OptionParser(usage="""usage: %prog [options] <in_file> <out_file> 

Creates a rectilinear grid with given dimensions for a given surface so that PRPEMESH can read. 

For example:
%prog input.vtp grid.vtr output.vtr

""")
p.add_option("-v", action="store_true", dest="verbose",  help="Verbose")
p.add_option("-a", action="store_true", dest="ascii",  help="ASCII, instead of binary, .vtr output")
(opts, args) = p.parse_args()

# get the com filename
if len(args) != 3:
   p.print_help()
   sys.exit(1)
(in_file, grid_file, out_file) = args

# work out readers and writers to use
if in_file.endswith(".vtu"):
   old_mesh_reader = vtk.vtkXMLUnstructuredGridReader()
elif in_file.endswith(".vtp"):
   old_mesh_reader = vtk.vtkXMLPolyDataReader()
else:
   print "Program not yet configured for type of input file", in_file
   sys.exit(2)

# to get rid of z values
flattener = vtk.vtkTransform()
flattener.Scale(1.0, 1.0, 0.0)

# read the input and the source
if opts.verbose: print "Reading", in_file
old_mesh_reader.SetFileName(in_file)
old_mesh_reader.Update()
data = old_mesh_reader.GetOutput()
data_flat = vtk.vtkTransformFilter()
data_flat.SetInput(data)
data_flat.SetTransform(flattener)
data_flat.Update()
data_flat = data_flat.GetOutput()


if opts.verbose: print "Reading mesh"
rect_grid_reader = vtk.vtkXMLRectilinearGridReader()
rect_grid_reader.SetFileName(grid_file)
rect_grid_reader.Update()


grid = vtk.vtkStructuredGrid()
x_coords = rect_grid_reader.GetOutput().GetXCoordinates()
y_coords = rect_grid_reader.GetOutput().GetYCoordinates()
z_coords = rect_grid_reader.GetOutput().GetZCoordinates()

dims = rect_grid_reader.GetOutput().GetDimensions()

xdiv = dims[0]
ydiv = dims[1]

points = vtk.vtkPoints()

for j in range(0,ydiv):
  for i in range(0,xdiv):
     points.InsertNextPoint(x_coords.GetValue(i),y_coords.GetValue(j),0)


grid.SetDimensions(dims)
grid.SetPoints(points)


mesh_flat = vtk.vtkTransformFilter()
mesh_flat.SetInput(grid)
mesh_flat.SetTransform(flattener)
mesh_flat.Update()
mesh_flat = mesh_flat.GetOutput()

# do the work
if opts.verbose: print "Interpolating"
probe = vtk.vtkProbeFilter()
probe.SetInput(mesh_flat)    # the input is the new mesh
probe.SetSource(data_flat)   # the source is the data
probe.Update()
probe = probe.GetOutput()
if opts.verbose: print "Done interpolating"

# warp output using elevation values
if probe.GetPointData().GetArray("Elevation"):
   to_output = vtk.vtkWarpScalar()
   to_output.SetInput(probe)
   to_output.SetInputArrayToProcess(0,0,0,0,"Elevation")
   to_output.Update()
   to_output = to_output.GetOutput()
else:
   to_output = probe

if opts.verbose: print "Adding provenance"
command_used = vtk.vtkStringArray()
command_used.SetName("provenance")
command_used.InsertNextValue(" ".join(sys.argv))
to_output.GetFieldData().AddArray(command_used)

# Converting to vtr PREPMESH format
data = to_output

if opts.verbose: print "Converting to PREPMESH vtr format"
dims = vtk.vtkDoubleArray()
dims = data.GetDimensions()

grid = vtk.vtkRectilinearGrid()
grid.SetDimensions(dims[0],dims[1],dims[2])

xArray = vtk.vtkDoubleArray()
xArray.SetNumberOfComponents(1)
yArray = vtk.vtkDoubleArray()
yArray.SetNumberOfComponents(1)
zArray = vtk.vtkDoubleArray()
zArray.SetNumberOfComponents(1)
elevation_Array = vtk.vtkDoubleArray()
elevation_Array.SetNumberOfComponents(1)

xArray.SetName("X_COORDINATES")
yArray.SetName("Y_COORDINATES")
zArray.SetName("Z_COORDINATES")

counter = 0

points = data.GetPoints()

p = points.GetData()
# Get elevation
for i in range(0,dims[0]*dims[1]*dims[2]):
  elevation_Array.InsertNextValue(p.GetValue(counter+2))
  counter = counter + 3

counter = 0
#X coordinates
for i in range(0,dims[0]):
  xArray.InsertNextValue(p.GetValue(counter))
  counter = counter + 3

counter = 0
#Y coordinates
for i in range(0,dims[1]):
  yArray.InsertNextValue(p.GetValue(counter+1))
  counter = counter + 3*dims[0]

# Set some arbitary value for Z
zArray.InsertNextValue(100.0)

elevation_Array.SetNumberOfTuples(dims[0]*dims[1]*dims[2])
elevation_Array.SetName("z_vals")

grid.SetXCoordinates(xArray)
grid.SetYCoordinates(yArray)
grid.SetZCoordinates(zArray)

grid.GetPointData().AddArray(elevation_Array)

if opts.verbose: print "Writing temp_vtrfile.vtr"
writer = vtk.vtkXMLRectilinearGridWriter()
writer.SetFileName("temp_vtrfile.vtr")
proceed = 0
if opts.ascii:
   proceed = 1
   writer.SetDataModeToAscii()
else:
   writer.SetDataModeToBinary()
writer.SetInputConnection(grid.GetProducerPort())
writer.Write()

startData = 0
#endData = 0

if opts.verbose: print "Converting to PREPMESH format"
# Arrange vtr file in Prepmesh format
if proceed == 1:

  in_file = open("temp_vtrfile.vtr", 'r')
  out_file = open(out_file, 'w')
  while True:
    line = in_file.readline()
    if len(line) == 0:
      break

    tuple = ()
    tuple = line.split()
    if(tuple[0] == "</DataArray>"):
      startData = 0

# Start dataArray rearrangement
    if (startData == 1):
      for i in range(len(tuple)):
        out_file.write("\t"+tuple[i]+"\n")

    else:
      out_file.write(line)

    if(tuple[0] == "<DataArray"):
      startData = 1
  os.system("rm temp_vtrfile.vtr")
  in_file.close()
  out_file.close()

else:
  print "Data is in Binary format, Cannot convert to PREPMESH format"
  sys.exit(0)

sys.exit(0)

