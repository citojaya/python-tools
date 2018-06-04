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
p = OptionParser(usage="""usage: %prog [options] <data> <mesh> <output>
Interpolates everything (viz Elevation values) in <data> vtk file onto
the new mesh given in <mesh> vtk file.  Then warps z values Elevation.
Writes out into <output> vtk file.

This kind of interpolation only really makes sense for vtk files describing surfaces.
The algorithm does the following:
 - flattens <data> and <mesh> by removing their z values
 - interpolates PointData (viz Elevation values) in <data> onto <mesh>
   using their (x,y) values only
 - Warps the resulting <data> by the Elevation PointData (if existant)

Another way to describe this process is that the mesh in <data> is
replaced with the mesh in <mesh>.

For example:
%prog elevation.vtp grid.vtp grid_with_elevation.vtp

""")
p.add_option("-v", action="store_true",        dest="verbose",  help="Verbose")
p.add_option("-a", action="store_true", dest="ascii",  help="ASCII, instead of binary, VTK output")
(opts, args) = p.parse_args()


# get arguments
if len(args) < 3:
   p.print_help()
   print " "
   print "Wrong number of arguments"
   sys.exit(1)
(datafn, meshfn, outputfn) = args




# work out readers and writers to use
if datafn.endswith(".vtu"):
   old_mesh_reader = vtk.vtkXMLUnstructuredGridReader()
elif datafn.endswith(".vtp"):
   old_mesh_reader = vtk.vtkXMLPolyDataReader()
else:
   print "Program not yet configured for type of input file", datafn
   sys.exit(2)

if meshfn.endswith(".vtu"):
   new_mesh_reader = vtk.vtkXMLUnstructuredGridReader()
elif meshfn.endswith(".vtp"):
   new_mesh_reader = vtk.vtkXMLPolyDataReader()
elif meshfn.endswith(".vts"):
   new_mesh_reader = vtk.vtkXMLStructuredGridReader()
else:
   print "Program not yet configured for type of input file", meshfn
   sys.exit(2)

if outputfn.endswith(".vtu"):
   writer = vtk.vtkXMLUnstructuredGridWriter()
elif outputfn.endswith(".vtp"):
   writer = vtk.vtkXMLPolyDataWriter()
elif outputfn.endswith(".vts"):
   writer = vtk.vtkXMLStructuredGridWriter()

else:
   print "Program not yet configured for type of output file", outputfn
   sys.exit(2)



# to get rid of z values
flattener = vtk.vtkTransform()
flattener.Scale(1.0, 1.0, 0.0)



# read the input and the source
if opts.verbose: print "Reading", datafn
old_mesh_reader.SetFileName(datafn)
old_mesh_reader.Update()
data = old_mesh_reader.GetOutput()
data_flat = vtk.vtkTransformFilter()
data_flat.SetInput(data)
data_flat.SetTransform(flattener)
data_flat.Update()
data_flat = data_flat.GetOutput()


if opts.verbose: print "Reading", meshfn
new_mesh_reader.SetFileName(meshfn)
new_mesh_reader.Update()
mesh = new_mesh_reader.GetOutput()
mesh_flat = vtk.vtkTransformFilter()
mesh_flat.SetInput(mesh)
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

# output
if opts.verbose: print "Writing", outputfn
writer.SetFileName(outputfn)
if opts.ascii:
   writer.SetDataModeToAscii()
else:
   writer.SetDataModeToBinary()
writer.SetInputConnection(to_output.GetProducerPort())
writer.Write()

