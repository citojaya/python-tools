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
p = OptionParser(usage="""usage: %prog [options] <vtp1> <vtp2> <vtpfile_out>
Mixes the z values from the two inputs and writes result, including an elevation pointdata array to output.

The two inputs must have a Delaunay triangulation in them.

For example:
%prog -p "0.2,0.8" bot.vtp top.vtp close_to_top.vtp

""")
p.add_option("-v", action="store_true",        dest="verbose",  help="Verbose")
p.add_option("-a", action="store_true", dest="ascii",  help="ASCII, instead of binary, .vtp output")
p.add_option("-p", "--parts", action="store", type="str", default="0.5,0.5", dest="prts", help="p1,p2: z_out = p1*z1 + p2*z2.  Ie, take p1 parts of infile1 and p2 parts of infile2.  Usually p1+p2=1.0.  Default=%default")
(opts, args) = p.parse_args()

# get the com filename
if len(args) != 3:
   p.print_help()
   sys.exit(1)
(bot_file, top_file, vtp_out) = args

if opts.verbose: print "Reading", bot_file
bot = vtk.vtkXMLPolyDataReader()
bot.SetFileName(bot_file)
bot.Update()
bot = bot.GetOutput()

if opts.verbose: print "Reading", top_file
top = vtk.vtkXMLPolyDataReader()
top.SetFileName(top_file)
top.Update()
top = top.GetOutput()


flattener = vtk.vtkTransform()
flattener.Scale(1.0, 1.0, 0.0)
if opts.verbose: print "Flattening", bot_file, "to remove z values of points"
flat_bot = vtk.vtkTransformFilter()
flat_bot.SetInput(bot)
flat_bot.SetTransform(flattener)
flat_bot.Update()
flat_bot = flat_bot.GetOutput()
if opts.verbose: print "Flattening", top_file, "to remove z values of points"
flat_top = vtk.vtkTransformFilter()
flat_top.SetInput(top)
flat_top.SetTransform(flattener)
flat_top.Update()
flat_top = flat_top.GetOutput()



# use the "top" elevation values and put them to "bot"
if opts.verbose: print "Putting elevation values from", top_file, "to", bot_file
bot_mesh_with_top_elev = vtk.vtkProbeFilter()
bot_mesh_with_top_elev.SetInput(flat_bot)
bot_mesh_with_top_elev.SetSource(flat_top)
bot_mesh_with_top_elev.Update()
bot_mesh_with_top_elev = bot_mesh_with_top_elev.GetOutput()



# Form the mix
if opts.verbose: print "Forming the mix", opts.prts
(p_bot, p_top) = map(float, opts.prts.split(","))
bot_elev_array = flat_bot.GetPointData().GetArray("Elevation")
top_elev_array = bot_mesh_with_top_elev.GetPointData().GetArray("Elevation")
final_z = vtk.vtkDoubleArray()
final_z.SetNumberOfValues(flat_bot.GetNumberOfPoints())
final_z.SetName("Elevation")
for ptid in range(flat_bot.GetNumberOfPoints()):
   final_z.SetValue(ptid, p_bot*bot_elev_array.GetValue(ptid) + p_top*top_elev_array.GetValue(ptid))
bot_mesh_with_top_elev.GetPointData().RemoveArray("Elevation")
bot_mesh_with_top_elev.GetPointData().AddArray(final_z)



if opts.verbose: print "Warping the elevation values to z"
bot_with_z = vtk.vtkWarpScalar()
bot_with_z.SetInput(bot_mesh_with_top_elev)
bot_with_z.SetInputArrayToProcess(0,0,0,0,"Elevation")
bot_with_z.Update()
bot_with_z = bot_with_z.GetOutput()


if bot_with_z.GetPointData().GetArray("vtkValidPointMask").GetRange()[0] < 1:
   # there are some invalid points
   if opts.verbose: print "**** Incompatible meshes produced invalid points ****"
   if opts.verbose: print "Removing invalid points"
   only_valid = vtk.vtkThreshold()
   only_valid.SetInput(bot_with_z)
   only_valid.SetInputArrayToProcess(0, 0, 0, 0, "vtkValidPointMask")
   only_valid.ThresholdBetween(0.5, 1.5)
   only_valid.Update()
   only_valid = only_valid.GetOutput()
   only_valid_poly = vtk.vtkGeometryFilter()
   only_valid_poly.SetInput(only_valid)
   only_valid_poly.Update()
   only_valid_poly = only_valid_poly.GetOutput()

   # now must re-mesh
   if opts.verbose: print "Performing Delaunay triangulation"
   delaunay = vtk.vtkDelaunay2D()
   delaunay.SetInput(only_valid_poly)
   delaunay.Update()
   if opts.verbose: print "Getting rid of unused points"
   no_orphaned_points = vtk.vtkCleanPolyData()
   no_orphaned_points.SetInput(delaunay.GetOutput())
   no_orphaned_points.Update()
   to_write = no_orphaned_points.GetOutput()
else:
   to_write = bot_with_z


to_write.GetPointData().RemoveArray("vtkValidPointMask")


if opts.verbose: print "Adding provenance"
command_used = vtk.vtkStringArray()
command_used.SetName("provenance")
command_used.InsertNextValue(" ".join(sys.argv))
to_write.GetFieldData().AddArray(command_used)


if opts.verbose: print "Writing", vtp_out
writer = vtk.vtkXMLPolyDataWriter()
writer.SetFileName(vtp_out)
if opts.ascii:
   writer.SetDataModeToAscii()
else:
   writer.SetDataModeToBinary()
writer.SetInputConnection(to_write.GetProducerPort())
writer.Write()


sys.exit(0)
