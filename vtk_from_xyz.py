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


def is_number(s):
   try:
      float(s)
      return True
   except ValueError:
      return False


# parse command line
p = OptionParser(usage="""usage: %prog [options] <xyzfile_in> <vtpfile_out>
Converts the input xyz points into vtp format.

xyzfile_in must have columns
x1 y1 z1
x2 y2 z2
....

All comments - that is, lines starting with # - are stored as fielddata in vtpfile_out

xyzfile_in can have additional columns, eg
x1 y1 z1 p1 q1
x2 y2 z2 p2 q2
....
and these data are stored as pointdata arrays in vtpfile_out (with name col4, col5, etc)


By default, a Delaunay triangulation of the input is performed.

For example:
%prog raw.xyz surface.vtp

""")
p.add_option("-v", action="store_true",        dest="verbose",  help="Verbose")
p.add_option("-a", action="store_true", dest="ascii",  help="ASCII, instead of binary, .vtp output")
p.add_option("--no_delaunay", action="store_true", dest="no_delaunay", help="By default, a Delaunay triangulation will be performed.  This prevents that process happening")
p.add_option("--tol", action="store", type="float", dest="tol", default=0.0, help="During the Delaunay trianulation, points spaced less than d apart will be regarded as a single point.  d = tol*(diagonal length of bounding box of points).   (Default=%default)")
(opts, args) = p.parse_args()

# get the com filename
if len(args) != 2:
   p.print_help()
   sys.exit(1)
(xyz_file, vtp_file) = args

provenance = vtk.vtkStringArray()
provenance.SetName("provenance")
provenance.InsertNextValue(" ".join(sys.argv))

comments = vtk.vtkStringArray()
comments.SetName("comments from " + xyz_file)

if opts.verbose: print "Reading", xyz_file
f = open(xyz_file, 'r')
xyz_data = []
extras = []
for line in f:
   if not line.split():
      # empty line
      continue
   if line.split()[0].startswith("#"):
      comments.InsertNextValue(line)
   else:
      xyz_data.append(map(float, line.split()[0:3]))
      extras.append(line.split()[3:])


# determine the number of extra columns and their type
num_extras = len(extras[0])
if num_extras > 0:
   if opts.verbose: print "Processing extra columns"
   extra_pd = []
   extra_type = []
   for extra in extras[0]:
      if is_number(extra):
         extra_pd.append(vtk.vtkDoubleArray())
         extra_type.append("Double")
      else:
         extra_pd.append(vtk.vtkStringArray())
         extra_type.append("String")
   if opts.verbose: 
      print "Found", num_extras, "extra data columns"
      print "Based on first line, i assume their types are:"
      print " ".join(extra_type)
   for extra in range(num_extras):
      extra_pd[extra].SetNumberOfValues(len(xyz_data))
      extra_pd[extra].SetName("col" + str(extra+4))
      if extra_type[extra] == "Double":
         # convert them to doubles now
         for data in extras:
            data[extra] = float(data[extra])


pts = vtk.vtkPoints()
pts.SetNumberOfPoints(len(xyz_data))
for ptid in range(len(xyz_data)):
   pts.InsertPoint(ptid, xyz_data[ptid])
   for extra in range(num_extras):
      extra_pd[extra].SetValue(ptid, extras[ptid][extra])


vtp_with_points = vtk.vtkPolyData()
vtp_with_points.SetPoints(pts)
vtp_with_points.GetFieldData().AddArray(provenance)
vtp_with_points.GetFieldData().AddArray(comments)
for extra in range(num_extras):
   vtp_with_points.GetPointData().AddArray(extra_pd[extra])
vtp_with_points.Update()


# prepare to put the Elevation PointData into the output
bounds = vtp_with_points.GetBounds()
elev = vtk.vtkElevationFilter()
elev.SetHighPoint(0, 0, bounds[5])
elev.SetLowPoint(0, 0, bounds[4])
elev.SetScalarRange(bounds[4], bounds[5])

if not opts.no_delaunay:
   if opts.verbose: print "Performing Delaunay triangulation"
   delaunay = vtk.vtkDelaunay2D()
   delaunay.SetTolerance(opts.tol)
   delaunay.SetInput(vtp_with_points)
   delaunay.Update()
   if opts.verbose: print "Getting rid of unused points"
   no_orphaned_points = vtk.vtkCleanPolyData()
   no_orphaned_points.SetInput(delaunay.GetOutput())
   no_orphaned_points.Update()
   elev.SetInput(no_orphaned_points.GetOutput())
else:
   elev.SetInput(vtp_with_points)



if opts.verbose: print "Writing", vtp_file
writer = vtk.vtkXMLPolyDataWriter()
writer.SetFileName(vtp_file)
if opts.ascii:
   writer.SetDataModeToAscii()
else:
   writer.SetDataModeToBinary()
writer.SetInputConnection(elev.GetOutput().GetProducerPort())
writer.Write()

sys.exit(0)
