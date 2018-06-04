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
p = OptionParser(usage="""usage: %prog <in_vtp>
Fill large triangles in <in_vtp> with points

Sometimes vtp surfaces with large holes create interpolation difficulties.  This script puts a new point at the centre of each cell which has greater area than twice the median cell area.

Often it is advantageous to run this program a few times to gradually fill in large holes with points.

By default, a Delaunay triangulation of the result is performed.

""")
p.add_option("-v", action="store_true", dest="verbose", help="Verbose output while proceeding")
p.add_option("-a", action="store_true", dest="ascii",  help="ASCII, instead of binary, .vtp output")
p.add_option("-o", action="store", type="str", dest="output", help="Output to this file instead of in_vtp")
p.add_option("--no_delaunay", action="store_true", dest="no_delaunay", help="By default, a Delaunay triangulation will be performed.  This prevents that process happening")
p.add_option("--tol", action="store", type="float", dest="tol", default=0.0, help="During the Delaunay trianulation, points spaced less than d apart will be regarded as a single point.  d = tol*(diagonal length of bounding box of points).   (Default=%default)")

(opts, args) = p.parse_args()
# check that arguments are correct
if len(args) != 1:
   p.print_help()
   sys.exit(1)
(in_vtp) = args[0]


if not in_vtp.endswith(".vtp"):
   print "This program only works with vtp"
   sys.exit(2)

# read the input
if opts.verbose: print "Reading", in_vtp
r = vtk.vtkXMLPolyDataReader()
r.SetFileName(in_vtp)
r.Update()
orig = r.GetOutput()




# find the area of the cells
if opts.verbose: print "Finding area of all", orig.GetNumberOfCells(), "cells"
areas = []
for cellid in range(orig.GetNumberOfCells()):
   try:
      areas.append((orig.GetCell(cellid).ComputeArea(), cellid))
   except:
      print "Cannot find the area of cell", cellid

areas.sort()
med_id = len(areas)/2
median_area = areas[med_id][0]
too_big = [area for area in areas if area[0]>=2*median_area]
if opts.verbose: print "Min area", areas[0][0], " Max area", areas[-1][0], " Median area", areas[med_id][0]
if len(too_big) == 0:
   if opts.verbose: print "No cells need filling"
   sys.exit(0)
else:
   if opts.verbose: print len(too_big), "cells need filling"



# find all the original points
if opts.verbose: print "Getting all original points"
xyz_data = []
for ptid in range(orig.GetNumberOfPoints()):
   xyz_data.append(orig.GetPoint(ptid))


if opts.verbose: print "Putting point at centroid of big cells"
xyz = []*3
for cell in too_big:
   # note, it must be a triangle
   # shove a point at the centroid
   # first, get the points
   the_cell = orig.GetCell(cell[1])
   xyz = [orig.GetPoint(the_cell.GetPointId(ver)) for ver in range(3)]
   centroid = [(xyz[0][coor] + xyz[1][coor] + xyz[2][coor])/3 for coor in range(3)]
   xyz_data.append(centroid)


if opts.verbose: print "Adding new points to new polydata"
pts = vtk.vtkPoints()
pts.SetNumberOfPoints(len(xyz_data))
for ptid in range(len(xyz_data)):
   pts.InsertPoint(ptid, xyz_data[ptid])


vtp_with_points = vtk.vtkPolyData()
vtp_with_points.SetPoints(pts)
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



if opts.verbose: print "Adding provenance"
elev.Update()
elev = elev.GetOutput()
if orig.GetFieldData().HasArray("provenance"):
   command_used = orig.GetFieldData().GetAbstractArray("provenance")
else:
   command_used = vtk.vtkStringArray()
   command_used.SetName("provenance")
command_used.InsertNextValue(" ".join(sys.argv))
elev.GetFieldData().AddArray(command_used)


out_vtp = in_vtp
if opts.output: out_vtp = opts.output
if opts.verbose: print "Writing", out_vtp
writer = vtk.vtkXMLPolyDataWriter()
writer.SetFileName(out_vtp)
if opts.ascii:
   writer.SetDataModeToAscii()
else:
   writer.SetDataModeToBinary()
writer.SetInputConnection(elev.GetProducerPort())
writer.Write()


sys.exit(0)
