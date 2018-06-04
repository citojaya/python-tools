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
p = OptionParser(usage="""usage: %prog <in_vtp> <radius> <out_vtp>
Extend (x,y,z) data found in in_vtp to given radius, by forming a planar approximation to the surface.

The in_vtp must contain a Delaunay triangulation.

By default, a Delaunay triangulation of the result is performed.

""")
p.add_option("-v", action="store_true", dest="verbose", help="Verbose output while proceeding")
p.add_option("-a", action="store_true", dest="ascii",  help="ASCII, instead of binary, .vtp output")
p.add_option("-r", action="store", dest="res", type=int, default=10, help="Radial resolution  (default=%default)")
p.add_option("--no_delaunay", action="store_true", dest="no_delaunay", help="By default, a Delaunay triangulation will be performed.  This prevents that process happening")
p.add_option("--tol", action="store", type="float", dest="tol", default=0.0, help="During the Delaunay trianulation, points spaced less than d apart will be regarded as a single point.  d = tol*(diagonal length of bounding box of points).   (Default=%default)")

(opts, args) = p.parse_args()
# check that arguments are correct
if len(args) != 3:
   p.print_help()
   sys.exit(1)
(in_vtp, outer_radius, out_vtp) = args
outer_radius = float(outer_radius)

# read the input
if opts.verbose: print "Reading", in_vtp
r = vtk.vtkXMLPolyDataReader()
r.SetFileName(in_vtp)
r.Update()
orig = r.GetOutput()


# find all the points
pts = []
for ptid in range(orig.GetNumberOfPoints()):
   pts.append(orig.GetPoint(ptid))

# form the plane
if opts.verbose:
   print "Forming the planar approximation"
cd = centroid(pts)
vxz = bestfit([[pt[0],pt[2]] for pt in pts])
if vxz[0] == 0.0:
   print "Error: surface has infinite slope when projected onto x-z plane"
a = vxz[1]/vxz[0]
vyz = bestfit([[pt[1],pt[2]] for pt in pts])
if vyz[0] == 0.0:
   print "Error: surface has infinite slope when projected onto y-z plane"
b = vyz[1]/vyz[0]

if opts.verbose: print "Slope params =", a, b


# generate points on a annulus
if opts.verbose: print "Generating points on an annulus around data"
dist_between_pts = outer_radius/opts.res
circumferential_resolution = int(2*3.1415*outer_radius/dist_between_pts)
ann = vtk.vtkDiskSource()
ann.SetInnerRadius(0.0)
ann.SetOuterRadius(outer_radius)
ann.SetRadialResolution(opts.res)
ann.SetCircumferentialResolution(circumferential_resolution)
ann.Update()

# transform this annulus to the data
(xmin, xmax, ymin, ymax, zmin, zmax) = orig.GetBounds()
(x_centre, y_centre) = (0.5*(xmax+xmin), 0.5*(ymax+ymin))
trans = vtk.vtkTransform()
trans.Translate(x_centre, y_centre, 0)
newpts = vtk.vtkTransformFilter()
newpts.SetInput(ann.GetOutput())
newpts.SetTransform(trans)
newpts.Update()
newpts = newpts.GetOutput()



# only include points which are not covered by original data
flattener = vtk.vtkTransform()
flattener.Scale(1.0, 1.0, 0.0)
flat = vtk.vtkTransformFilter()
flat.SetInput(orig)
flat.SetTransform(flattener)
flat.Update()
loc = vtk.vtkCellLocator()
loc.SetDataSet(flat.GetOutput())
loc.BuildLocator()
extension_pts = []
for ptid in range(newpts.GetNumberOfPoints()):
   xyz = newpts.GetPoint(ptid)
   if loc.FindCell(xyz) < 0:
      extension_pts.append(xyz)

if len(extension_pts) == 0:
   print "Your original data entirely covers the area to be extrapolated to"
   print "No extrapolation will be done"
   print "Perhaps your radius is set too small"
   sys.exit(2)

# do the extrapolation
if opts.verbose: print "Extrapolating z values to those points by limiting to a plane"
loc = vtk.vtkPointLocator()
loc.SetDataSet(flat.GetOutput())
loc.BuildLocator()
vtk_pts = vtk.vtkPoints()
num_origs = orig.GetNumberOfPoints()
num_exts = len(extension_pts)
vtk_pts.SetNumberOfPoints(num_origs + num_exts)
for ptid in range(num_origs):
   vtk_pts.InsertPoint(ptid, pts[ptid])
for ptid in range(num_exts):
   xyz = extension_pts[ptid]
   z_plane = cd[2] + a*(xyz[0]-cd[0]) + b*(xyz[1]-cd[1])
   closest_ptid = loc.FindClosestPoint(xyz)
   if closest_ptid > 0:
      (x_closest, y_closest, z_closest) = pts[closest_ptid]
   else:
      # can't find a closest point for some reason
      (x_closest, y_closest, z_closest) = (xyz[0], xyz[1], z_plane)
   dist_from_closest = ((xyz[0]-x_closest)**2 + (xyz[1]-y_closest)**2)**0.5
   dist_from_edge = outer_radius - ((xyz[0]-x_centre)**2 + (xyz[1]-y_centre)**2)**0.5
   z_to_use = (dist_from_closest*z_plane + dist_from_edge*z_closest)/(dist_from_closest + dist_from_edge)
   vtk_pts.InsertPoint(ptid + num_origs, xyz[0], xyz[1], z_to_use)


# insert points into polydata and do delaunay
if opts.verbose: print "Creating new polydata"
ext = vtk.vtkPolyData()
ext.SetPoints(vtk_pts)
ext.Update()

# prepare to put the Elevation PointData into the output
bounds = vtk_pts.GetBounds()
elev = vtk.vtkElevationFilter()
elev.SetHighPoint(0, 0, bounds[5])
elev.SetLowPoint(0, 0, bounds[4])
elev.SetScalarRange(bounds[4], bounds[5])

if not opts.no_delaunay:
   if opts.verbose: print "Performing Delaunay triangulation"
   delaunay = vtk.vtkDelaunay2D()
   delaunay.SetTolerance(opts.tol)
   delaunay.SetInput(ext)
   delaunay.Update()
   if opts.verbose: print "Getting rid of unused points"
   no_orphaned_points = vtk.vtkCleanPolyData()
   no_orphaned_points.SetInput(delaunay.GetOutput())
   no_orphaned_points.Update()
   elev.SetInput(no_orphaned_points.GetOutput())
else:
   elev.SetInput(ext)


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
