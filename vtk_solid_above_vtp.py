#!/usr/bin/env python

import os
import sys
from optparse import OptionParser, OptionValueError

import numpy as np
import scipy.linalg

import vtk
from vtk.util import numpy_support

# so we can find our ../lib no matter how we are called
findbin = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(findbin + '/../lib')

# parse command line
p = OptionParser(usage="""usage: %prog [options] <in_vtp> <height> <out_vtu>
Builds a solid above a .vtp

The .vtp must contain triangles (eg from a Delaunay triangulation).

For example:
%prog surface.vtp 10.0 solid.vtu
Builds a prismatic solid with surface.vtp as its base and 10 high.

""")
p.add_option("-v", action="store_true",        dest="verbose",  help="Verbose")
p.add_option("-a", action="store_true", dest="ascii",  help="ASCII, instead of binary, .vtu output")
(opts, args) = p.parse_args()

if len(args) != 3:
   p.print_help()
   sys.exit(1)
# get the output file and remove it from list
(in_vtp, height, out_vtu) = args
height = float(height)


# read the input
if opts.verbose: print "Reading", in_vtp
base = vtk.vtkXMLPolyDataReader()
base.SetFileName(in_vtp)
base.Update()
base = base.GetOutput()

num_pts = base.GetNumberOfPoints()
if opts.verbose: print "Found", num_pts, "points: extracting them"
# Extract the point coordinates, and include the ones shifted vertically
pts = vtk.vtkPoints()
pts.SetNumberOfPoints(2*num_pts)
for ptid in range(num_pts):
   xyz = base.GetPoint(ptid)
   xyz_vertical = [xyz[0], xyz[1], xyz[2] + height]
   pts.InsertPoint(ptid, xyz)
   pts.InsertPoint(ptid + num_pts, xyz_vertical)


if opts.verbose: print "Creating solid unstructured grid and adding points"
solid = vtk.vtkUnstructuredGrid()
solid.SetPoints(pts)
solid.Update()


if opts.verbose: print "Inserting", base.GetNumberOfPolys(), "wedges into", out_vtu
for polyid in range(base.GetNumberOfPolys()):
   if base.GetCellType(polyid) != 5:
      print "Cannot form a mesh with non-triangular polygons in the base"
      sys.exit(2)
   vtk_cell = base.GetCell(polyid)
   wedge = vtk.vtkWedge().GetPointIds()
   for pt in range(3):
      poly_pt = vtk_cell.GetPointId(pt)
      wedge.SetId(pt, poly_pt)
      wedge.SetId(pt+3, poly_pt+num_pts)
   solid.InsertNextCell(13, wedge)



if opts.verbose: print "Adding provenance"   
if base.GetFieldData().HasArray("provenance"):
   command_used = base.GetFieldData().GetAbstractArray("provenance")
else:
   command_used = vtk.vtkStringArray()
   command_used.SetName("provenance")
command_used.InsertNextValue(" ".join(sys.argv))
solid.GetFieldData().AddArray(command_used)




if opts.verbose: print "Writing", out_vtu
writer = vtk.vtkXMLUnstructuredGridWriter()
writer.SetFileName(out_vtu)
if opts.ascii:
   writer.SetDataModeToAscii()
else:
   writer.SetDataModeToBinary()
writer.SetInputConnection(solid.GetProducerPort())
writer.Write()

sys.exit(0)
