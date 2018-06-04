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

def is_number(s):
   try:
      float(s)
      return True
   except ValueError:
      return False



# parse command line
p = OptionParser(usage="""usage: %prog [options] <output_vtp>
Builds a .vtp file containing polygons from plaintext files


Recursively searches a directory (specified using -d option) for all files ending in .poly.
Extracts the polygons described by these files.

Each file should have the format
x1 y1
x2 y2
...
which describes a polygon, ordered either anticlockwise or clockwise.
The last point is assumed to be joined to the first point.

It is possible to define celldata for the polygons by including lines of the form
#celldata name value
in the poly files.
The tag #celldata is special: all other lines beginning with # are ignored.
The "name" of celldata cannot be poly_name or full_name.

Example polygon file:
# this is a comment line
#comment excav_date 2008/7/26
0 0
300 100
0 900
-300 800
""")
p.add_option("-v", action="store_true",        dest="verbose",  help="Verbose")
p.add_option("-a", action="store_true",        dest="ascii",  help="Write the .vtk vile in ascii format")
p.add_option("-d", action="store", dest="search_dir", default=".", help="Search this directory for .poly files.   (Default=%default)")
p.add_option("--bare", action="store_true", dest="bare", help="Usually polygon names and cell data are included in the .vtp file.  This option turns that feature off")
(opts, args) = p.parse_args()


if len(args) != 1:
   p.print_help()
   print "Incorrect number of arguments"
   sys.exit(1)
(output_vtp) = args[0]



if opts.verbose: print "Finding all files ending with .poly in", opts.search_dir
poly_files = []
for dirname, dirnames, filenames in os.walk(opts.search_dir):
   poly_files += [os.path.join(dirname, filename) for filename in filenames if filename.endswith(".poly")]



vtk_point_number = 0
vtk_cell_number = 0
vtk_all_points = vtk.vtkPoints()
vtk_all_polys = vtk.vtkCellArray()
poly_names = vtk.vtkStringArray()
full_names = vtk.vtkStringArray()
cell_data = {}
cell_data_names = []
vtk_cell_data = {}
for poly_file in poly_files:
   if opts.verbose: print "Reading file", poly_file
   poly_name = os.path.splitext(os.path.basename(poly_file))[0]
   cell_data[vtk_cell_number] = {}
   poly_reader = open(poly_file, "r")
   poly_pts = []
   num_points_in_poly = 0
   for line in poly_reader:
      if not line.strip():
         continue
      if not opts.bare and line.split()[0].lower().startswith("#celldata") and len(line.split()) >= 3:
         (nm, val) = line.split()[1:3]
         cell_data_names.append(nm)
         if is_number(val):
            vtk_cell_data[nm] = vtk.vtkDoubleArray()
            val = float(val)
         else:
            vtk_cell_data[nm] = vtk.vtkStringArray()
         cell_data[vtk_cell_number][nm] = val
         vtk_cell_data[nm].SetName(nm)
         continue
      if line.split()[0].startswith("#"):
         # a comment, or a malformed #celldata line, or opts.bare is turned on
         continue
      num_points_in_poly += 1
      vtk_all_points.InsertNextPoint(float(line.split()[0]), float(line.split()[1]), 0.0)
   poly_reader.close()

   polygon = vtk.vtkPolygon()
   polygon.GetPointIds().SetNumberOfIds(num_points_in_poly)
   for pt in range(num_points_in_poly):
      polygon.GetPointIds().SetId(pt, pt + vtk_point_number)
   vtk_point_number += num_points_in_poly
   vtk_all_polys.InsertNextCell(polygon)
   poly_names.InsertNextValue(poly_name)
   full_names.InsertNextValue(poly_file)
   vtk_cell_number += 1



# build the polygons into polydata
vtk_polygons = vtk.vtkPolyData()
vtk_polygons.SetPoints(vtk_all_points)
vtk_polygons.SetPolys(vtk_all_polys)



# check no name conflict   
cell_data_names = set(cell_data_names)
if "poly_name" in cell_data_names or "full_name" in cell_data_names:
   print "The string poly_name and full_name cannot appear in the #celldata for any poly"
   p.print_help()
   print "The string poly_name and full_name cannot appear in the #celldata for any poly"
   sys.exit(9)



if not opts.bare:
   if opts.verbose: print "Found the following celldata in the polys:"
   for nm in cell_data_names:
      if opts.verbose: print " ", nm
      vtk_cell_data[nm].SetNumberOfValues(len(poly_files))

   # fill the celldata arrays
   for cellid in range(len(poly_files)):
      for (nm, val) in cell_data[cellid].items():
         vtk_cell_data[nm].SetValue(cellid, val)

   # put them into the polydata
   poly_names.SetName("poly_name")
   vtk_polygons.GetCellData().AddArray(poly_names)
   full_names.SetName("full_name")
   vtk_polygons.GetCellData().AddArray(full_names)
   for (nm, ar) in vtk_cell_data.items():
      vtk_polygons.GetCellData().AddArray(ar)


command_used = vtk.vtkStringArray()
command_used.SetName("command_used")
command_used.InsertNextValue(" ".join(sys.argv))
vtk_polygons.GetFieldData().AddArray(command_used)


if opts.verbose: print "Writing", output_vtp
writer = vtk.vtkXMLPolyDataWriter()
writer.SetFileName(output_vtp)
if opts.ascii:
   writer.SetDataModeToAscii()
else:
   writer.SetDataModeToBinary()
writer.SetInputConnection(vtk_polygons.GetProducerPort())
writer.Write()

sys.exit(0)
