#!/usr/bin/env python

import os
import sys
from optparse import OptionParser, OptionValueError
import subprocess
import tempfile

import numpy as np
import scipy.linalg
import glob
import vtk
import math
from vtk.util import numpy_support


# so we can find our ../lib no matter how we are called
findbin = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(findbin + '/../lib')


def flux_area(verbose, poly_data):
   """ Calculate cosflow flux area for both eles and nodes in a surface

   @type  verbose:    bool
   @param verbose:    whether to print progress to screen
   @type  poly_data:  vtkPolyData
   @param poly_data:  poly data describing the surface
   @rtype:            tuple of vtkFloatArray
   @return:           (nodal flux area, elemental flux area)
   """

   num_pts = poly_data.GetNumberOfPoints()
   num_cells = poly_data.GetNumberOfCells()

   #if verbose: print " Flux Area: Recording cell ids"
   cell_ids = vtk.vtkIdFilter()
   cell_ids.SetInputConnection(poly_data.GetProducerPort())
   cell_ids.CellIdsOn()
   cell_ids.Update()
   cell_ids = cell_ids.GetOutput() # vtkPolyData


   #if verbose: print " Flux Area: Triangulating"
   tri_os = vtk.vtkTriangleFilter()
   tri_os.SetInputConnection(cell_ids.GetProducerPort())
   tri_os.Update()
   tri_os = tri_os.GetOutput() # vtkPolyData
   orig_cell_ids = tri_os.GetCellData().GetScalars("vtkIdFilter_Ids")

   #if verbose: print " Flux Area: Calculating areas"
   cell_areas = [0.0 for cell in range(num_cells)]
   nodal_areas = [0.0 for nod in range(num_pts)]
   for cellid in range(tri_os.GetNumberOfCells()):
      the_cell = tri_os.GetCell(cellid)
      #vtk_points = the_cell.GetPoints()
      orig_cell_id = orig_cell_ids.GetValue(cellid)
      tri_area = the_cell.ComputeArea()

      xyz1 = the_cell.GetPoints().GetPoint(0)
      xyz2 = the_cell.GetPoints().GetPoint(1)
      xyz3 = the_cell.GetPoints().GetPoint(2)
    
      p1p2 = math.sqrt(math.pow((xyz1[0]-xyz2[0]),2)+math.pow((xyz1[1]-xyz2[1]),2)+math.pow((xyz1[2]-xyz2[2]),2))
      p1p3 = math.sqrt(math.pow((xyz1[0]-xyz3[0]),2)+math.pow((xyz1[1]-xyz3[1]),2)+math.pow((xyz1[2]-xyz3[2]),2))
      p2p3 = math.sqrt(math.pow((xyz2[0]-xyz3[0]),2)+math.pow((xyz2[1]-xyz3[1]),2)+math.pow((xyz2[2]-xyz3[2]),2))
      
      theta1 = math.acos((p1p3*p1p3 + p1p2*p1p2 - p2p3*p2p3)/(2.0*p1p3*p1p2))
      theta2 = math.acos((p1p2*p1p2 + p2p3*p2p3 - p1p3*p1p3)/(2.0*p1p2*p2p3))
      theta3 = math.acos((p1p3*p1p3 + p2p3*p2p3 - p1p2*p1p2)/(2.0*p1p3*p2p3))

      ptid = the_cell.GetPointId(0)
      nodal_areas[ptid] += tri_area*theta1/math.pi
      ptid = the_cell.GetPointId(1)
      nodal_areas[ptid] += tri_area*theta2/math.pi
      ptid = the_cell.GetPointId(2)
      nodal_areas[ptid] += tri_area*theta3/math.pi
      cell_areas[orig_cell_id] += tri_area

      #for pt in range(3):
      #   ptid = the_cell.GetPointId(pt)
         # equally proportion triangle areas to the points
      #   nodal_areas[ptid] += tri_area/3.0
      #cell_areas[orig_cell_id] += tri_area

   #if verbose: print " Flux Area: Building cell and point area vtk array"
   na = vtk.vtkFloatArray()
   na.SetNumberOfValues(num_pts)
   na.SetName("cosflow_nodal_flux_area")
   for ptid in range(num_pts):
      na.SetValue(ptid, nodal_areas[ptid])
   ca = vtk.vtkFloatArray()
   ca.SetNumberOfValues(num_cells)
   ca.SetName("cosflow_ele_flux_area")
   for cellid in range(num_cells):
      ca.SetValue(cellid, cell_areas[cellid])

   return (na, ca)

def generate_vtp():
  pts = vtk.vtkPoints()
  pts.SetNumberOfPoints(len(xyz_data))
  for ptid in range(len(xyz_data)):
    pts.InsertPoint(ptid, xyz_data[ptid])

  vtp_with_points = vtk.vtkPolyData()
  vtp_with_points.SetPoints(pts)
  vtp_with_points.Update()

  if opts.verbose: print "Performing Delaunay triangulation"
  delaunay = vtk.vtkDelaunay2D()
  #delaunay.SetTolerance(opts.tol)
  delaunay.SetInput(vtp_with_points)
  delaunay.Update()
  if opts.verbose: print "Getting rid of unused points"
  no_orphaned_points = vtk.vtkCleanPolyData()
  no_orphaned_points.SetInput(delaunay.GetOutput())
  no_orphaned_points.Update()

  return (no_orphaned_points.GetOutput())

def write_flux_area_and_seepage():
  to_write_area = []
  to_write_seepage = []

  for i in range(len(required_nodes_ptids)):
    line = 'SET FLUX_AREA '+str(cf_nodes.GetValue(required_nodes_ptids[i]))+' '+str(area[0].GetValue(i))
    to_write_area.append(line+'\n')

    line = 'SET FLUIDCONSTRAINT '+str(cf_nodes.GetValue(required_nodes_ptids[i]))+' '+seepage_phase+' '+seepage_value
    to_write_seepage.append(line+'\n')
    
  f_name = casename+'.'+side+'.'+time+'.flux_area'
  f = open(f_name, 'w')
  f.writelines(to_write_area)
  f.close()
 
  f_name = casename+'.'+side+'.'+time+'.seepage'
  f = open(f_name, 'w')
  f.writelines(to_write_seepage)
  f.close()

  return (1)

# parse command line
p = OptionParser(usage="""usage: %prog [options] <solid_vtu> <casename> <when> <side> <seepage_phase> <seepage_value>

Calculate flux area of nodes given in "side" and write to *.flux_area (COSFLOW bdy files) 

solid_vtu        - complete model (Ex:final_solid.vtu)
when             - before or after open cut
side             - zmax, zmin, xmin, xmax, ymin, ymax 



""")
p.add_option("-a", action="store_true", dest="ascii",  help="ASCII, instead of binary, .vtp output")
p.add_option("-v", action="store_true", dest="verbose",  help="Verbose")
(opts, args) = p.parse_args()

# get the com filename
if len(args) != 6:
   p.print_help()
   sys.exit(1)
(casename, solid_vtu, time, side, seepage_phase, seepage_value) = args

if opts.verbose: print "Reading solid"
r = vtk.vtkXMLUnstructuredGridReader()
r.SetFileName(solid_vtu)
r.Update()
solid = r.GetOutput()

if opts.verbose: print "Getting outer surface"
outer_surface = vtk.vtkDataSetSurfaceFilter()
outer_surface.PassThroughCellIdsOn() # useful later on
outer_surface.PassThroughPointIdsOn() # useful later on
outer_surface.SetInputConnection(solid.GetProducerPort())
outer_surface.Update()
outer_surface = outer_surface.GetOutput() # vtkPolyData

# Initial values
xmin = 10000000.0
ymin = 10000000.0
xmax = -10000000.0
ymax = -10000000.0
node_dict = {}
# Reading master node file
f = open(casename+'.nod', 'r')
for line in f:
  if not line.split():
    continue
  tuple = line.split()
  node_dict[tuple[0]] = tuple[0]
  if float(tuple[1]) < xmin:
    xmin = float(tuple[1])
  if float(tuple[2]) < ymin:
    ymin = float(tuple[2])
 
  if float(tuple[1]) > xmax:
    xmax = float(tuple[1])
  if float(tuple[2]) > ymax:
    ymax = float(tuple[2])

f.close()

if opts.verbose: print "Transforming"
rotated_surf = vtk.vtkTransformFilter()
trans = vtk.vtkTransform()
trans.PostMultiply()
  #trans.Translate(opts.xtrans, opts.ytrans, opts.ztrans)
if side == "ymin":
  trans.RotateX(-90)
  z_shift = abs(ymax-ymin) - 0.2
elif side == "ymax":
  trans.RotateX(90)
  z_shift = abs(ymax-ymin) - 0.2
elif side == "xmin":
  trans.RotateY(90)
  z_shift = abs(xmax-xmin) - 0.2
elif side == "xmax":
  trans.RotateY(-90)
  z_shift = abs(xmax-xmin) - 0.2
elif side == "zmin" or side == "zmax":
  z_shift = 0.2
else:
  print side+" is not a valid side"
  exit(0)
  

rotated_surf.SetInput(outer_surface)
rotated_surf.SetTransform(trans)
rotated_surf.Update()
rotated_surf = rotated_surf.GetOutput()
#if opts.verbose: print "Writing", vtp_file
#writer = vtk.vtkXMLPolyDataWriter()
#writer.SetFileName('xxx.vtp')
#if opts.ascii:
#   writer.SetDataModeToAscii()
#else:
#   writer.SetDataModeToBinary()
#writer.SetInputConnection(rotated_surf.GetProducerPort())
#writer.Write()
outdata = vtk.vtkTransformFilter()

outdata.SetInput(solid)
outdata.SetTransform(trans)
outdata.Update()

#if side == "zmin":
#  trans.RotateY(180)
#  outdata.SetTransform(trans)
#  outdata.Update()

shifter = vtk.vtkTransform()
shifter.Translate(0, 0, z_shift)
shifted_solid = vtk.vtkTransformFilter()
shifted_solid.SetInputConnection(outdata.GetOutput().GetProducerPort())
shifted_solid.SetTransform(shifter)
shifted_solid.Update()
shifted_solid = shifted_solid.GetOutput() # vtkUnstructuredData

#if opts.verbose: print "Writing", vtp_file
#writer = vtk.vtkXMLUnstructuredGridWriter()
#writer.SetFileName('yyy.vtu')
#if opts.ascii:
#   writer.SetDataModeToAscii()
#else:
#   writer.SetDataModeToBinary()
#writer.SetInputConnection(shifted_solid.GetProducerPort())
#writer.Write()

cf_nodes = rotated_surf.GetPointData().GetArray("cfnode")

# Cell locator for the solid
in_vol = vtk.vtkCellLocator()
in_vol.SetDataSet(shifted_solid)
in_vol.BuildLocator()

extracted_nodes_ptids = []

for ptid in range(rotated_surf.GetNumberOfPoints()):
  xyz = rotated_surf.GetPoint(ptid)
  if side == "zmin" or side == "zmax":
    if in_vol.FindCell(xyz) < 0:
      extracted_nodes_ptids.append(ptid)
  else:
    if in_vol.FindCell(xyz) >= 0:
      extracted_nodes_ptids.append(ptid)

required_nodes_ptids = []
for i in range(len(extracted_nodes_ptids)):
  val = str(cf_nodes.GetValue(extracted_nodes_ptids[i]))
  if val in node_dict:
    required_nodes_ptids.append(extracted_nodes_ptids[i])

# Create a required points list
xyz_data = []
for i in range(len(required_nodes_ptids)):
  xyz_data.append(rotated_surf.GetPoints().GetPoint(required_nodes_ptids[i]))

# Get poly data
if opts.verbose: print "Getting poly data"
clean_poly_data = generate_vtp()

# Get area
if opts.verbose: print "Calculating nodal areas"
area = flux_area(opts.verbose, clean_poly_data)

# Write flux area and seepage for COSFLOW
if opts.verbose: print "Writing COSFLOW flux_area and seepage files"
error = write_flux_area_and_seepage()

#if opts.verbose: print "Writing", vtp_file
#writer = vtk.vtkXMLPolyDataWriter()
#writer.SetFileName(casename+'.'+side+'.vtp')
#if opts.ascii:
#   writer.SetDataModeToAscii()
#else:
#   writer.SetDataModeToBinary()
#writer.SetInputConnection(clean_poly_data.GetProducerPort())
#writer.Write()

print "DONE"

sys.exit(0)

