#!/usr/bin/env python

import os
import sys
from optparse import OptionParser, OptionValueError
import subprocess
import tempfile

import numpy as np
import scipy.linalg

import vtk
from vtk.util import numpy_support


# so we can find our ../lib no matter how we are called
findbin = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(findbin + '/../lib')


def points_in_volume(verbose, mesh, vol, name, type):
   """ Determine whether points in "mesh" lie in "vol"

   @type  verbose:    bool
   @param verbose:    whether to print progress to screen
   @type  mesh:       vtkDataSet
   @param mesh:       Data containing points
   @type  vol:        vtkUnstructuredGrid
   @param vol:        the volume
   @type  type:       int
   @param type:       dictates type of output
   @type  name:       char
   @param name:       name given to the output
   @rtype:            tuple of vtkUnsignedCharArray (if type=0) or vtkFloatArray (if type=1), of length equal to number of points in mesh
   @return:           tuple value is 0 if outside volume, 1 if inside volume.  name is "name"
   """

   if opts.verbose: print "  Getting", name, "points"
   in_vol = vtk.vtkCellLocator()
   in_vol.SetDataSet(vol)
   in_vol.BuildLocator()

   if type == 0:
      pts = vtk.vtkUnsignedCharArray() # try to reduced file size, as pts is only 0 or 1.  So don't use IntArray.  Would be nice to use BitArray, but XMLUnstructuredGridWriter doesn't handle it.
   elif type == 1:
      pts = vtk.vtkFloatArray()
   pts.SetNumberOfValues(mesh.GetNumberOfPoints())
   pts.SetName(name)
   for ptid in range(mesh.GetNumberOfPoints()):
      xyz = mesh.GetPoint(ptid)
      if in_vol.FindCell(xyz) >= 0:
         pts.SetValue(ptid, 1)
      else:
         pts.SetValue(ptid, 0)

   return pts





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

   if verbose: print " Flux Area: Recording cell ids"
   cell_ids = vtk.vtkIdFilter()
   cell_ids.SetInputConnection(poly_data.GetProducerPort())
   cell_ids.CellIdsOn()
   cell_ids.Update()
   cell_ids = cell_ids.GetOutput() # vtkPolyData


   if verbose: print " Flux Area: Triangulating"
   tri_os = vtk.vtkTriangleFilter()
   tri_os.SetInputConnection(cell_ids.GetProducerPort())
   tri_os.Update()
   tri_os = tri_os.GetOutput() # vtkPolyData
   orig_cell_ids = tri_os.GetCellData().GetScalars("vtkIdFilter_Ids")


   if verbose: print " Flux Area: Calculating areas"
   cell_areas = [0.0 for cell in range(num_cells)]
   nodal_areas = [0.0 for nod in range(num_pts)]
   for cellid in range(tri_os.GetNumberOfCells()):
      the_cell = tri_os.GetCell(cellid)
      orig_cell_id = orig_cell_ids.GetValue(cellid)
      tri_area = the_cell.ComputeArea()
      for pt in range(3):
         ptid = the_cell.GetPointId(pt)
         # equally proportion triangle areas to the points
         nodal_areas[ptid] += tri_area/3.0
      cell_areas[orig_cell_id] += tri_area

   if verbose: print " Flux Area: Building cell and point area vtk array"
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



def metis_reorder(verbose, vtk_data):
   """ Reorder a mesh using METIS

   @type  verbose:    bool
   @param verbose:    whether to print progress to screen
   @type  vtk_data:   vtkDataSet
   @param vtk_data:   contains a mesh to be reordered using METIS.
   @rtype:            vtkDataSet
   @return:           the mesh which has been reordered
   """

   if verbose: print " Metis Reorder: Finding all adjacent nodes"
   atn = {} # stands for "attached to node".  note, in a finite element all nodes are attached to each other
   for ptid in range(vtk_data.GetNumberOfPoints()):
      atn[ptid+1] = [] # metis expects point ids starting at 1
   for cellid in range(vtk_data.GetNumberOfCells()):
      the_cell = vtk_data.GetCell(cellid)
      num_pts = the_cell.GetNumberOfPoints()
      noc = [the_cell.GetPointId(ptid)+1 for ptid in range(num_pts)] # stands for "nodes of cell" - the "+1" is so metis gets point ids starting at 1
      for ptid in range(num_pts):
         atn[noc[ptid]] += noc
      else:
         print "Reordering has not yet been implemented for this type of mesh"
         # it's easy enough to do: just need to build the atn for each node

   if verbose: print " Metis Reorder: Forming sets of links"
   num_links = 0 # this double-counts links, as a link between ptid A and B is also a link between ptid B and A
   for ptid in range(vtk_data.GetNumberOfPoints()):
      atn[ptid+1] = map(str, set(atn[ptid+1]) - set([ptid+1]))
      num_links += len(atn[ptid+1])

   if verbose: print " Metis Reorder: Outputting to METIS graph"
   f = open("for_metis.txt", 'w')
   f.write(str(vtk_data.GetNumberOfPoints()) + " " + str(num_links/2) + "\n")
   for ptid in range(vtk_data.GetNumberOfPoints()):
      f.write(" ".join(atn[ptid+1]) + "\n")
   f.close()
   if verbose: print " Metis Reorder: Doing the reordering"
   cmd = os.path.join("/apps/geomodel/svn/branch5/source/metis-4.0.3/", "onmetis")
   cmd += " " + "for_metis.txt"
   subprocess.call(cmd, shell=True)
   f = open("for_metis.txt.iperm", 'r')
   iperm = [int(nn) for nn in f.readlines()] # note, that .iperm file starts with ordering 0, not 1.
   f.close()

   if verbose: print "Building re-ordered mesh"
   sld = vtk.vtkUnstructuredGrid()
   rpts = vtk.vtkPoints()
   rpts.SetNumberOfPoints(vtk_data.GetNumberOfPoints())
   for ptid in range(vtk_data.GetNumberOfPoints()):
      xyz = vtk_data.GetPoint(ptid)
      rpts.InsertPoint(iperm[ptid], xyz)
   sld.SetPoints(rpts)
   sld.Update()
   ptids = vtk.vtkIdList()
   for cellid in range(vtk_data.GetNumberOfCells()):
      hex = vtk.vtkHexahedron().GetPointIds() # generalise this stuff for other cell types
      vtk_data.GetCellPoints(cellid, ptids)
      for pt in range(8):
         hex.SetId(pt, iperm[ptids.GetId(pt)])
      sld.InsertNextCell(12, hex)

   return sld



def metis_partition(verbose, parts, vtk_data):
   """ Partition a mesh using METIS

   @type  verbose:    bool
   @param verbose:    whether to print progress to screen
   @type  method:     parts
   @param method:     number of parts to split the mesh into
   @type  vtk_data:   vtkDataSet
   @param poly_data:  contains a mesh to be split using METIS.  This object needs GetNumberOfCells and GetCell methods (so might have to use GetOutput() before passing vtk_data to this function)
   @rtype:            vtkIntArray
   @return:           describes what processor (zero to num) the elements belong to.  name=cosflow_processor_number
   """

   num_cells = vtk_data.GetNumberOfCells()

   cpn = vtk.vtkIntArray()
   cpn.SetNumberOfValues(num_cells)
   cpn.SetName("cosflow_processor_number")


   if parts < 2 or num_cells < 1:
      if verbose: print " Metis Partition: Not partitioning: all cosflow_processor_number will be 0"
      for cellid in range(num_cells):
         cpn.SetValue(cellid, 0)

   else:
      if verbose: print " Metis Partition: Writing file containing elements"
      # check that all cells are of the same type
      cell_type = vtk_data.GetCellType(0)
      for cellid in range(1, num_cells):
         if cell_type != vtk_data.GetCellType(cellid):
            print " Metis version 4 cannot work on meshes with elements of more than one type"
            sys.exit(1)
      if cell_type == 5:
         metis_type = 1 # triangle
         num_nodes = 3
      elif cell_type == 10:
         metis_type = 2 # tet
         num_nodes = 4
      elif cell_type == 12:
         metis_type = 3 # hex
         num_nodes = 8
      elif cell_type == 9:
         metis_type = 4 # quad
         num_nodes = 4
      else:
         print " Metis version 4 can only work on tri, tet, hex or quad meshes"
         sys.exit(1)
         
      (f, fname) = tempfile.mkstemp(dir = ".")
      os.write(f, str(num_cells) + " " + str(metis_type) + "\n")
      for cellid in range(num_cells):
         the_cell = vtk_data.GetCell(cellid)
         os.write(f, ' '.join(map(str, [the_cell.GetPointId(ptid)+1 for ptid in range(num_nodes)]))+ "\n") # have to "+1" because metis needs nodes numbered from 1 up
      os.close(f)
         
      if verbose: print " Metis Partition: Doing the parition"
      cmd = os.path.join("/apps/geomodel/svn/branch5/source/metis-4.0.3/", "partnmesh")
      cmd += " " + fname + " " + str(parts)
      subprocess.call(cmd, shell=True)
         
      if verbose: print " Metis Partition: Making vtk array of processor number"
      f = open(fname + ".epart." + str(parts), 'r')
      cellid = 0
      for line in f:
         cpn.SetValue(cellid, int(line.split()[0]))
         cellid += 1
      f.close()

      os.unlink(fname)
      os.unlink(fname + ".epart." + str(parts))
      os.unlink(fname + ".npart." + str(parts))

   return cpn
      

   
   


   
# parse command line
p = OptionParser(usage="""usage: %prog [options] <basemesh> <casename>
Generates 3D mesh based on 2D <basemesh> (an exodus file).  Decorates it with pressures, materials, processor numbers, etc, and writes to casename.vtu (by default)

""")
p.add_option("-v", action="store_true",        dest="verbose",  help="Verbose")
p.add_option("-a", action="store_true", dest="ascii",  help="ASCII, instead of binary, .vtu output")
p.add_option("-o", action="store", dest="output", type="str", help="Output filename.  Default=casename.vtu.   Note that the output .vtu file will always have a fielddata variable called cosflow_casename with value \"casename\".")
p.add_option("-p", action="store", dest="partitions",  type=int, default=1, help="Partition the base mesh into this many partitions using Metis (default=%default)")
p.add_option("-m", action="store", dest="multikey", type="str", default="h", help="Multikey of the mesh partition (default=%default)")
p.add_option("--min_thick", action="store", dest="min_thickness",  type=float, default=-10.0, help="Minimum thickness of layers.  If a negative value is passed, then use the internal hard-coded minimum thickness tables.  (default=%default)")
p.add_option("--surf_dir", action="store", dest="surf_dir", type="str", default="../../preprocessing/surface_data", help="Directory to get the surface data from, relative to current dir (default=%default)")
p.add_option("--solid_dir", action="store", dest="solid_dir", type="str", default="../../preprocessing/solid_model", help="Directory to get the solids data from, relative to current dir (default=%default)")
p.add_option("-e", action="store", dest="elevation",  type=float, default=1000.0, help="Elevation at which pressure=0.  This is to form the initial pressure distribution (default=%default)")
p.add_option("--chopped", action="store_true", dest="chopped",  help="Chop off the mesh at MYC")
p.add_option("--subdivided", action="store_true", dest="subdivided",  help="Subdivide AQ3 and AQ4")
p.add_option("--extended_east", action="store_true", dest="extended_east",  help="Use the topography that has been extended to the east")
p.add_option("--split_weath", action="store_true", dest="split_weath",  help="Split the weath layer into weath and sp5")
p.add_option("--ys6", action="store_true", dest="ys6",  help="Use YS6 layer (with roof MYC_roof+100) instead of standard SP4.  In this case --subdivided is assumed, and any --chopped is ignored")
p.add_option("--ys6_1", action="store_true", dest="ys6_1",  help="Use YS6 layer (with roof MYC_roof+90) instead of standard SP4.  In this case --subdivided is assumed, and any --chopped is ignored")
p.add_option("--ys6_1_sp4", action="store_true", dest="ys6_1_sp4",  help="Use YS6 layer (with roof MYC_roof+90) and standard SP4.  In this case --subdivided is assumed, and any --chopped is ignored")
p.add_option("--ys6_true_sp4", action="store_true", dest="ys6_true_sp4",  help="Use YS6 layer with surfaces given by Chris on 7 Jan 2013, and standard SP4.  In this case --subdivided is assumed, and any --chopped is ignored")
p.add_option("-x", action="store", dest="xtrans", type="float", default=0.0, help="Translate the in .exo file by this amount in x direction")
p.add_option("-y", action="store", dest="ytrans", type="float", default=0.0, help="Translate the in .exo file by this amount in y direction")
p.add_option("-r", action="store_true", dest="reorder", help="Use onmetis to reorder 3D mesh")
(opts, args) = p.parse_args()

# get the com filename
if len(args) != 2:
   p.print_help()
   sys.exit(1)
(base_mesh_file, casename) = args

if opts.output:
   out_vtu = opts.output
else:
   out_vtu = casename + ".vtu"

if len(opts.multikey) != 1:
   print "The multikey must be a single letter"
   sys.exit(2)
   

surface_names = ["LTH_ext_minus_250", "LTH_ext_minus_5", "LTH_ext", "LTH_ext_plus_3.2", "LTH_ext_plus_6.5", "two_IRD_ext_plus_MDR_floor_ext", "half_IRD_ext_plus_MDR_floor_ext", "MDR_roof_ext", "KAT_floor_ext", "KAT_roof_ext", "MYC_floor_ext", "MYC_roof_ext", "MYC_roof_ext_plus_135", "MYC_roof_ext_plus_145", "all_100m_minus_10", "all_100m"] # real thing
default_solid_nums = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15] # if, when a hex is created, its nodes don't need to be shifted to avoid intersections, it is given this solid number
solid_names = ["aboveLTH_ext_minus_250", "aboveLTH_ext_minus_5", "aboveLTH_ext", "aboveLTH_ext_plus_3.2", "aboveLTH_ext_plus_6.5", "abovetwo_IRD_ext_plus_MDR_floor_ext", "abovehalf_IRD_ext_plus_MDR_floor_ext", "aboveMDR_roof_ext", "aboveKAT_floor_ext", "aboveKAT_roof_ext", "aboveMYC_floor_ext", "aboveMYC_roof_ext", "aboveMYC_roof_ext_plus_135", "aboveMYC_roof_ext_plus_145", "aboveall_100m_minus_10", "aboveall_100m"] # real thing
material_names = ["sp0", "floor", "devel", "roof", "aq1", "sp1", "aq2", "sp2", "kat", "aq3", "sp3", "aq4", "sp4", "aq5", "weath"] # real thing
if opts.chopped:
   surface_names = ["LTH_ext_minus_250", "MYC_roof_ext", "MYC_roof_ext_plus_135", "MYC_roof_ext_plus_145", "all_100m_minus_10", "all_100m"] # chopped
   default_solid_nums = [1, 2, 3, 4, 5]
   solid_names = ["aboveLTH_ext_minus_250", "aboveMYC_roof_ext", "aboveMYC_roof_ext_plus_135", "aboveMYC_roof_ext_plus_145", "aboveall_100m_minus_10", "aboveall_100m"] # chopped
   material_names = ["low_perm", "aq4", "sp4", "aq5", "weath"] # chopped
if opts.subdivided:
   surface_names = ["LTH_ext_minus_250", "LTH_ext_minus_5", "LTH_ext", "LTH_ext_plus_3.2", "LTH_ext_plus_6.5", "two_IRD_ext_plus_MDR_floor_ext", "half_IRD_ext_plus_MDR_floor_ext", "MDR_roof_ext", "KAT_floor_ext", "KAT_roof_ext", "two_KAT_roof_ext_plus_MYC_floor_ext", "KAT_roof_ext_plus_two_MYC_floor_ext", "MYC_floor_ext", "MYC_roof_ext", "MYC_roof_ext_plus_45", "MYC_roof_ext_plus_90", "MYC_roof_ext_plus_135", "MYC_roof_ext_plus_145", "all_100m_minus_10", "all_100m"] # real thing
   default_solid_nums = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 10, 10, 11, 12, 12, 12, 13, 14, 15] # if, when a hex is created, its nodes don't need to be shifted to avoid intersections, it is given this solid number
if opts.chopped and opts.subdivided:
   surface_names = ["LTH_ext_minus_250", "MYC_roof_ext", "MYC_roof_ext_plus_45", "MYC_roof_ext_plus_90", "MYC_roof_ext_plus_135", "MYC_roof_ext_plus_145", "all_100m_minus_10", "all_100m"]
   default_solid_nums = [1, 2, 2, 2, 3, 4, 5]
   solid_names = ["aboveLTH_ext_minus_250", "aboveMYC_roof_ext", "aboveMYC_roof_ext_plus_135", "aboveMYC_roof_ext_plus_145", "aboveall_100m_minus_10", "aboveall_100m"]
   material_names = ["low_perm", "aq4", "sp4", "aq5", "weath"]
if opts.ys6:
   surface_names = ["LTH_ext_minus_250", "LTH_ext_minus_5", "LTH_ext", "LTH_ext_plus_3.2", "LTH_ext_plus_6.5", "two_IRD_ext_plus_MDR_floor_ext", "half_IRD_ext_plus_MDR_floor_ext", "MDR_roof_ext", "KAT_floor_ext", "KAT_roof_ext", "two_KAT_roof_ext_plus_MYC_floor_ext", "KAT_roof_ext_plus_two_MYC_floor_ext", "MYC_floor_ext", "MYC_roof_ext", "MYC_roof_ext_plus_45", "MYC_roof_ext_plus_90", "MYC_roof_ext_plus_100", "MYC_roof_ext_plus_135", "all_100m_minus_10", "all_100m"] # real thing
   default_solid_nums = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 10, 10, 11, 12, 12, 13, 14, 14, 15] # if, when a hex is created, its nodes don't need to be shifted to avoid intersections, it is given this solid number
   solid_names = ["aboveLTH_ext_minus_250", "aboveLTH_ext_minus_5", "aboveLTH_ext", "aboveLTH_ext_plus_3.2", "aboveLTH_ext_plus_6.5", "abovetwo_IRD_ext_plus_MDR_floor_ext", "abovehalf_IRD_ext_plus_MDR_floor_ext", "aboveMDR_roof_ext", "aboveKAT_floor_ext", "aboveKAT_roof_ext", "aboveMYC_floor_ext", "aboveMYC_roof_ext", "aboveMYC_roof_ext_plus_90", "aboveMYC_roof_ext_plus_100", "aboveall_100m_minus_10", "aboveall_100m"]
if opts.ys6_1:
   surface_names = ["LTH_ext_minus_250", "LTH_ext_minus_5", "LTH_ext", "LTH_ext_plus_3.2", "LTH_ext_plus_6.5", "two_IRD_ext_plus_MDR_floor_ext", "half_IRD_ext_plus_MDR_floor_ext", "MDR_roof_ext", "KAT_floor_ext", "KAT_roof_ext", "two_KAT_roof_ext_plus_MYC_floor_ext", "KAT_roof_ext_plus_two_MYC_floor_ext", "MYC_floor_ext", "MYC_roof_ext", "MYC_roof_ext_plus_45", "MYC_roof_ext_plus_80", "MYC_roof_ext_plus_90", "MYC_roof_ext_plus_135", "all_100m_minus_10", "all_100m"] # real thing
   default_solid_nums = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 10, 10, 11, 12, 12, 13, 14, 14, 15] # if, when a hex is created, its nodes don't need to be shifted to avoid intersections, it is given this solid number
   solid_names = ["aboveLTH_ext_minus_250", "aboveLTH_ext_minus_5", "aboveLTH_ext", "aboveLTH_ext_plus_3.2", "aboveLTH_ext_plus_6.5", "abovetwo_IRD_ext_plus_MDR_floor_ext", "abovehalf_IRD_ext_plus_MDR_floor_ext", "aboveMDR_roof_ext", "aboveKAT_floor_ext", "aboveKAT_roof_ext", "aboveMYC_floor_ext", "aboveMYC_roof_ext", "aboveMYC_roof_ext_plus_80", "aboveMYC_roof_ext_plus_90", "aboveall_100m_minus_10", "aboveall_100m"]
if opts.ys6_1_sp4:
   surface_names = ["LTH_ext_minus_250", "LTH_ext_minus_5", "LTH_ext", "LTH_ext_plus_3.2", "LTH_ext_plus_6.5", "two_IRD_ext_plus_MDR_floor_ext", "half_IRD_ext_plus_MDR_floor_ext", "MDR_roof_ext", "KAT_floor_ext", "KAT_roof_ext", "two_KAT_roof_ext_plus_MYC_floor_ext", "KAT_roof_ext_plus_two_MYC_floor_ext", "MYC_floor_ext", "MYC_roof_ext", "MYC_roof_ext_plus_45", "MYC_roof_ext_plus_80", "MYC_roof_ext_plus_90", "MYC_roof_ext_plus_135", "MYC_roof_ext_plus_145", "all_100m_minus_10", "all_100m"] # real thing
   default_solid_nums = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 10, 10, 11, 12, 12, 13, 14, 15, 16, 17] # if, when a hex is created, its nodes don't need to be shifted to avoid intersections, it is given this solid number
   solid_names = ["aboveLTH_ext_minus_250", "aboveLTH_ext_minus_5", "aboveLTH_ext", "aboveLTH_ext_plus_3.2", "aboveLTH_ext_plus_6.5", "abovetwo_IRD_ext_plus_MDR_floor_ext", "abovehalf_IRD_ext_plus_MDR_floor_ext", "aboveMDR_roof_ext", "aboveKAT_floor_ext", "aboveKAT_roof_ext", "aboveMYC_floor_ext", "aboveMYC_roof_ext", "aboveMYC_roof_ext_plus_80", "aboveMYC_roof_ext_plus_90", "aboveMYC_roof_ext_plus_135", "aboveMYC_roof_ext_plus_145", "aboveall_100m_minus_10", "aboveall_100m"]
   material_names = ["sp0", "floor", "devel", "roof", "aq1", "sp1", "aq2", "sp2", "kat", "aq3", "sp3", "aq4", "ys6", "aq4_up", "sp4", "aq5", "weath"] # real thing
if opts.ys6_true_sp4:
   surface_names = ["LTH_ext_minus_250", "LTH_ext_minus_5", "LTH_ext", "LTH_ext_plus_3.2", "LTH_ext_plus_6.5", "two_IRD_ext_plus_MDR_floor_ext", "half_IRD_ext_plus_MDR_floor_ext", "MDR_roof_ext", "KAT_floor_ext", "KAT_roof_ext", "two_KAT_roof_ext_plus_MYC_floor_ext", "KAT_roof_ext_plus_two_MYC_floor_ext", "MYC_floor_ext", "MYC_roof_ext", "MYC_roof_ext_plus_45", "ys6_floor_ext", "ys6_roof_ext", "MYC_roof_ext_plus_135", "MYC_roof_ext_plus_145", "all_100m_minus_10", "all_100m"] # real thing
   default_solid_nums = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 10, 10, 11, 12, 12, 13, 14, 15, 16, 17] # if, when a hex is created, its nodes don't need to be shifted to avoid intersections, it is given this solid number
   solid_names = ["aboveLTH_ext_minus_250", "aboveLTH_ext_minus_5", "aboveLTH_ext", "aboveLTH_ext_plus_3.2", "aboveLTH_ext_plus_6.5", "abovetwo_IRD_ext_plus_MDR_floor_ext", "abovehalf_IRD_ext_plus_MDR_floor_ext", "aboveMDR_roof_ext", "aboveKAT_floor_ext", "aboveKAT_roof_ext", "aboveMYC_floor_ext", "aboveMYC_roof_ext", "aboveys6_floor_ext", "aboveys6_roof_ext", "aboveMYC_roof_ext_plus_135", "aboveMYC_roof_ext_plus_145", "aboveall_100m_minus_10", "aboveall_100m"]
   material_names = ["sp0", "floor", "devel", "roof", "aq1", "sp1", "aq2", "sp2", "kat", "aq3", "sp3", "aq4", "ys6", "aq4_up", "sp4", "aq5", "weath"] # real thing
if opts.extended_east:
   surface_names[-1] = "topo_extended_east"
   surface_names[-2] = "topo_extended_east_minus_10"
   solid_names[-1] = "abovetopo_extended_east"
   solid_names[-2] = "abovetopo_extended_east_minus_10"
if opts.split_weath:
   final_surf_name = surface_names[-1]
   surface_names[-1] = final_surf_name + "_minus_5"
   surface_names.append(final_surf_name)
   final_material_name = material_names[-1]
   material_names[-1] = "sp5"
   material_names.append(final_material_name)
   final_solid_name = solid_names[-1]
   solid_names[-1] = final_solid_name + "_minus_5"
   solid_names.append(final_solid_name)
   default_solid_nums.append(default_solid_nums[-1]+1)


num_surfaces = len(surface_names)
num_solids = len(solid_names)
if opts.min_thickness < 0:
   min_thick_table = [50, 4.5, 3.15, 3.25, 10, 6.0, 10, 4.0, 4.0, 10, 4.0, 10, 9.5, 10, 9.5, 11.111] # real thing, last entry is a dummy - must be positive
   if opts.chopped:
      min_thick_table = [50, 10, 9.5, 10, 9.5, 11.111] # chopped last entry is a dummy - must be positive
   if opts.subdivided:
      min_thick_table = [50, 4.5, 3.15, 3.25, 10, 6.0, 10, 3.0, 2.5, 3.3, 3.3, 3.3, 3.0, 3.3, 3.3, 3.3, 9.5, 10, 9.5, 11.111] # real thing, last entry is a dummy - must be positive
   if opts.chopped and opts.subdivided:
      min_thick_table = [50, 3.3, 3.3, 3.3, 9.5, 10, 9.5, 11.111] # last entry is a dummy - must be positive
   if opts.ys6 or opts.ys6_1:
      min_thick_table = [50, 4.5, 3.15, 3.25, 10, 6.0, 10, 3.0, 2.5, 3.3, 3.3, 3.3, 3.0, 3.3, 3.3, 9.5, 3, 3, 9.5, 11.111] # real thing, last entry is a dummy - must be positive
   if opts.ys6_1_sp4:
      min_thick_table = [50, 4.5, 3.15, 3.25, 10, 6.0, 10, 3.0, 2.5, 3.3, 3.3, 3.3, 3.011, 3.3012, 3.3012, 6.013, 3.014, 3.015, 3.016, 9.5017, 11.111] # last entry is a dummy - must be positive
   if opts.ys6_true_sp4:
      min_thick_table = [50, 4.5, 3.15, 3.25, 10, 6.0, 10, 3.0, 2.5, 3.3, 3.3, 3.3, 3.011, 3.3012, 3.3012, 2.0, 3.014, 3.015, 3.016, 9.5017, 11.111] # last entry is a dummy - must be positive
   if opts.split_weath:
      final_thickness = min_thick_table[-2]
      min_thick_table[-2] = final_thickness/2
      min_thick_table[-1] = final_thickness/2
      min_thick_table.append(11.111)
else:
   min_thick_table = [opts.min_thickness]*num_surfaces
surface_flat = [vtk.vtkTransformFilter() for x in range(num_surfaces)]
base_mesh_on_surface = [vtk.vtkWarpScalar() for x in range(num_surfaces)]
loc = [vtk.vtkCellLocator() for x in range(num_surfaces)]

if len(material_names) != num_solids - 1:
   print "material_names and solid_names have incompatable lengths"
   sys.exit(3)
if num_surfaces != len(default_solid_nums) + 1:
   print "default_solid_nums and number of surfaces are incompatable"
   sys.exit(4)

# this is a transform to get rid of the z values
# the data must be flat for the interpolation to work
# otherwise vtk thinks the xyz points don't lie on
# the vtp surface and no interpolation occurs
flattener = vtk.vtkTransform()
flattener.Scale(1.0, 1.0, 0.0)


if opts.verbose: print "Reading", base_mesh_file
base_mesh_f = vtk.vtkExodusReader()
base_mesh_f.SetFileName(base_mesh_file)
base_mesh_f.Update()
base_mesh_f = base_mesh_f.GetOutput()   # now base is vtkUnstructuredGrid

if opts.xtrans != 0.0 or opts.ytrans != 0.0:
   if opts.verbose: print "Translating by", str(opts.xtrans), str(opts.ytrans)
   base = vtk.vtkTransformFilter()
   transe = vtk.vtkTransform()
   transe.Translate(opts.xtrans, opts.ytrans, 0.0)
   base.SetInput(base_mesh_f)
   base.SetTransform(transe)
   base.Update()
   base = base.GetOutput()
else:
   base = base_mesh_f




if opts.verbose: print "Checking orientation"
clockwise = {}
clockw = [True]*3
for cellid in range(base.GetNumberOfCells()):
   this_ele = base.GetCell(cellid)
   xyz = []
   for pt in range(4):
      xyz.append(base.GetPoint(this_ele.GetPointId(pt)))
   for pt in range(1,4):
      a = [xyz[pt-1][0] - xyz[pt][0], xyz[pt-1][1] - xyz[pt][1]]
      b = [xyz[(pt+1)%4][0] - xyz[pt][0], xyz[(pt+1)%4][1] - xyz[pt][1]]
      clockw[pt-1] = np.cross(a,b) <= 0
   if all(clockw):
      clockwise[cellid] = True
   elif not any(clockw):
      clockwise[cellid] = False
   else:
      print "Cannot determine orientation of cell", cellid
      print "This has points", xyz
      if cellid > 1:
         print "Cosflow may give crazy results for this mesh"
         print "Using orientation of", cellid-1
         clockwise[cellid] = clockwise[cellid-1]
      else:
         sys.exit(5)


if opts.verbose: print "Flattening", base_mesh_file, "to remove z values of points"
base_flat = vtk.vtkTransformFilter()
base_flat.SetInputConnection(base.GetProducerPort())
base_flat.SetTransform(flattener)
base_flat.Update()
base_flat = base_flat.GetOutput() # now base_flat is vtkUnstructuredGrid




# read vtp files 
for surf_num in range(num_surfaces):
   surf_name = surface_names[surf_num] + ".vtp"
   if opts.verbose: print "Reading", surf_name
   r = vtk.vtkXMLPolyDataReader()
   r.SetFileName(os.path.join(opts.surf_dir, surf_name))
   r.Update()

   if opts.verbose: print "Flattening", surf_name, "to remove z values of points"
   surface_flat[surf_num].SetInputConnection(r.GetOutput().GetProducerPort())
   surface_flat[surf_num].SetTransform(flattener)
   surface_flat[surf_num].Update()
   surface_flat[surf_num] = surface_flat[surf_num].GetOutput() # vtkPolyData


   if opts.verbose: print "Interpolating", base_mesh_file, "onto", surf_name
   probe = vtk.vtkProbeFilter()
   probe.SetInput(base_flat)
   probe.SetSource(surface_flat[surf_num])
   probe.Update()

   if opts.verbose: print "Putting elevation into", base_mesh_file
   base_mesh_on_surface[surf_num].SetInputConnection(probe.GetOutput().GetProducerPort())
   base_mesh_on_surface[surf_num].SetInputArrayToProcess(0, 0, 0, 0, "Elevation")
   base_mesh_on_surface[surf_num].Update()
   base_mesh_on_surface[surf_num] = base_mesh_on_surface[surf_num].GetOutput() # now this is vtkUnstructuredGrid




num_pts = base.GetNumberOfPoints()
if opts.verbose: print "Found", num_pts, "points: extracting them"
pts = vtk.vtkPoints()
pts.SetNumberOfPoints(num_surfaces*num_pts)


if opts.verbose:
   if opts.min_thickness > 0:
      print "Putting unentangled points on surfaces into vtk points.  Assuming that", surface_names[num_surfaces-1], "is correct and using min layer thickness", opts.min_thickness
   else:
      print "Putting unentangled points on surfaces into vtk points.  Assuming that", surface_names[num_surfaces-1], "is correct and using min layer thickness table"
z_max = base_mesh_on_surface[num_surfaces-1].GetBounds()[5] + 2*min_thick_table[-1]
pt_shifted = [False]*num_surfaces*num_pts # record if a point has been shifted.
for ptid in range(num_pts):
   z_current = z_max
   for surf_num in range(num_surfaces-1, -1, -1):
      xyz = base_mesh_on_surface[surf_num].GetPoint(ptid)
      if xyz[2] >= z_current - min_thick_table[surf_num]:
         pt_shifted[ptid + surf_num*num_pts] = True
         pts.InsertPoint(ptid + surf_num*num_pts, xyz[0], xyz[1], z_current - min_thick_table[surf_num])
         z_current -= min_thick_table[surf_num]
      else:
         pts.InsertPoint(ptid + surf_num*num_pts, xyz)
         z_current = xyz[2]



if opts.verbose: print "Creating solid unstructured grid and adding points"
usolid = vtk.vtkUnstructuredGrid()
usolid.SetPoints(pts)
usolid.Update()



material_number = vtk.vtkIntArray()
material_number.SetNumberOfValues(base.GetNumberOfCells()*(num_surfaces-1))
material_number.SetName("material")
if opts.verbose: print "Inserting", base.GetNumberOfCells()*(num_surfaces-1), "hexahedra in", out_vtu
ptids = vtk.vtkIdList()
current_hexid = 0
for cellid in range(base.GetNumberOfCells()):
   if base.GetCellType(cellid) != 9:
      print "Cannot form a valid cosflow mesh with non-quads in the base"
      print cellid, base.GetCellType(cellid)
      sys.exit(2)
   base.GetCellPoints(cellid, ptids)
   quad_pts = [ptids.GetId(pt) for pt in range(4)]
   for surf_num in range(num_surfaces-1):
      hex = vtk.vtkHexahedron().GetPointIds()
      contains_shifted_points = False
      if clockwise[cellid]:
         for pt in range(4):
            hex.SetId(pt, quad_pts[pt]+surf_num*num_pts)
            hex.SetId(pt+4, quad_pts[pt]+(surf_num+1)*num_pts)
            contains_shifted_points = contains_shifted_points or pt_shifted[quad_pts[pt]+surf_num*num_pts] or pt_shifted[quad_pts[pt]+(surf_num+1)*num_pts]
      else:
         hex.SetId(0, quad_pts[3]+surf_num*num_pts)
         hex.SetId(1, quad_pts[2]+surf_num*num_pts)
         hex.SetId(2, quad_pts[1]+surf_num*num_pts)
         hex.SetId(3, quad_pts[0]+surf_num*num_pts)
         hex.SetId(4, quad_pts[3]+(surf_num+1)*num_pts)
         hex.SetId(5, quad_pts[2]+(surf_num+1)*num_pts)
         hex.SetId(6, quad_pts[1]+(surf_num+1)*num_pts)
         hex.SetId(7, quad_pts[0]+(surf_num+1)*num_pts)
      usolid.InsertNextCell(12, hex)
      contains_shifted_points = any([pt_shifted[qpt+surf_num*num_pts] or pt_shifted[qpt+(surf_num+1)*num_pts] for qpt in quad_pts])
      if contains_shifted_points:
         material_number.SetValue(current_hexid, 0)
      else:
         material_number.SetValue(current_hexid, default_solid_nums[surf_num])
      current_hexid += 1


if not opts.reorder:
   solid = usolid
else:
   solid = metis_reorder(opts.verbose, usolid)
   


   
solid.GetCellData().AddArray(metis_partition(opts.verbose, opts.partitions, solid))




if opts.verbose: print "Calculating the centroid of each element"
cc = vtk.vtkCellCenters()
cc.SetInputConnection(solid.GetProducerPort())
cc.Update()
cc = cc.GetOutput()  # vtkUnstructuredGrid


for sol in range(num_solids-1):
   sol_name = solid_names[sol] + ".vtu"
   if opts.verbose: print "Reading", sol_name, "and building the locator"
   r = vtk.vtkXMLUnstructuredGridReader()
   r.SetFileName(os.path.join(opts.solid_dir, sol_name))
   r.Update()
   loc[sol].SetDataSet(r.GetOutput())
   loc[sol].BuildLocator()



if opts.verbose: print "Prescribing material numbers to the mesh"
five_percent = cc.GetNumberOfPoints()/20
for cellid in range(cc.GetNumberOfPoints()):
   if (cellid+1)%five_percent == 0:
      print " ", str(int((100.0*cellid)/cc.GetNumberOfPoints())) + "%"
   if material_number.GetValue(cellid) == 0:
      xyz = cc.GetPoint(cellid)
      for sol in range(num_solids-2, -1, -1):
         if loc[sol].FindCell(xyz) >= 0:
            material_number.SetValue(cellid, sol + 1)
            break
      else:
         # must have been really bad insersections, or min_thickness set too large - ele is below the bottom surface.  put to material=0
         material_number.SetValue(cellid, 1)
solid.GetCellData().AddArray(material_number)




(x_min, x_max, y_min, y_max, z_min, z_max) = solid.GetBounds()[0:6]
buffer = 1  # to allow for roundoff
(p_min, p_max, q_min, q_max) = (x_min + buffer, x_max - buffer, y_min + buffer, y_max - buffer)
if opts.verbose: print "Getting outer perimeter"
perim = vtk.vtkUnsignedCharArray()
perim.SetNumberOfValues(solid.GetNumberOfPoints())
perim.SetName("cosflow_perim")
for ptid in range(solid.GetNumberOfPoints()):
   xyz = solid.GetPoint(ptid)
   if xyz[0] <= p_min or xyz[0] >= p_max or xyz[1] <= q_min or xyz[1] >= q_max:
      perim.SetValue(ptid, 1)
   else:
      perim.SetValue(ptid, 0)
solid.GetPointData().AddArray(perim)


##############################################################
#
# Now work with the outer surface stuff.
# Don't add things to solid till the very end
#
##############################################################


if opts.verbose: print "Getting outer surface"
outer_surface = vtk.vtkDataSetSurfaceFilter()
outer_surface.PassThroughCellIdsOn() # useful later on
outer_surface.PassThroughPointIdsOn() # useful later on
outer_surface.SetInputConnection(solid.GetProducerPort())
outer_surface.Update()
outer_surface = outer_surface.GetOutput() # vtkPolyData



if opts.verbose: print "Reading top volume"
r = vtk.vtkXMLUnstructuredGridReader()
r.SetFileName(os.path.join(opts.solid_dir, "above" + surface_names[-1] + ".vtu"))
r.Update()
r = r.GetOutput()
shifter = vtk.vtkTransform()
shifter.Translate(0.0, 0.0, -min(min_thick_table)/2)
vol_at_top = vtk.vtkTransformFilter()
vol_at_top.SetInputConnection(r.GetProducerPort())
vol_at_top.SetTransform(shifter)
vol_at_top.Update()
vol_at_top = vol_at_top.GetOutput() # vtkUnstructuredGrid



if opts.verbose: print "Reading bottom volume"
r = vtk.vtkXMLUnstructuredGridReader()
r.SetFileName(os.path.join(opts.solid_dir, "below" + surface_names[0] + ".vtu"))
r.Update()
r = r.GetOutput()
shifter = vtk.vtkTransform()
shifter.Translate(0.0, 0.0, min(min_thick_table)/2)
vol_at_bot = vtk.vtkTransformFilter()
vol_at_bot.SetInputConnection(r.GetProducerPort())
vol_at_bot.SetTransform(shifter)
vol_at_bot.Update()
vol_at_bot = vol_at_bot.GetOutput() # vtkUnstructuredGrid


if opts.verbose: print "Getting top surface"
outer_surface.GetPointData().AddArray(points_in_volume(opts.verbose, outer_surface, vol_at_top, "in_top", 1))
if opts.verbose: print "Getting bottom surface"
outer_surface.GetPointData().AddArray(points_in_volume(opts.verbose, outer_surface, vol_at_bot, "in_bot", 1))
p2c = vtk.vtkPointDataToCellData()
p2c.SetInputConnection(outer_surface.GetProducerPort())
p2c.PassPointDataOn()
p2c.Update()
on_top = vtk.vtkThreshold()
on_top.SetInputConnection(p2c.GetOutput().GetProducerPort())
on_top.SetInputArrayToProcess(0, 0, 0, 1, "in_top")
on_top.ThresholdBetween(0.9, 1.1) 
on_top.Update()
on_top_poly = vtk.vtkGeometryFilter()
on_top_poly.SetInputConnection(on_top.GetOutput().GetProducerPort())
on_top_poly.Update()
on_top_poly = on_top_poly.GetOutput() # vtkPolyData
on_bot = vtk.vtkThreshold()
on_bot.SetInputConnection(p2c.GetOutput().GetProducerPort())
on_bot.SetInputArrayToProcess(0, 0, 0, 1, "in_bot")
on_bot.ThresholdBetween(0.9, 1.1) 
on_bot.Update()
on_bot_poly = vtk.vtkGeometryFilter()
on_bot_poly.SetInputConnection(on_bot.GetOutput().GetProducerPort())
on_bot_poly.Update()
on_bot_poly = on_bot_poly.GetOutput() # vtkPolyData





if opts.verbose: print "Getting flux areas of top nodes"
areas = flux_area(opts.verbose, on_top_poly)
on_top_poly.GetPointData().AddArray(areas[0])
on_top_poly.GetCellData().AddArray(areas[1])




if opts.verbose: print "Writing on_top_poly.vtp"
writer = vtk.vtkXMLPolyDataWriter()
writer.SetFileName("on_top_poly.vtp")
if opts.ascii:
   writer.SetDataModeToAscii()
else:
   writer.SetDataModeToBinary()
writer.SetInputConnection(on_top_poly.GetProducerPort())
writer.Write()


##############################################################
#
# Now prepare to record various things into solid
#
##############################################################



if opts.verbose: print "Recording top nodes to put in solid"
uca_zmax = vtk.vtkUnsignedCharArray()
uca_zmax.SetNumberOfValues(solid.GetNumberOfPoints())
uca_zmax.SetName("cosflow_zmax_points")
for ptid in range(solid.GetNumberOfPoints()):
   # initialise them all to zero
   uca_zmax.SetValue(ptid, 0)
orig_ptids = on_top_poly.GetPointData().GetArray("vtkOriginalPointIds")
for ptid in range(on_top_poly.GetNumberOfPoints()):
   uca_zmax.SetValue(orig_ptids.GetValue(ptid), 1)



if opts.verbose: print "Recording bot nodes to put in solid"
uca_zmin = vtk.vtkUnsignedCharArray()
uca_zmin.SetNumberOfValues(solid.GetNumberOfPoints())
uca_zmin.SetName("cosflow_zmin_points")
for ptid in range(solid.GetNumberOfPoints()):
   # initialise them all to zero
   uca_zmin.SetValue(ptid, 0)
orig_ptids = on_bot_poly.GetPointData().GetArray("vtkOriginalPointIds")
for ptid in range(on_bot_poly.GetNumberOfPoints()):
   uca_zmin.SetValue(orig_ptids.GetValue(ptid), 1)



if opts.verbose: print "Recording nodal areas to put into solid"
fa_nodal_area = vtk.vtkFloatArray()
fa_nodal_area.SetNumberOfValues(solid.GetNumberOfPoints())
fa_nodal_area.SetName("cosflow_nodal_flux_area")
for ptid in range(solid.GetNumberOfPoints()):
   # initialise them all to zero
   fa_nodal_area.SetValue(ptid, 0.0)
orig_ptids = on_top_poly.GetPointData().GetArray("vtkOriginalPointIds")
nodal_flux_areas = on_top_poly.GetPointData().GetArray("cosflow_nodal_flux_area")
for ptid in range(on_top_poly.GetNumberOfPoints()):
   fa_nodal_area.SetValue(orig_ptids.GetValue(ptid), nodal_flux_areas.GetValue(ptid))


if opts.verbose: print "Recording elemental areas to put into solid"
fa_ele_area = vtk.vtkFloatArray()
fa_ele_area.SetNumberOfValues(solid.GetNumberOfCells())
fa_ele_area.SetName("cosflow_ele_flux_area")
for cellid in range(solid.GetNumberOfCells()):
   # initialise them all to zero
   fa_ele_area.SetValue(cellid, 0.0)
orig_cellids = on_top_poly.GetCellData().GetArray("vtkOriginalCellIds")
ele_flux_areas = on_top_poly.GetCellData().GetArray("cosflow_ele_flux_area")
for cellid in range(on_top_poly.GetNumberOfCells()):
   fa_ele_area.SetValue(orig_cellids.GetValue(cellid), ele_flux_areas.GetValue(cellid))



if opts.verbose: print "Building initial pressure distribution"
ini_pp = vtk.vtkFloatArray()
ini_pp.SetNumberOfValues(solid.GetNumberOfPoints())
ini_pp.SetName("cosflow_initial_porepressure")
for ptid in range(solid.GetNumberOfPoints()):
   xyz = solid.GetPoint(ptid)
   ini_pp.SetValue(ptid, max( (opts.elevation - xyz[2])*10*1000, -1.0E5)) # pressure=-1.0E+5 corresponds to saturation=0.1, which is immobile saturation


solid.GetPointData().AddArray(uca_zmax)
solid.GetPointData().AddArray(uca_zmin)
solid.GetPointData().AddArray(fa_nodal_area)
solid.GetPointData().AddArray(ini_pp)
solid.GetCellData().AddArray(fa_ele_area)

command_used = vtk.vtkStringArray()
command_used.SetName("provenance")
command_used.InsertNextValue(" ".join(sys.argv))
solid.GetFieldData().AddArray(command_used)

cn = vtk.vtkStringArray()
cn.SetName("cosflow_casename")
cn.InsertNextValue(casename)
solid.GetFieldData().AddArray(cn)

mk = vtk.vtkStringArray()
mk.SetName("cosflow_multikey")
mk.InsertNextValue(opts.multikey)
solid.GetFieldData().AddArray(mk)

mats = vtk.vtkStringArray()
mats.SetName("material_names")
mats.InsertNextValue("None")
for mat in material_names:
   mats.InsertNextValue(mat)
solid.GetFieldData().AddArray(mats)


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

