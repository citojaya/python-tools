#!/tools/python/2.6-x86_64/bin/python

import os
import sys

import re
from optparse import OptionParser, OptionValueError

import numpy as np

import glob
import vtk
from vtk.util import numpy_support

# so we can find our ../lib no matter how we are called
findbin = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(findbin + '/../lib')


# parse command line
p = OptionParser(usage="""usage: %prog [options] <casename> <out_vtu>
Generates 3D mesh based on <casename>.nod and <casename>.ele files

Materials are included.
""")
p.add_option("-v", action="store_true", dest="verbose",  help="Verbose")
p.add_option("-a", action="store_true", dest="ascii",  help="ASCII, instead of binary, .vtu output")
p.add_option("-m", action="store", dest="multikey", type="str", help="Multikey of the mesh partition.  If given, then <out_vtu> will contain cosflow_multikey")
(opts, args) = p.parse_args()


# Get the arguments
if len(args) != 2:
  p.print_help()
  print "\nProgram must be called with 2 arguments"
  sys.exit(1)
(casename, out_file) = args

if not out_file.endswith(".vtu"):
  print "Out_vtu must end with .vtu"
  sys.exit(2)



# Read the .nod file, create points, and cfnode PointData
cfnode = vtk.vtkIntArray()
cfnode.SetName("cfnode")
points = vtk.vtkPoints()
ptid_num = {} # given a cosflow node number (as a string), this returns the vtk point ID

if opts.verbose: print "Reading", casename + ".nod and building vtk points"
f = open(casename + ".nod", 'r')
ptid = 0
for line in f:
  nxyz = line.split()
  if not nxyz:
    continue
  points.InsertNextPoint(float(nxyz[1]), float(nxyz[2]), float(nxyz[3]))
  cfnode.InsertNextValue(int(nxyz[0]))
  ptid_num[nxyz[0]] = ptid
  ptid += 1
f.close()


# add points and cfnode to output
output = vtk.vtkUnstructuredGrid()
output.SetPoints(points)
output.GetPointData().AddArray(cfnode)

 



# Read .ele file, create elements, create materials
largest_mat_number = 0
mat_number_dict = {} # mat_number[material_string] = number
mat_name = vtk.vtkStringArray()
mat_name.SetName("cosflow_material_name")
mat_number = vtk.vtkIntArray()
mat_number.SetName("cosflow_material")
cfele = vtk.vtkIntArray()
cfele.SetName("cfele")
cellid_num = {} # given a cosflow ele number (as a string), this returns the vtk cell ID

if opts.verbose: print "Reading", casename+".ele and building the mesh"
f = open(casename+'.ele', 'r')
cellid = 0
for line in f:
  tuple = line.split()
  if not tuple:
    continue

  # get the material number corresponding to the material name
  if tuple[0] not in mat_number_dict:
    mat_number_dict[tuple[0]] = largest_mat_number
    mat_name.InsertNextValue(tuple[0])
    largest_mat_number += 1
  mat_number.InsertNextValue(mat_number_dict[tuple[0]])

  # set the cfele number
  cfele.InsertNextValue(int(tuple[1]))

  # create the element
  if len(tuple) == 6:
    quad = vtk.vtkQuad().GetPointIds()
    quad.SetId(0, ptid_num[tuple[2]])
    quad.SetId(1, ptid_num[tuple[3]])
    quad.SetId(2, ptid_num[tuple[4]])
    quad.SetId(3, ptid_num[tuple[5]])
    output.InsertNextCell(9, quad)
  elif len(tuple) == 10:
    hex = vtk.vtkHexahedron().GetPointIds()
    hex.SetId(0, ptid_num[tuple[2]])
    hex.SetId(1, ptid_num[tuple[3]])
    hex.SetId(2, ptid_num[tuple[4]])
    hex.SetId(3, ptid_num[tuple[5]])
    hex.SetId(4, ptid_num[tuple[6]])
    hex.SetId(5, ptid_num[tuple[7]])
    hex.SetId(6, ptid_num[tuple[8]])
    hex.SetId(7, ptid_num[tuple[9]])
    output.InsertNextCell(12, hex)

  cellid_num[tuple[1]] = cellid
  cellid += 1
     
f.close()


# add cfele and mat_number arrays to output mesh
output.GetCellData().AddArray(cfele)
output.GetCellData().AddArray(mat_number)
output.GetCellData().AddArray(mat_name)




# If required, read the multiproc .ele files and assign CPU numbers
if opts.multikey:
  if opts.verbose: print "Reading", casename + opts.multikey + "*.ele files and assigining CPU numbers to the elements"

  cpu_array = vtk.vtkIntArray()
  cpu_array.SetNumberOfValues(output.GetNumberOfCells())
  cpu_array.SetName("cosflow_cpu")

  for cellid in range(output.GetNumberOfCells()):
    cpu_array.SetValue(cellid, -1)

  list_of_files = glob.glob(casename+opts.multikey+'*.ele')
  ncpus = len(list_of_files)
  for file_name in list_of_files:
    if opts.verbose: print "  ", file_name
    cpu_no = int(file_name[len(casename + opts.multikey):-4])
    f = open(file_name, 'r')
    for line in f:
      tuple = line.split()
      if not tuple:
        continue
      cpu_array.SetValue(cellid_num[tuple[1]], cpu_no)
    f.close()

  output.GetCellData().AddArray(cpu_array)





# add the FieldData
command_used = vtk.vtkStringArray()
command_used.SetName("provenance")
command_used.InsertNextValue(" ".join(sys.argv))
output.GetFieldData().AddArray(command_used)

cn = vtk.vtkStringArray()
cn.SetName("cosflow_casename")
cn.InsertNextValue(casename)
output.GetFieldData().AddArray(cn)

if opts.multikey:
  mk = vtk.vtkStringArray()
  mk.SetName("cosflow_multikey")
  mk.InsertNextValue(opts.multikey)
  output.GetFieldData().AddArray(mk)
  ncpu = vtk.vtkIntArray()
  ncpu.SetName("cosflow_ncpus")
  ncpu.InsertNextValue(ncpus)
  output.GetFieldData().AddArray(ncpu)

output.GetFieldData().AddArray(mat_name)



# Write output file
if opts.verbose: print "Writing", out_file
writer = vtk.vtkXMLUnstructuredGridWriter()
writer.SetFileName(out_file)
if opts.ascii:
   writer.SetDataModeToAscii()
else:
   writer.SetDataModeToBinary()
writer.SetInputConnection(output.GetProducerPort())
writer.Write()
