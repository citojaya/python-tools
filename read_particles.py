#!/usr/bin/python


import numpy as np
import os
import sys
from optparse import OptionParser, OptionValueError

import numpy as np
#import scipy.linalg
import glob
import vtk
from vtk.util import numpy_support


# so we can find our ../lib no matter how we are called
findbin = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(findbin + '/../lib')



# parse command line
p = OptionParser(usage="""usage: %prog [options] <inputfile> <outputfile>
Reads binary file and written to ascii format


""")
p.add_option("-v", action="store_true",        dest="verbose",  help="Verbose")
p.add_option("-a", action="store_true", dest="ascii",  help="ASCII, instead of binary, .vtu output")

(opts, args) = p.parse_args()
# get the com filename
if len(args) != 2:
   p.print_help()
   sys.exit(1)
(in_file, out_file) = args


f = open(in_file+".dat", 'r')
#output = vtk.vtkPolyData()
output = vtk.vtkStructuredPoints()

nodes = []
nodes_array = vtk.vtkIntArray()
nodes_array.SetName("nodes")
node_x = []
node_y = []
node_z = []
counter = 0
point = vtk.vtkPoints()
field_data = vtk.vtkFieldData()

for line in f:
  if not line.split():
    continue

  tuple = ()
  tuple = line.split()

  nodes_array.InsertNextValue(counter)
  nodes.append(tuple[0])
  node_x.append(tuple[1])
  node_y.append(tuple[2])
  node_z.append(tuple[3])
  field_data.AddArray(nodes_array)
  counter += 1


  #print tuple


for i in range(0,counter):
  point.InsertNextPoint(float(node_x[i]),float(node_y[i]),float(node_z[i]))


#output.SetPoints(point)
#output.GetPointData().AddArray(nodes_array)
output.SetFieldData(field_data)

# Write output file
if opts.verbose: print "Writing", out_file+".poly"
writer = vtk.vtkStructuredPointsWriter()
writer.SetFileName(out_file+".vtp")


#if opts.ascii:
#   writer.SetDataModeToAscii()
#else:
#   writer.SetDataModeToBinary()

writer.SetInputConnection(output.GetProducerPort())
#writer.SetInputConnection(solid.GetProducerPort())
writer.Write()


print "DONE"

f.close()
