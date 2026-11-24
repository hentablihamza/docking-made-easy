from pymol.cgo import *
from pymol import cmd
minX, minY, minZ = 9.697, -6.072, 44.405
maxX, maxY, maxZ = 33.697, 6.678, 59.780
box = [LINEWIDTH, 2.0, BEGIN, LINES, COLOR, 1.0, 1.0, 0.0]
corners = [(minX, minY, minZ), (maxX, minY, minZ), (maxX, maxY, minZ), (minX, maxY, minZ),
           (minX, minY, maxZ), (maxX, minY, maxZ), (maxX, maxY, maxZ), (minX, maxY, maxZ)]
edges = [(0,1),(1,2),(2,3),(3,0),(4,5),(5,6),(6,7),(7,4),(0,4),(1,5),(2,6),(3,7)]
for i, j in edges:
    box += [VERTEX, *corners[i], VERTEX, *corners[j]]
box.append(END)
cmd.load_cgo(box, 'docking_box')
