from __future__ import (absolute_import, division,
                        print_function, unicode_literals)

import numpy
import gdspy
from scipy.interpolate import CubicSpline
import math

# ------------------------------------------------------------------ #
#      S Bend
# ------------------------------------------------------------------ #
'''
Returns a cell with an s bend as described in the following paper:
http://onlinelibrary.wiley.com/doi/10.1002/(SICI)1098-2760(199811)19:4%3C289::AID-MOP13%3E3.0.CO;2-Y/full

Zs .................. Length of the S Bend
Xs .................. Height of the S Bend
waveguideWidth ...... Width of the waveguid
layerNumber    ...... Silicon layer number

'''
'''
class SBend(gdspy.Cell):

    def __init__(self,Zs=10,Xs=10,waveguideWidth=0.5,layerNumber=1):
        # Initialize super class
        super(gdspy.Cell, self).__init__()

        self.Zs = Zs
        self.Xs = Xs
        self.waveguideWidth = waveguideWidth
        self.layerNumber = layerNumber

        z = numpy.linspace(0,Zs,1000);
        y = self.Xs / self.Zs * z - self.Xs / (2*numpy.pi) * numpy.sin(2*numpy.pi/self.Zs * z);
        points = numpy.stack([z,y],axis=1)

        self.name = 'sBend_Length='+str(self.Zs)+'_Height='+str(self.Xs)
        poly = gdspy.PolyPath(points,self.waveguideWidth,layer=self.layerNumber)

        super(gdspy.Cell, self).add(poly)

'''
def SBend(Zs=10,Xs=10,waveguideWidth=0.5,layerNumber=1):

    z = numpy.linspace(0,Zs,1000);

    y = Xs / Zs * z - Xs / (2*numpy.pi) * numpy.sin(2*numpy.pi/Zs * z);

    points = numpy.stack([z,y],axis=1)

    poly_cell = gdspy.Cell('sBend_Length='+str(Zs)+'_Height='+str(Xs))
    poly_cell.add(gdspy.PolyPath(points,waveguideWidth,layer=layerNumber))

    return poly_cell

# ------------------------------------------------------------------ #
#      Y-Branch (from Lukas's paper)
# ------------------------------------------------------------------ #

'''
Returns a YBranch cell as described in Lukas's paper:
https://www.osapublishing.org/oe/abstract.cfm?uri=oe-21-1-1310

layerNumber    ...... Silicon layer number

'''

def YBranch(layerNumber=1):

    # X location of spline points
    taperLength = 2.0;
    x = numpy.linspace(0,taperLength,13)

    # Width of waveguides
    waveguideWidth = 0.5

    # Widths of each spline length as described in the paper
    w = numpy.array([0.5, 0.5, 0.6, 0.7, 0.9, 1.26, 1.4, 1.4, 1.4, 1.4, 1.31, 1.2, 1.2]);

    # Top and bottom vectors for splines
    topY   = w/taperLength

    # Spline measurements
    numSpline    = 90
    cs           = CubicSpline(x, topY)
    xs           = numpy.linspace(0,taperLength,numSpline)
    topSpline    = cs(xs)
    bottomSpline = -topSpline[::-1]
    xx    = numpy.concatenate([xs,xs[::-1]])
    yy    = numpy.concatenate([topSpline,bottomSpline])

    vertices = numpy.stack([xx,yy],axis=1)

    # Create the taper polygon and cell
    yTaperCell = gdspy.Cell('yTaper')
    taperPoly = gdspy.Polygon(vertices, layerNumber)
    taperPoly.fracture()
    yTaperCell.add(taperPoly)

    # S Bends
    Zs = 13
    Xs = 2

    sBendCellTop = SBend(Zs,Xs,waveguideWidth,layerNumber)
    translationTop = [taperLength,topY[-1] - waveguideWidth/2.0]
    translationBottom = [taperLength,-(topY[-1] - waveguideWidth/2.0)]

    # Output branch (single waveguide)
    #L = 1.5
    #outputBranchCell = gdspy.Cell('outputBranch')
    #outputBranchCell.add(gdspy.PolyPath([[0, 0],[L, 0]],waveguideWidth,layer=layerNumber))

    # First we need a cell to add the polygons to.
    yBranchCell = gdspy.Cell('Ybranch')

    # Add all cells to master cell
    yBranchCell.add(gdspy.CellReference(yTaperCell, (0, 0), x_reflection=False))
    yBranchCell.add(gdspy.CellReference(sBendCellTop, translationTop, x_reflection=False))
    yBranchCell.add(gdspy.CellReference(sBendCellTop, translationBottom, x_reflection=True))
    #yBranchCell.add(gdspy.CellReference(outputBranchCell, (-L,0), x_reflection=False))

    return yBranchCell
