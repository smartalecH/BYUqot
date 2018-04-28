# ------------------------------------------------------------------ #
# objectLibrary.py
# ------------------------------------------------------------------ #
#
# A series of useful silicon photonic devices
#
# ------------------------------------------------------------------ #
# VERSION HISTORY
# 10 Apr 2018 - AMH - Initialization
# 24 Apr 2018 - AMH - Added vernier pattern, directional couplers
# 27 Apr 2018 - AMH - Added Bragg Grating
# ------------------------------------------------------------------ #
# List of components:
# ------------------------------------------------------------------ #
#
# S-Bend ................................................... 04/10/2018
# Y-Branch ................................................. 04/10/2018
# Branch Coupler ........................................... 04/24/2018
# Vernier Pattern .......................................... 04/24/2018
# Coupling Taper ........................................... 04/24/2018
# MZI ...................................................... 04/24/2018
# Bragg Grating ............................................ 04/27/2018

# ------------------------------------------------------------------ #
#      Import libraries
# ------------------------------------------------------------------ #

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)

import numpy
import gdspy
from scipy.interpolate import CubicSpline
import math

# ------------------------------------------------------------------ #
#      S-Bend
# ------------------------------------------------------------------ #
#
# Returns a cell with an s bend as described in the following paper:
# A. Kumar and S. Aditya, “Performance of S-bends for integrated-optic
# waveguides,” Microwave and Optical Technology Letters,
# vol. 19, no. 4, pp. 289–292, Nov. 1998.
# http://onlinelibrary.wiley.com/doi/10.1002/(SICI)1098-2760(199811)19:4%3C289::AID-MOP13%3E3.0.CO;2-Y/full
#
# Zs .................. Length of the S Bend (default = 10 microns)
# Xs .................. Height of the S Bend (default = 10 microns)
# waveguideWidth ...... Width of the waveguid (default = 0.5 microns)
# layerNumber    ...... Silicon layer number (default = 1)

# Returns: gds cell with the sbend structure

def SBend(Zs=10,Xs=10,waveguideWidth=0.5,layerNumber=1):

    z = numpy.linspace(0,Zs,1000);

    y = Xs / Zs * z - Xs / (2*numpy.pi) * numpy.sin(2*numpy.pi/Zs * z);

    points = numpy.stack([z,y],axis=1)

    poly_cell = gdspy.Cell('sBend_Length='+str(Zs)+'_Height='+str(Xs),exclude_from_current=True)
    poly_cell.add(gdspy.PolyPath(points,waveguideWidth,layer=layerNumber))

    return poly_cell

# ------------------------------------------------------------------ #
#      Y-Branch (from Lukas's paper)
# ------------------------------------------------------------------ #


# Returns a YBranch cell as described in Lukas's paper:
# https://www.osapublishing.org/oe/abstract.cfm?uri=oe-21-1-1310

# layerNumber .......... Silicon layer number (default = 1)

# Returns: gds cell with the Y Branch structure

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

# ------------------------------------------------------------------ #
#      Branch Coupler
# ------------------------------------------------------------------ #
# A simple broadband branch coupler as implemented in Lukas's toolbox.
# Paper references:
# Z. Lu et al., “Broadband silicon photonic directional coupler
# using asymmetric-waveguide based phase control,”
# Optics Express, vol. 23, no. 3, p. 3795, Feb. 2015.


def branchCoupler(layerNumber = 1, polarization = 'TE'):
    # Intialize cells
    symmetricCoupler = gdspy.Cell('symmetricCoupler',exclude_from_current=True)
    phaseControl     = gdspy.Cell('phaseControl',exclude_from_current=True)
    branchCoupler    = gdspy.Cell('branchCoupler_'+polarization,exclude_from_current=True)
    couplerTaper     = gdspy.Cell('couplerTaper',exclude_from_current=True)

    # Determine the curavture based on the mode of interest
    if polarization == 'TE':
        R = 5
        L1 = 12.4
        L2 = 4.6
    else:
        R = 10
        L1 = 2.2
        L2 = 6.1

    #  ------------  Assemble the symmetric coupler ----------------- #
    symmWaveguideWidth = 0.5;
    symmWaveguideSep   = 0.2;

    # Top waveguide
    point1 = (-L1/2,symmWaveguideSep/2)
    point2 = (L1/2,symmWaveguideSep/2 + symmWaveguideWidth)
    symmetricCoupler.add(gdspy.Rectangle(point1, point2, layer=layerNumber))

    # Bottom waveguide
    point1 = (-L1/2,-symmWaveguideSep/2)
    point2 = (L1/2,-(symmWaveguideSep/2 + symmWaveguideWidth))
    symmetricCoupler.add(gdspy.Rectangle(point1, point2, layer=layerNumber))

    #  ------------  Assemble the phaseControl ---------------------- #
    phaseControlTopWidth  = 0.6;
    phaseControlBotWidth  = 0.4;
    phaseControlTopSep = 0.1;
    phaseControlBotSep = 0.2;

    # Top waveguide
    point1 = (-L2/2,phaseControlTopSep)
    point2 = (L2/2,phaseControlTopSep + phaseControlTopWidth)
    phaseControl.add(gdspy.Rectangle(point1, point2, layer=layerNumber))

    # Bottom waveguide
    point1 = (-L2/2,-phaseControlBotSep)
    point2 = (L2/2,-(phaseControlBotSep + phaseControlBotWidth))
    phaseControl.add(gdspy.Rectangle(point1, point2, layer=layerNumber))

    #  ------------  Assemble the tapers ---------------------------- #
    taperLength = 1

    # top taper
    couplerTaper.add(gdspy.Polygon(
        [(-taperLength/2,symmWaveguideSep/2),
        (-taperLength/2,symmWaveguideSep/2 + symmWaveguideWidth),
        (taperLength/2,(phaseControlTopSep + phaseControlTopWidth)),
        (taperLength/2,(phaseControlTopSep))],
        layer=layerNumber
        ))

    # bottom taper
    couplerTaper.add(gdspy.Polygon(
        [(-taperLength/2,-symmWaveguideSep/2),
        (-taperLength/2,-(symmWaveguideSep/2 + symmWaveguideWidth)),
        (taperLength/2,-(phaseControlBotSep + phaseControlBotWidth)),
        (taperLength/2,-(phaseControlBotSep))],
        layer=layerNumber
        ))

    #  ------------  Assemble the entire structure ------------------ #
    symCouplerLeft = gdspy.CellReference(symmetricCoupler)
    symCouplerLeft.translate(-(L1/2 + L2/2 + taperLength),0)

    symCouplerRight = gdspy.CellReference(symmetricCoupler)
    symCouplerRight.translate((L1/2 + L2/2 + taperLength),0)

    taperLeft = gdspy.CellReference(couplerTaper)
    taperLeft.translate(-L2/2 - taperLength/2,0)

    taperRight = gdspy.CellReference(couplerTaper,x_reflection=True,rotation=180)
    taperRight.translate(L2/2 + taperLength/2,0)

    sBendCell = SBend(Zs=2*R,Xs=R/2,waveguideWidth=0.5,layerNumber=layerNumber)

    branchCoupler.add(gdspy.CellReference(phaseControl))
    branchCoupler.add(symCouplerLeft)
    branchCoupler.add(symCouplerRight)
    branchCoupler.add(taperLeft)
    branchCoupler.add(taperRight)
    branchCoupler.add(gdspy.CellReference(
        sBendCell).translate(L2/2 + L1 + taperLength,symmWaveguideSep/2 + symmWaveguideWidth/2))
    branchCoupler.add(gdspy.CellReference(
        sBendCell,x_reflection=True).translate(L2/2 + L1 + taperLength,-(symmWaveguideSep/2 + symmWaveguideWidth/2)))
    branchCoupler.add(gdspy.CellReference(
        sBendCell,x_reflection=True,rotation=180).translate(-(L2/2 + L1 + taperLength),symmWaveguideSep/2 + symmWaveguideWidth/2))
    branchCoupler.add(gdspy.CellReference(
        sBendCell,rotation=180).translate(-(L2/2 + L1 + taperLength),-(symmWaveguideSep/2 + symmWaveguideWidth/2)))

    branchCoupler.flatten();
    return branchCoupler

# ------------------------------------------------------------------ #
#      Vernier Pattern
# ------------------------------------------------------------------ #

def vernier(layerNumber = 1):

    # Vernier mask design parameters (all values in microns)
    numFingers        = 10       # Number of fingers to have on top and bottom
    fingerWidth       = 30       # Width of each finger
    fingerSpacing     = 40       # Spacing between fingers
    longFingerLength  = 200;     # Length of the long, middle finger
    shortFingerLength = 150;     # Length of the short, outer fingers
    baseThickness     = 76;      # Thickness of edge border of design

    separationDistance = 380     # distance from edge of pattern to origin


    # Intialize cells
    vernierCell   = gdspy.Cell('vernier')
    vernier2DCell = gdspy.Cell('vernier2D')

    # Calculate properties
    vernierWidth  = (longFingerLength + baseThickness)
    vernierHeight = (2*numFingers + 1) * fingerWidth + fingerSpacing * 2 * numFingers
    xCenter       = -(vernierWidth/2 - baseThickness)

    # First, place the long finger in the middle
    vernierCell.add(gdspy.Rectangle([xCenter, -fingerWidth/2],
    [xCenter+longFingerLength, fingerWidth/2],layer=layerNumber))

    # Next, iterate through and place the other fingers
    for k in range(1,numFingers+1):
        # Add top fingers
        vernierCell.add(gdspy.Rectangle(
        [xCenter, fingerWidth/2 + fingerSpacing*k + fingerWidth*(k-1)],
        [xCenter+shortFingerLength,fingerWidth/2 + fingerSpacing*k +fingerWidth*k],
        layer=layerNumber))

        # Add bottom fingers
        vernierCell.add(gdspy.Rectangle(
        [xCenter, -(fingerWidth/2 + fingerSpacing*k + fingerWidth*(k-1))],
        [xCenter+shortFingerLength,-(fingerWidth/2 + fingerSpacing*k +fingerWidth*k)],
        layer=layerNumber))

        # Finally, add the edge
        baseHeight = (2*numFingers + 1) * fingerWidth + fingerSpacing * 2 * numFingers
        vernierCell.add(gdspy.Rectangle([-vernierWidth/2, baseHeight/2],
        [xCenter, -baseHeight/2],layer=layerNumber))

    # Now let's flatten
    vernierCell.flatten()

    # Get vernier dimensions
    vernierDims   = vernierCell.get_bounding_box()
    vernierWidth  = abs(vernierDims[0,0] - vernierDims[1,0])
    vernierHeight = abs(vernierDims[0,1] - vernierDims[1,1])

    # Place one Vernier pattern in the x direction
    xCell = gdspy.CellReference(vernierCell,rotation=-90)
    xCell.translate(-(vernierHeight/2 + separationDistance),-vernierWidth/2)
    vernier2DCell.add(xCell)

    # Place another Vernier pattern in the y direction
    yCell = gdspy.CellReference(vernierCell,rotation=180)
    yCell.translate(-vernierWidth/2,-(vernierHeight/2 + separationDistance))
    vernier2DCell.add(yCell)

    # Return final cell
    return vernier2DCell

# ------------------------------------------------------------------ #
# Coupling Taper
# ------------------------------------------------------------------ #

def couplingTaper(layerNumber = 1, taperWidth = 0.180, waveguideWidth = 0.5, taperLength = 50):

    # Intialize cells
    couplingTaperCell   = gdspy.Cell('couplingTaper')
    pts = [(-taperLength/2, -taperWidth/2), (-taperLength/2, taperWidth/2),
           (taperLength/2, waveguideWidth/2),(taperLength/2, -waveguideWidth/2)]
    couplingTaperCell.add(gdspy.Polygon(pts,layer=layerNumber))

    return couplingTaperCell

# ------------------------------------------------------------------ #
# MZI
# ------------------------------------------------------------------ #

def MZI(deltaL = 40, Lref = 40, gapLength = 20,
        waveguideWidth = 0.5, bendRadius = 5,
        coupleType = "Y",polarization="TM",layerNumber = 1,
        maxH = 100):

    # Decide which coupler to use
    if coupleType == "Y":
        Coupler = YBranch()
    elif coupleType == "C":
        Coupler = branchCoupler(layerNumber = 1, polarization = polarization)
    else:
        raise Exception('Invalid Coupler Type Specified; must be Y or C')

    # Calculate how many, and how long each leg of the coupling arm needs to be
    numLegs = 2
    if deltaL/2 > maxH:
        numLegs = numpy.floor((deltaL) / maxH + 1)
        if numLegs % 2 == 1:
            numLegs = numLegs + 1
    leg = (deltaL/numLegs) + (4*bendRadius - numpy.pi*bendRadius)  # compensate with half circumf of circle

    # Connect bottom branch
    bottomBranch = gdspy.Cell('referenceBranch_'+(coupleType),exclude_from_current=True)
    bottomBranch.add(gdspy.Rectangle([-Lref/2 , waveguideWidth/2],[Lref/2, -waveguideWidth/2],layer=1))

    # Connect top branch
    gapLength = Lref / (numLegs+1)

    #outerArm = (Lref-gapLength)/2
    basicTurn = [1,-1,-1,1]
    basicLength = [gapLength,leg,gapLength,leg]
    turn   = []
    length = []
    for k in range(0,int(numLegs/2)):
        length.extend(basicLength)
        turn.extend(basicTurn)
    length.append(gapLength)
    topBranch = gdspy.Cell('deltaBranch_'+(coupleType),exclude_from_current=True)
    topBranchPoly = gdspy.L1Path((0, 0), '+x', waveguideWidth, length, turn,layer=1);
    topBranchPoly.fillet(radius=bendRadius)
    leftSquare    = gdspy.Rectangle([0,-waveguideWidth/2],[waveguideWidth,waveguideWidth/2],layer=1)
    rightSquare   = gdspy.Rectangle([Lref-waveguideWidth,-waveguideWidth/2],[Lref,waveguideWidth/2],layer=1)

    union = gdspy.fast_boolean(topBranchPoly,[leftSquare,rightSquare],'or',layer=layerNumber)

    topBranch.add(union)

    # Get Branch dimensions
    CouplerDims   = Coupler.get_bounding_box()
    CouplerWidth  = abs(CouplerDims[0,0] - CouplerDims[1,0])
    CouplerHeight = abs(CouplerDims[0,1] - CouplerDims[1,1])

    # Since the coupler isn't centered, let's compensate for the later translation
    if coupleType == "Y":
        CouplerWidth = 2*CouplerWidth

    MZIcell = gdspy.Cell('MZI_deltaL='+str(deltaL)+'_Lref='+str(Lref)+'_coupleType='+coupleType)
    MZIcell.add(gdspy.CellReference(Coupler, (-Lref/2 - CouplerWidth/2,0)))
    MZIcell.add(gdspy.CellReference(Coupler, (Lref/2 + CouplerWidth/2,0),rotation=180))
    MZIcell.add(gdspy.CellReference(bottomBranch, (0,-(CouplerHeight/2 - waveguideWidth/2.0))))
    MZIcell.add(gdspy.CellReference(topBranch, (-Lref/2,+(CouplerHeight/2 - waveguideWidth/2.0))))

    MZIcell.flatten()

    return MZIcell

# ------------------------------------------------------------------ #
#      Bragg Grating
# ------------------------------------------------------------------ #

def braggGrating(period = .310, NG = 400,
                 waveguideWidth = 0.5, dwidth = 0.05,layerNumber=1):
    # Intialize cells
    name = 'braggCell=' + str(int(period*1e3)) + '_NG='  + str(int(NG)) + '_dwidth=' + str(int(dwidth*1e3))
    braggCell = gdspy.Cell(name)

    # Calculate Parameters
    L = NG * period

    # Generate main waveguide
    waveguide = gdspy.Rectangle(
        [-L/2,waveguideWidth/2 -  dwidth/2],
        [L/2,-waveguideWidth/2 +  dwidth/2],
        layer=layerNumber)

    # Add side strips
    strips = []
    startRect = -L/2
    stopRect  = -L/2 + period/2
    for k in range(0,NG):
        strips.append(gdspy.Rectangle(
            [startRect,waveguideWidth/2 + dwidth/2],
            [stopRect,-waveguideWidth/2 - dwidth/2],
            layer=layerNumber))

        startRect = startRect + period
        stopRect = stopRect + period

    finalShape = gdspy.fast_boolean(waveguide,strips,'or',layer=layerNumber)
    braggCell.add(finalShape)
    return braggCell
