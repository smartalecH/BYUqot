
import objectLibrary as obLib
import gdspy

# ------------------------------------------------------------------ #
#      Create single MZI
# ------------------------------------------------------------------ #

def MZI(deltaL = 40, Lref = 20, gapLength = 10, waveguideWidth = 0.5, bendRadius = 5):

    layerNumber = 1

    LTop = Lref + deltaL

    leg = (LTop - Lref) / 2

    YBranch = obLib.YBranch()

    # Connect bottom branch
    bottomBranch = gdspy.Cell('referenceBranch')
    bottomBranch.add(gdspy.Rectangle([-Lref/2 , waveguideWidth/2],[Lref/2, -waveguideWidth/2],layer=1))

    # Do top branch
    bendRadius = 5;

    outerArm = (Lref-gapLength)/2

    length = [outerArm,leg,gapLength,leg,outerArm]
    print(length)
    turn = [1,-1,-1,1]
    topBranch = gdspy.Cell('deltaBranch')

    topBranchPoly = gdspy.L1Path((0, 0), '+x', waveguideWidth, length, turn,layer=1);
    topBranchPoly.fillet(bendRadius)
    leftSquare    = gdspy.Rectangle([0,-waveguideWidth/2],[waveguideWidth,waveguideWidth/2],layer=1)
    rightSquare   = gdspy.Rectangle([Lref-waveguideWidth,-waveguideWidth/2],[Lref,waveguideWidth/2],layer=1)

    union = gdspy.fast_boolean(topBranchPoly,[leftSquare,rightSquare],'or',layer=layerNumber)

    topBranch.add(union)


    # Half the hieght of the entire splitter - used to align branches with splitters
    cen2Branch = 2.0 + 1.2/2 - waveguideWidth/2.0

    MZIcell = gdspy.Cell('MZI')
    MZIcell.add(gdspy.CellReference(YBranch, (-Lref/2 - 15,0)))
    MZIcell.add(gdspy.CellReference(YBranch, (Lref/2 + 15,0),rotation=180))
    MZIcell.add(gdspy.CellReference(bottomBranch, (0,-(cen2Branch))))
    MZIcell.add(gdspy.CellReference(topBranch, (-Lref/2,+(cen2Branch))))

    return MZIcell



# ------------------------------------------------------------------ #
#      OUTPUT
# ------------------------------------------------------------------ #

MZIcell = MZI(deltaL = 30, Lref = 20, waveguideWidth = 0.5);

# Output the layout to a GDSII file (default to all created cells).
# Set the units we used to micrometers and the precision to nanometers.
gdspy.write_gds('tutorial.gds', unit=1.0e-6, precision=1.0e-9)
