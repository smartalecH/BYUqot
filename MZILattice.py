
import objectLibrary as obLib
import gdspy

# ------------------------------------------------------------------ #
#      Create single MZI
# ------------------------------------------------------------------ #


waveguideWidth = 0.5

L = 20
LTop = 30

legs = (LTop - L) / 2

YBranch = obLib.YBranch()

# Connect bottom branch
bottomBranch = gdspy.Cell('outputBranch')
bottomBranch.add(gdspy.Rectangle([-L/2 , waveguideWidth/2],[L/2, -waveguideWidth/2],layer=1))

# Do top branch





cen2Branch = 2 + 1.2/2 - 5/2.0

MZI = gdspy.Cell('MZI')
MZI.add(gdspy.CellReference(YBranch, (-L/2 - 15,0)))
MZI.add(gdspy.CellReference(YBranch, (L/2 + 15,0),rotation=180))
MZI.add(gdspy.CellReference(bottomBranch, (0,-(cen2Branch))))



# ------------------------------------------------------------------ #
#      OUTPUT
# ------------------------------------------------------------------ #

# Output the layout to a GDSII file (default to all created cells).
# Set the units we used to micrometers and the precision to nanometers.
gdspy.write_gds('tutorial.gds', unit=1.0e-6, precision=1.0e-9)
