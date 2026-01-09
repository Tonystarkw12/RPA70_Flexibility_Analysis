# PyMOL Script: Color RPA70 by Flexibility (B-factor)
# ===========================================================
# This script loads the 1JMC structure and colors residues
# from blue (rigid/low B-factor) to red (flexible/high B-factor)
#
# Usage: pymol color_by_flexibility.pml
# ===========================================================

# Load the structure
load 1JMC.pdb, rpa70

# Hide everything except protein
hide everything
show cartoon, rpa70

# Color by B-factor (spectrum coloring)
# Blue = low B-factor (rigid), Red = high B-factor (flexible)
spectrum b, blue_red, minimum=0, maximum=100, rpa70

# Set background
bg_color white

# Set representation
set cartoon_smooth_loops, 0
set cartoon_cylindrical_helices, 1
set cartoon_cylindrical_sheets, 1

# Show as cartoon with transparent surface
show surface, rpa70
set transparency, 0.5, rpa70

# Zoom to fit
zoom rpa70

# Add a color key (B-factor scale)
# Note: This creates a simple representation
# Blue  = Rigid (low B-factor, ~0-30 Å²)
# White = Medium (medium B-factor, ~30-60 Å²)
# Red   = Flexible (high B-factor, ~60+ Å²)

# Print statistics
print "=========================================="
print "RPA70 Flexibility Visualization"
print "=========================================="
print "Color scheme: B-factor (flexibility)"
print "  Blue  = Rigid (low B-factor)"
print "  White = Medium flexibility"
print "  Red   = Flexible (high B-factor)"
print "=========================================="
