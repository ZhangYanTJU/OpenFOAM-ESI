#!/bin/sh
cd "${0%/*}" || exit                        # Run from this directory
#------------------------------------------------------------------------------

cp -f rawSurfaces/*.stl.gz .

# Vessel surface
surfaceAdd outlet.stl vessel.stl vessel.stl

# Sparger
surfaceCheck sparger.stl
surfaceAdd gasInlet.stl sparger_0.obj spargerInlet.stl
surfaceConvert sparger_1.obj spargerShaft.stl
surfaceOrient -inside spargerInlet.stl "(0 0.1 1)" spargerInlet.stl
surfaceOrient -inside spargerShaft.stl "(0 0.1 -0.1)" spargerShaft.stl

# Rotating shaft
surfaceOrient -inside shaftRotating.stl "(0 0.1 1)" shaftRotating.stl

# Static shaft
surfaceOrient -inside shaft.stl "(0 0.1 1)" shaft.stl

# Stirrer
surfaceSplitByTopology stirrer.stl stirrer.stl
surfaceOrient -inside stirrer.stl "(0 0.1 1)" stirrer.stl

mv -f stirrer_bafflePart_0.stl stirrer_baffles.stl

surfaceCheck stirrer_baffles.stl
mv -f stirrer_baffles_0.obj stirrer_baffles_plate.obj

#------------------------------------------------------------------------------
