==Global30arcsec==

The purpose of this code is to import sea-level model outputs, ice-sheet reconstructions, and paleoclimate model outputs, into a GRASS GIS location.

This folder does not include Jerry Mitrovica's Fortran code to convert spherical harmonic model outputs onto a grid, or any of the model outputs themselves. It needs these in order to work.

Global30arcsec was the original folder name, because I hard-coded the Fortran outputs to provide global grids at 30 arcseconds meridionally (I think) and then interpolated to 30 arcseconds zonally. The interpolation is valid because the spherical harmonics were extracted to a high enough order that everything finer is flexurally smoothed. I also use these codes to convert the model output topographies into anomolies from present-day topography.
