# ice-age-rivers
Computes drainage basin history and river discharge on a continental scale over a glacial cycle.

To run the code, use "run\_code.py". It's a bit clunky and needs to be done by hand, but usually works well. There are just a lot of files involved, and it takes some processing time!

Subfolders include:
* Model outputs for the paper ("...\_2016")
* Code to prepare gridded inputs of meltwater discharge to the coast for GCMs (GCM\_input)
* Code to compute ice-sheet volume in terms of sea-level equivalent (SLE\_out)
* Code to generate the figures for the ESurf paper (GRASSplot)
* Code to generate input files for this program (gridImport)
* Input files with time-steps (input)

### Reference:

**Wickert, A. D. (2016), [Reconstruction of North American drainage basins and river discharge since the Last Glacial Maximum](https://www.earth-surf-dynam.net/4/831/2016/esurf-4-831-2016.html), *Earth Surf. Dynam.*, *4*, 831-869, doi:10.5194/esurf-4-831-2016.**
