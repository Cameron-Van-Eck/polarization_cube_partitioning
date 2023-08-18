# polarization_cube_partitioning
Some small scripts for partitioning polarization cubes into small pieces. 
Not intended to be particularly general.

I wrote these for POSSUM cubes, as a quick and dirty way to divide them
into smaller pieces to improve the disk IO speed of the pieplines.

In the initial version, is splits FITS cubes into 4 separate pieces
(splitting in half vertically and horizontally).

The partition_fits_cube.py script divides a set of Stokes I, Q, U cubes
into 4 parts, as well as a Selavy source list.

The merge_quarters script assembles the 1D pipeline outputs 
(catalog, spectra table, FDF table)
back into one file. It does not currently reassemble the FITS files!
