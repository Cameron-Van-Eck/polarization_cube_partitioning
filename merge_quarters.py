import numpy as np
import astropy.io.fits as pf
import rmtable
import polspectra
import glob
import astropy.table as at

import sys

SB=sys.argv[-1]

cat_files=glob.glob(f'pipeline/rms1d/possum/*{SB}*.csv')
tabs=[]
for file in cat_files:
    tabs.append(at.Table.read(file))
alltab=at.vstack(tabs,join_type='exact')
alltab.write(f'/home/cameron/SB{SB}_pipeline.csv')

spectra_files=glob.glob(f'pipeline/rms1d/possum/*{SB}*_polspectra.fits')
spectra=polspectra.from_FITS(spectra_files[0])
for file in spectra_files[1:]:
    spectra.merge_tables(polspectra.from_FITS(file),merge_type='exact',source_numbers='concat')
spectra.write_FITS(f'/home/cameron/SB{SB}_pipeline_polspectra.fits')


FDF_files=glob.glob(f'pipeline/rms1d/possum/*{SB}*_FDFtable.fits')
FDFs=[]
for file in FDF_files:
    FDFs.append(pf.getdata(file))
all_FDFs=at.vstack([ at.Table(x) for x in FDFs],join_type='exact')
all_FDFs.write(f'/home/cameron/SB{SB}_pipeline_FDFtable.fits')
