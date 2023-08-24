#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  3 10:06:00 2021
@author: cvaneck

Breaks large cube into 4 parts.

Call sequence:

script.py incube source_list outfile_basename

Version 1.1.0 2023-08-24
 - Uses argparse to handle command line arguments
 - previous hard-coded quartering and naming bound into a loop, allowing different partitioning
 - size of partitions controlled by n_div; n_div = 2 gives original behaviour
 - NOTE: n_div > 2 untested *****
 - partition_fits_cube.py -h gives help
"""

import astropy.io.fits as pf
import sys
import astropy.wcs as wcs
import astropy.table as at
import argparse as ap

__version__ = '1.1.0'


def file_name_augment(filename, incr):
    # Filename of form a.b.ext is returned as a.b.incr.ext
    # The routine adds '.' separators to incr

    parts = filename.split('.')
    parts_a = parts[:-1] + [incr] + parts[-1:]
    return '.'.join(parts_a)


def arg_init():
    parser = ap.ArgumentParser(prog='partition_fits_cube', formatter_class=ap.RawDescriptionHelpFormatter,
                               description='Divide I,Q,U cubes into partitions.',
                               epilog='See -x for more explanation')
    # noinspection PyTypeChecker
    parser.add_argument('incube', nargs=1, type=str, help="Input cube name")
    parser.add_argument('-c', '--cat', type=str, help="Catalogue name")
    parser.add_argument('-o', '--out_dir', type=str, help='Output directory')

    parser.add_argument('-n', '--n_div', type=int, default=2, help="Number of partitions per axis [%(default)d]")
    parser.add_argument('-d', '--do_num', type=str, default='all', help=f"Number of partitions to do [%(default)s]")
    parser.add_argument('-p', '--padding', type=int, default=38,
                        help="Width of overlap (padding) between partitions [%(default)d]")
    parser.add_argument('-v', '--verbose', action='store_true')
    parser.add_argument('-x', '--explain', action='store_true', help="give an expanded explanation")

    return parser


def main():
    args = arg_init().parse_args()

    basefile = args.incube[0].replace('.i.', '._.')
    outfile_base_dir = args.out_dir
    cat_file = args.cat
    n_div = args.n_div
    do_num = n_div ** 2
    if args.do_num.isnumeric():
        do_num = int(args.do_num)

    padding = args.padding
    pd = padding

    # basefile = basefile.replace('.i.','._.').replace('.q.','._.').replace('.u.','._.')
    fits_in = {s: basefile.replace('._.', f'.{s}.') for s in ['i', 'q', 'u']}

    # Open I cube, get size:
    hdu = pf.open(fits_in['i'], memmap=True)

    header = hdu[0].header
    xsize = header['NAXIS1']
    ysize = header['NAXIS2']

    hdu.close()

    print(f'Input image {xsize:d} x {ysize:d}')

    print('Prepare for catalogue division')

    cat = at.Table.read(cat_file)
    extension = cat_file[-4:]
    if extension == '.csv':
        fileformat = 'csv'
    elif extension == '.xml':
        fileformat = 'votable'
    else:
        raise Exception('This script currently only supports CSV or VOTable catalogs.')

    cs = wcs.WCS(header, fix=False, naxis=(1, 2))
    # Compute the catalogue coordinates to pixel positions in the image
    dx, dy = cs.all_world2pix(cat['col_ra_deg_cont'], cat['col_dec_deg_cont'], 0)

    bx, by = [], []
    sx, sy = [], []

    for i in range(n_div + 1):
        bx.append(i * xsize // n_div)
        by.append(i * ysize // n_div)

    for i in range(n_div):
        sx.append(slice(max(0, bx[i] - pd), min(xsize, bx[i + 1] + pd)))
        sy.append(slice(max(0, by[i] - pd), min(ysize, by[i + 1] + pd)))

    for ir in range(n_div):
        for ic in range(n_div):
            rc = f'R{ir:d}C{ic:d}'
            cat_out = file_name_augment(cat_file, rc)
            # cat_out = f'{outfile_base}.{rc}{extension}'
            sub_cat = (bx[ic] <= dx) & (dx <= bx[ic + 1]) & (by[ir] <= dy) & (dy <= by[ir + 1])
            at.Table(cat[sub_cat]).write(cat_out, format=fileformat, overwrite=True)
            print(rc, cat_out, sub_cat.sum())

    if do_num > 0:
        n_done = 0
        hdu_d = {}
        for stokes in ['i', 'q', 'u']:
            hdu = pf.open(fits_in[stokes], memmap=True)
            hdu_d[stokes] = hdu

        for ir in range(n_div):
            for ic in range(n_div):
                rc = f'R{ir:d}C{ic:d}'
                for stokes in ['i', 'q', 'u']:
                    file_base = outfile_base_dir + '/' + fits_in[stokes].split('/')[-1]
                    fits_out = file_name_augment(file_base, rc)

                    subimage = hdu_d[stokes][0].data[:, :, sy[ir], sx[ic]]

                    outheader = hdu_d[stokes][0].header.copy()
                    outheader['CRPIX1'] = hdu_d[stokes][0].header['CRPIX1'] - sx[ic].start
                    outheader['CRPIX2'] = hdu_d[stokes][0].header['CRPIX2'] - sy[ir].start
                    prim = pf.PrimaryHDU(header=outheader, data=subimage)
                    out = pf.HDUList([prim] + hdu_d[stokes][1:])
                    print(stokes, rc, fits_out, sy[ir], sx[ic])
                    out.writeto(fits_out)
                n_done += 1
                if n_done == do_num:
                    break
            if n_done == do_num:
                break

        for stokes in ['i', 'q', 'u']:
            hdu_d[stokes].close()


if __name__ == "__main__":
    sys.exit(main())
