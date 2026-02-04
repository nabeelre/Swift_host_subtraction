import argparse
import os
import shutil
import SwiftPhotom.help as SH
import SwiftPhotom.uvot as up


def run_swift_photom(infile, sn_reg='sn.reg', bg_reg='snbkg.reg',
                     det_limit=3.0, ab=False, filter='ALL',
                     no_combine=False, obj='test'):
    """
    Run Swift photometry analysis on science and template images.

    Parameters
    ----------
    infile : list of str
        List of input file paths. First file(s) are science images,
        remaining are template images.
    sn_reg : str, optional
        Source region file path (default: 'sn.reg')
    bg_reg : str, optional
        Background region file path (default: 'snbkg.reg')
    det_limit : float, optional
        Detection limit in sigma (default: 3.0)
    ab : bool, optional
        Use AB magnitudes (default: False)
    filter : str, optional
        Filter to process, or 'ALL' for all filters (default: 'ALL')
    no_combine : bool, optional
        Do not combine images (default: False)
    obj : str, optional
        Object name for output (default: 'test')

    Returns
    -------
    dict
        Dictionary containing photometry results with aperture sizes as keys
    """
    ap_size = up.get_aperture_size(sn_reg)
    user_ap = ap_size + '_arcsec'

    obj_file_list, tem_file_list = up.interpret_infile(infile)

    obj_file_list = up.sort_file_list(obj_file_list)
    tem_file_list = up.sort_file_list(tem_file_list)

    if os.path.isdir('reduction'):
        shutil.rmtree('reduction')

    os.mkdir('reduction')

    mag = {user_ap: [], '5_arcsec': []}

    filt_list = up.sort_filters(filter)

    for filter_name in filt_list:

        if filter_name not in obj_file_list:
            continue

        filter_dir = os.path.join('reduction', filter_name)
        if not os.path.isdir(filter_dir):
            os.mkdir(filter_dir)

        print('Working on filter ' + filter_name)
        template_exists = 1

        print('Creating product file for the object.')
        prod_file = up.create_product(obj_file_list[filter_name], filter_name,
                                      no_combine=no_combine)

        print('Running uvotmaghist on the object image.\n')
        phot_file = up.run_uvotmaghist(prod_file, sn_reg, bg_reg, filter_name)

        if filter_name not in tem_file_list:
            print('No template provided for filter ' + filter_name +
                  '. Simple aperture photometry will be performed.\n\n')

            template_exists = 0

        if template_exists:
            print('Creating product file for the template.')
            prod_file = up.create_product(tem_file_list[filter_name], filter_name,
                                          template=1, no_combine=no_combine)

            print('Running uvotmaghist on the template image.\n')
            templ_file = up.run_uvotmaghist(prod_file, sn_reg, bg_reg,
                                            filter_name)

            mag_filter = up.extract_photometry(phot_file, ab, det_limit,
                                               ap_size, templ_file)

        else:
            mag_filter = up.extract_photometry(phot_file, ab, det_limit,
                                               ap_size)

        for ap in mag_filter:
            mag[ap] = mag[ap] + mag_filter[ap]

        print('\n')

    up.output_mags(mag, ap_size, obj=obj)

    return mag


def parse_args():
    """
    Parse command-line arguments for Swift photometry.

    Returns
    -------
    argparse.Namespace
        Parsed command-line arguments
    """
    parser = argparse.ArgumentParser(
        description=SH.description,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("infile", nargs='+', help=SH.infile_help)
    parser.add_argument("-s", "--sn", dest="sn_reg", type=str, default='sn.reg',
                        help=SH.sn_reg)
    parser.add_argument("-b", "--bg", dest="bg_reg", type=str, default='snbkg.reg',
                        help=SH.bg_reg)
    parser.add_argument("-d", "--detection", dest="det_limit", type=float, default=3.,
                        help=SH.det_limit)
    parser.add_argument("-a", "--ab", dest="ab", default=0, action='store_true',
                        help=SH.ab_mag)
    parser.add_argument("-f", "--filter", dest="filter", default='ALL', help=SH.filter)
    parser.add_argument("--no_combine", dest="no_combine", default=0, action='store_true',
                        help=SH.no_combine)
    parser.add_argument("--obj", dest="obj", default='test', help=SH.obj)

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    run_swift_photom(
        infile=args.infile,
        sn_reg=args.sn_reg,
        bg_reg=args.bg_reg,
        det_limit=args.det_limit,
        ab=bool(args.ab),
        filter=args.filter,
        no_combine=bool(args.no_combine),
        obj=args.obj
    )
