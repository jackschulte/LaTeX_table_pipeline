"""Generates the table of stellar properties measured in the literature."""

import numpy as np
import os
import pandas as pd
from astroquery.vizier import Vizier
from astropy.coordinates import Angle
from grab_tres_vsini import grab_tres_vsini
import re
import logging
from table_utils import round_sig_figs, write

def gen_lit_str(array, value, error=None):
    '''
    Appends a string version of the literature value obtained via query to the appropriate list of strings

    Parameters
    -----------
    array: the list/array to append the value string to
    value: the value to append to the array
    error: the error on that value, if there is one
    '''

    if (value == None) or (isinstance(value, np.ma.core.MaskedConstant)):
        array.append('& --- ') # if the value doesn't exist in the query, add in the filler
    elif (error == None) or (isinstance(error, np.ma.core.MaskedConstant)):
        array.append(r'& ' + str(value) + ' ') # if there's no error, only append the value
    else:
        array.append(r'& $' + str(value) + r' \pm ' + str(error) + '$ ') # if there's a value and associated error, append them both

def add_source(array, source):
    '''
    Appends a number corresponding to one of the sources to the final column of the table.

    Parameters
    -----------
    array: the array corresponding to the row that you are adding the source to
    source: an integer or string corresponding to the source that you wish to append
    '''

    array.append(f' & {source} ')

def lit_table(target_list, path, file_prefix=None, outputpath='.', vsini_type='gaia', vsini_external=None, tres_username=None, tres_password=None,
              add_source_column=False, grab_mags_from_sedfile=True, max_targets_per_table=5, MNRAS=False):
    '''
    Generates a 'literature' table, using photometric and astrometric parameters from Gaia, 2MASS, and WISE. Optionally
    grabs vsini measurements from TRES. WARNING: Collecting TRES vsini measurements will increase runtime by ~4 min.

    Parameters
    -----------
    target_list: an array of strings containing the names of each target. Nominally, these should be TOI IDs. Ex: ['TOI-1855', 'TOI-2107']
    path: path of EXOFASTv2 output files
    file_prefix: an array of strings containing the file prefix used in each EXOFASTv2 fit. Only necessary if grab_mags_from_sedfile=True
    outputpath: the folder in which the table should be generated. Current working directory by default
    vsini_type: Accepts 'gaia' to use Gaia's vbroad, 'tres' to scrape the TRES/CHIRON site (tess.exoplanets.dk), or 'external' to provide an external array of vsini values
    vsini_external: N by 2 array containing vsini (first column) and vsini_err (second column) from an external source
    tres_username: String including the user's TRES website username if vsini_type='tres'
    tres_passworde: String including the user's TRES website password if vsini_type='tres'
    add_source_column: a boolean to determine whether a source column is added.
    grab_mags_from_sedfile: a boolean to determine whether the magnitudes are collected from an EXOFAST SED file instead of astroquery
    max_targets_per_table: the maximum number of targets shown in a single table. If len(target_list) exceeds this value, the
        targets are split across multiple lit_table.tex files. Every table after the first is captioned "\\textit{(Continued)}"
        and every table except the last gets "\\addtocounter{table}{-1}" so that all pieces share one table number. The
        \\begin{minipage} notes block is only written in the last table. Set to None (or 0) to force a single table.
    MNRAS: set to True to use the MNRAS class's \\contcaption for the continuation tables, in place of the
        "\\caption{\\textit{(Continued)}}" and "\\addtocounter{table}{-1}" pair written otherwise. \\contcaption
        supplies the continuation wording and holds the table number itself, so neither is needed alongside it.
    '''

    # Setting up to save the table as a .tex file

    if os.path.exists(outputpath) == False:
        os.mkdir(outputpath)

    # Turning TOIs into TIC IDs

    url="https://exofop.ipac.caltech.edu/tess/download_toi.php?sort=toi&output=pipe"
    TOI_df=pd.read_csv(url, delimiter='|', index_col=1)
    TIC_IDs = [] # initializing list
    TESS_mags = [] # TESS mags are in exofop, not Vizier
    TESS_mags_err = []
    for toi in target_list:
        toi_id = float(toi[4:]) + 0.01
        TIC_ID = TOI_df.loc[toi_id]['TIC ID']
        TESS_mag = TOI_df.loc[toi_id]['TESS Mag']
        TESS_mag_err = TOI_df.loc[toi_id]['TESS Mag err']
        TIC_IDs.append(TIC_ID)
        TESS_mags.append(TESS_mag)
        TESS_mags_err.append(TESS_mag_err)
    
    # grabbing vsini from the TRES/CHIRON site

    vsini_tres = []
    vsini_tres_err = []
    if vsini_type == 'tres':
        for ticid in TIC_IDs:
            vsini, vsini_err = grab_tres_vsini(tres_username, tres_password, ticid)
            if vsini and not np.isnan(vsini):
                vsini_3sigfig = round_sig_figs(vsini, 3)
                vsini_tres.append(vsini_3sigfig)
                if '.' in str(vsini_3sigfig):
                    decimal_places = len(str(vsini_3sigfig).split('.')[1])
                else:
                    decimal_places = 0
                vsini_tres_err.append(round(vsini_err, decimal_places)) # round the vsini error to the same number of decimal places as the vsini
                
            else:
                vsini_tres.append(None)
                vsini_tres_err.append(None)

    # initializing rows
    ra_arr=[r'$\alpha_{J2000}\ddagger$ & Right Ascension (h:m:s) ']
    dec_arr=[r'$\delta_{J2000}\ddagger$ & Declination (d:m:s) ']
    '''
    \ddagger for the following note: Right Ascension and Declination are in epoch J2000. 
    The coordinates come from Vizier where the Gaia RA and Dec have been precessed and corrected to J2000 from epoch J2016.
    '''
    gaia_g_arr=[r'${\rm G}$ & Gaia $G$ mag. ']
    gaia_bp_arr=[r'$G_{\rm BP}$ & Gaia $G_{\rm BP}$ mag. ']
    gaia_rp_arr=[r'$G_{\rm RP}$ & Gaia $G_{\rm RP}$ mag. ']
    tmag_arr=[r'${\rm T}$ & TESS mag. ']
    j_2mass_arr=[r'$J$ & 2MASS $J$ mag. ']
    h_2mass_arr=[r'$H$ & 2MASS $H$ mag. ']
    k_2mass_arr=[r'$K$ & 2MASS $K$ mag. ']
    wise1_arr=[r'$W1$ & WISE $W1$ mag. ']
    wise2_arr=[r'$W2$ & WISE $W2$ mag. ']
    wise3_arr=[r'$W3$ & WISE $W3$ mag. ']
    wise4_arr=[r'$W4$ & WISE $W4$ mag. ']
    pmra_arr=[r'$\mu_{\alpha}$ & Gaia DR3 proper motion in RA (mas yr$^{-1}$)']
    pmdec_arr=[r'$\mu_{\delta}$ & Gaia DR3 proper motion in Dec (mas yr$^{-1}$)']
    parallax_arr=[r'$\pi$ & Gaia DR3 Parallax (mas) ']
    vsini_arr=[r'$v\sin{i_\star}$ & Projected rotational velocity (km s$^{-1}$) ']
    # add note that vbroad includes other effects such as macroturbulence, template mismatch, and other instrumental effects

    # initializing additional identifier arrays (one ' & ...' entry per target so they can be sliced per table)
    tic_id_list=[]
    tycho_id_list=[]
    twomass_id_list=[]
    gaia_id_list=[]

    wise4_targets = 0 # how many targets have a WISE4 magnitude; the row is dropped if none do

    # query for data and build rows
    for i in range(len(TIC_IDs)):
        # querying for data
        gaia_columns=['_r', 'RAJ2000', 'DEJ2000', 'Gmag', 'e_Gmag', 'BPmag', 'e_BPmag', 'RPmag', 'e_RPmag', 'pmRA', \
            'e_pmRA', 'pmDE', 'e_pmDE', 'Plx', 'e_Plx', 'Vbroad', 'e_Vbroad', 'TYC2', '2MASS', 'DR3Name']
        vgaia = Vizier(columns=gaia_columns, catalog='I/355/gaiadr3')
        data_gaia = vgaia.query_region('TIC ' + str(TIC_IDs[i]), radius=Angle(6, "arcsec"))
        data_gaia = data_gaia[0] # Vizier returns a list of tables, but we only want the first one
        data_gaia.sort('_r') # sort by distance from the target star
        if len(data_gaia) > 1:
            logging.warning(f'Multiple sources in Vizier Gaia query of TIC {TIC_IDs[i]}. Selecting closest source.')
        data_gaia = data_gaia[0]

        twomass_columns=['_r', 'Jmag', 'e_Jmag', 'Hmag', 'e_Hmag', 'Kmag', 'e_Kmag']
        v2mass = Vizier(columns=twomass_columns, catalog='II/246/out')
        data_2MASS = v2mass.query_region('TIC ' + str(TIC_IDs[i]), radius=Angle(6, "arcsec"))
        data_2MASS = data_2MASS[0] # Vizier returns a list of tables, but we only want the first one
        data_2MASS.sort('_r') # sort by distance from the target star
        if len(data_2MASS) > 1:
            logging.warning(f'Multiple sources in Vizier 2MASS query of TIC {TIC_IDs[i]}. Selecting closest source.')
        data_2MASS = data_2MASS[0]
        

        WISE_columns=['_r', 'W1mag', 'e_W1mag', 'W2mag', 'e_W2mag', 'W3mag', 'e_W3mag', 'W4mag', 'e_W4mag']
        vwise = Vizier(columns=WISE_columns, catalog='II/311/wise')
        data_WISE = vwise.query_region('TIC ' + str(TIC_IDs[i]), radius=Angle(6, "arcsec"))
        data_WISE = data_WISE[0] # Vizier returns a list of tables, but we only want the first one
        data_WISE.sort('_r') # sort by distance from the target star
        if len(data_WISE) > 1:
            logging.warning(f'Multiple sources in Vizier WISE query of TIC {TIC_IDs[i]}. Selecting closest source.')
        data_WISE = data_WISE[0]

        # store other identifiers
        tic_id_list.append(' & TIC ' + str(TIC_IDs[i]))
        tycho_id = str(data_gaia['TYC2'])
        if len(tycho_id) > 0:
            tycho_id_list.append(' & TYC ' + tycho_id)
        else:
            tycho_id_list.append(' & ---')
        twomass_id = str(data_gaia['2MASS'])
        if len(twomass_id) > 0:
            twomass_id_list.append(' & J' + twomass_id)
        else:
            twomass_id_list.append(' & ---')
        gaia_id = str(data_gaia['DR3Name'])
        gaia_id = re.sub('Gaia DR3 ', '', gaia_id) # removing prefix
        gaia_id_list.append(f' & {gaia_id}')

        # Grabbing and formatting RA/Dec
        ra = data_gaia['RAJ2000']
        ra_angle = Angle(ra, 'deg')
        ra_hr = int(ra_angle.hms[0])
        if ra_hr < 10: # filling with zeroes to match hh:mm:ss format
            ra_hr = f'0{ra_hr}'
        ra_min = int(ra_angle.hms[1])
        if ra_min < 10:
            ra_min = f'0{ra_min}'
        ra_sec = round(float(ra_angle.hms[2]), 3)
        if ra_sec < 10:
            ra_sec = f'0{ra_sec}'
        ra_str = f'{ra_hr}:{ra_min}:{ra_sec}'
        dec = data_gaia['DEJ2000']
        dec_angle = Angle(dec, 'deg')
        dec_deg = int(dec_angle.dms[0])
        if abs(dec_deg) < 10:
            dec_deg = f'0{dec_deg}'
        dec_min = abs(int(dec_angle.dms[1]))
        if dec_min < 10:
            dec_min = f'0{dec_min}'
        dec_sec = abs(round(float(dec_angle.dms[2]), 3))
        if dec_sec < 10:
            dec_sec = f'0{dec_sec}'
        dec_str = f'{dec_deg}:{dec_min}:{dec_sec}'

        # Grabbing the used magnitudes from SED files
        if grab_mags_from_sedfile == True:
            columns = ['bandname', 'magnitude', 'used_errors', 'catalog_errors', 'star_index']
            sedtable = pd.read_csv(path + file_prefix[i] + '.sed', sep=r'\s+', skiprows=1, header=None, names=columns, comment='#', dtype=str)
            if sedtable.bandname.isin(['Gaia_G_EDR3']).any():
                gaia_g = sedtable.magnitude[sedtable.bandname == 'Gaia_G_EDR3'].iloc[0]
                gaia_g_err = sedtable.used_errors[sedtable.bandname == 'Gaia_G_EDR3'].iloc[0]
                gaia_bp = sedtable.magnitude[sedtable.bandname == 'Gaia_BP_EDR3'].iloc[0]
                gaia_bp_err = sedtable.used_errors[sedtable.bandname == 'Gaia_BP_EDR3'].iloc[0]
                gaia_rp = sedtable.magnitude[sedtable.bandname == 'Gaia_RP_EDR3'].iloc[0]
                gaia_rp_err = sedtable.used_errors[sedtable.bandname == 'Gaia_RP_EDR3'].iloc[0]
            else:
                gaia_g = sedtable.magnitude[sedtable.bandname == 'Gaia'].iloc[0]
                gaia_g_err = sedtable.used_errors[sedtable.bandname == 'Gaia'].iloc[0]
                gaia_bp = sedtable.magnitude[sedtable.bandname == 'GaiaBP'].iloc[0]
                gaia_bp_err = sedtable.used_errors[sedtable.bandname == 'GaiaBP'].iloc[0]
                gaia_rp = sedtable.magnitude[sedtable.bandname == 'GaiaRP'].iloc[0]
                gaia_rp_err = sedtable.used_errors[sedtable.bandname == 'GaiaRP'].iloc[0]

            j_2mass = sedtable.magnitude[sedtable.bandname == 'J2M'].iloc[0]
            j_2mass_err = sedtable.used_errors[sedtable.bandname == 'J2M'].iloc[0]
            h_2mass = sedtable.magnitude[sedtable.bandname == 'H2M'].iloc[0]
            h_2mass_err = sedtable.used_errors[sedtable.bandname == 'H2M'].iloc[0]
            k_2mass = sedtable.magnitude[sedtable.bandname == 'K2M'].iloc[0]
            k_2mass_err = sedtable.used_errors[sedtable.bandname == 'K2M'].iloc[0]

            wise1 = sedtable.magnitude[sedtable.bandname == 'WISE1'].iloc[0]
            wise1_err = sedtable.used_errors[sedtable.bandname == 'WISE1'].iloc[0]
            wise2 = sedtable.magnitude[sedtable.bandname == 'WISE2'].iloc[0]
            wise2_err = sedtable.used_errors[sedtable.bandname == 'WISE2'].iloc[0]
            wise3 = sedtable.magnitude[sedtable.bandname == 'WISE3'].iloc[0]
            wise3_err = sedtable.used_errors[sedtable.bandname == 'WISE3'].iloc[0]
            if sedtable.bandname.isin(['WISE4']).any():
                wise4 = sedtable.magnitude[sedtable.bandname == 'WISE4'].iloc[0]
                wise4_err = sedtable.used_errors[sedtable.bandname == 'WISE4'].iloc[0]
            else:
                wise4 = None
                wise4_err = None
        else:
            gaia_g = data_gaia['Gmag']
            gaia_g_err = data_gaia['e_Gmag']
            gaia_bp = data_gaia['BPmag']
            gaia_bp_err = data_gaia['e_BPmag']
            gaia_rp = data_gaia['RPmag']
            gaia_rp_err = data_gaia['e_RPmag']

            j_2mass = data_2MASS['Jmag']
            j_2mass_err = data_2MASS['e_Jmag']
            h_2mass = data_2MASS['Hmag']
            h_2mass_err = data_2MASS['e_Hmag']
            k_2mass = data_2MASS['Kmag']
            k_2mass_err = data_2MASS['e_Kmag']

            wise1 = data_WISE['W1mag']
            wise1_err = data_WISE['e_W1mag']
            wise2 = data_WISE['W2mag']
            wise2_err = data_WISE['e_W2mag']
            wise3 = data_WISE['W3mag']
            wise3_err = data_WISE['e_W3mag']
            wise4 = data_WISE['W4mag']
            wise4_err = data_WISE['e_W4mag']

        # grabbing astrometric parameters
        pmra = data_gaia['pmRA']
        pmra_err = data_gaia['e_pmRA']
        pmdec = data_gaia['pmDE']
        pmdec_err = data_gaia['e_pmDE']
        parallax = data_gaia['Plx']
        parallax_err = data_gaia['e_Plx']
        vbroad = data_gaia['Vbroad']
        vbroad_err = data_gaia['e_Vbroad']

        gen_lit_str(ra_arr, ra_str)
        gen_lit_str(dec_arr, dec_str)
        gen_lit_str(gaia_g_arr, round(float(gaia_g), 3), round(float(gaia_g_err), 3))
        gen_lit_str(gaia_bp_arr, round(float(gaia_bp), 3), round(float(gaia_bp_err), 3))
        gen_lit_str(gaia_rp_arr, round(float(gaia_rp), 3), round(float(gaia_rp_err), 3))
        gen_lit_str(tmag_arr, round(float(TESS_mags[i]), 4), round(float(TESS_mags_err[i]), 4))
        gen_lit_str(j_2mass_arr, round(float(j_2mass), 3), round(float(j_2mass_err), 3))
        gen_lit_str(h_2mass_arr, round(float(h_2mass), 3), round(float(h_2mass_err), 3))
        gen_lit_str(k_2mass_arr, round(float(k_2mass), 3), round(float(k_2mass_err), 3))
        gen_lit_str(wise1_arr, round(float(wise1), 3), round(float(wise1_err), 3))
        gen_lit_str(wise2_arr, round(float(wise2), 3), round(float(wise2_err), 3))
        gen_lit_str(wise3_arr, round(float(wise3), 3), round(float(wise3_err), 3))
        if (wise4 is None) or isinstance(wise4, np.ma.core.MaskedConstant):
            gen_lit_str(wise4_arr, None) # a filler cell, so the row stays as wide as the others
        else:
            gen_lit_str(wise4_arr, round(float(wise4), 3), round(float(wise4_err), 3))
            wise4_targets += 1
        gen_lit_str(pmra_arr, round(float(pmra), 3), round(float(pmra_err), 3))
        gen_lit_str(pmdec_arr, round(float(pmdec), 3), round(float(pmdec_err), 3))
        gen_lit_str(parallax_arr, round(float(parallax), 4), round(float(parallax_err), 4))
        if vsini_type == 'gaia':
            gen_lit_str(vsini_arr, vbroad, vbroad_err)
        elif vsini_type == 'tres':
            gen_lit_str(vsini_arr, vsini_tres[i], vsini_tres_err[i])
        elif vsini_type == 'external':
            gen_lit_str(vsini_arr, vsini_external[i][0], vsini_external[i][1])

    if add_source_column == True:
        # adding sources to the rows
        add_source(ra_arr, 1) # 1 corresponds to Gaia
        add_source(dec_arr, 1)
        add_source(gaia_g_arr, 1)
        add_source(gaia_bp_arr, 1)
        add_source(gaia_rp_arr, 1)
        add_source(pmra_arr, 1)
        add_source(pmdec_arr, 1)
        add_source(parallax_arr, 1)
        add_source(vsini_arr, 2) # 2 corresponds to TRES
        add_source(tmag_arr, 3) # 3 corresponds to the TIC
        add_source(j_2mass_arr, 4) # 4 corresponds to 2MASS
        add_source(h_2mass_arr, 4)
        add_source(k_2mass_arr, 4)
        add_source(wise1_arr, 5) # 5 corresponds to WISE
        add_source(wise2_arr, 5)
        add_source(wise3_arr, 5)
        if wise4_targets > 0:
            add_source(wise4_arr, 5)

    # Generating the preamble

    preamble = (r'\providecommand{\bjdtdb}{\ensuremath{\rm {BJD_{TDB}}}}'+'\n'+
        r'\providecommand{\feh}{\ensuremath{\left[{\rm Fe}/{\rm H}\right]}}'+'\n'+
        r'\providecommand{\teff}{\ensuremath{T_{\rm eff}}}'+'\n'+
        r'\providecommand{\teq}{\ensuremath{T_{\rm eq}}}'+'\n'+
        r'\providecommand{\ecosw}{\ensuremath{e\cos{\omega_*}}}'+'\n'+
        r'\providecommand{\esinw}{\ensuremath{e\sin{\omega_*}}}'+'\n'+
        r'\providecommand\msun{M$_\odot$\xspace}'+'\n'+
        r'\providecommand{\rsun}{R$_\odot$\xspace}'+'\n'+
        r'\providecommand{\lsun}{L$_\odot$\xspace}'+'\n'+
        r'\providecommand{\mj}{\ensuremath{\,M_{\rm J}}}'+'\n'+
        r'\providecommand{\rj}{\ensuremath{\,R_{\rm J}}}'+'\n'+
        r'\providecommand{\me}{\ensuremath{\,M_{\rm E}}}'+'\n'+
        r'\providecommand{\re}{\ensuremath{\,R_{\rm E}}}'+'\n'+
        r'\providecommand{\fave}{\langle F \rangle}'+'\n'+
        r'\providecommand{\fluxcgs}{10$^9$ erg s$^{-1}$ cm$^{-2}$}'+'\n'+
        r'\providecommand{\tess}{\textit{TESS}\xspace}'+'\n')

    notes = (r'\vspace{2mm}' + '\n' +
             r'\begin{minipage}{\textwidth}' + '\n' +
             r'\textbf{Notes:}' + '\n' +
             r'\footnotesize' + '\n' +
             r'The uncertainties of the photometric measurements have a systematic floor applied that is usually larger than the reported catalog errors.\\' + '\n' +
             r'$\ddagger$ Right Ascension and Declination are in epoch J2000. Coordinates are from Vizier where Gaia RA and Dec have been precessed and corrected from epoch J2016.\\' + '\n' +
             r'Sources: (1) \cite{GaiaDR3}; (2) \S\ref{subsubsec:tres} \& \S\ref{subsubsec:chiron}; (3) \cite{Stassun:2019}; (4) \cite{Cutri:2003, Skrutskie:2006}; (5) \cite{Wright:2010, Cutri:2012}' + '\n' +
             r'\end{minipage}' + '\n')

    # Number of targets that actually made it into the rows
    n_targets = len(TIC_IDs)

    # Deciding how to split the targets across tables
    if not max_targets_per_table or max_targets_per_table >= n_targets:
        chunk_size = n_targets
    else:
        chunk_size = int(max_targets_per_table)
    n_chunks = max(1, int(np.ceil(n_targets / chunk_size)))

    # Deriving the output filename(s). The first table keeps the classic 'lit_table.tex'
    # name (bumping a numeric suffix if it already exists); continuation tables append '_2', '_3', ...
    first_name = 'lit_table.tex'
    suffix = 2
    while os.path.exists(f'{outputpath}/{first_name}'):
        first_name = f'lit_table_{suffix}.tex'
        suffix += 1
    stem = first_name[:-len('.tex')]
    filenames = [first_name] + [f'{stem}_{k}.tex' for k in range(2, n_chunks + 1)]
    print('Saving this table as ' + ', '.join(filenames) + '...')

    extra_cols = 3 if add_source_column else 2

    def _row_slice(arr, sl):
        '''Return a table row for the targets in slice `sl`: label cell(s) + sliced target
        cells + any trailing source cell.'''
        label = arr[:1]
        entries = arr[1:1 + n_targets]
        tail = arr[1 + n_targets:]  # the source cell, if add_source_column added one
        return label + entries[sl] + tail

    def _write_chunk(fname, sl, is_first, is_last):
        chunk_targets = list(target_list[sl])
        n_chunk = len(chunk_targets)
        colstring = 'c' * n_chunk + ('c' if add_source_column else '')
        namestring = ''.join(' & ' + t for t in chunk_targets)
        tic_id_str = ''.join(tic_id_list[sl])
        tycho_id_str = ''.join(tycho_id_list[sl])
        twomass_id_str = ''.join(twomass_id_list[sl])
        gaia_id_str = ''.join(gaia_id_list[sl])

        if is_first:
            caption = r'\caption{Measured Properties from Literature}'
        elif MNRAS:
            caption = r'\contcaption{}'
        else:
            caption = r'\caption{\textit{(Continued)}}'

        with open(f'{outputpath}/{fname}', 'w') as fout:
            fout.write(preamble)
            fout.write(r'\begin{table*}' + '\n' +
                       r'\centering' + '\n' +
                       caption + '\n')
            if is_first:
                fout.write(r'\label{tab:lit}' + '\n')
            fout.write(r'\resizebox{\textwidth}{!}{' + '\n')

            if add_source_column == True:
                fout.write(r'\begin{tabular}{ll' + colstring + '}'+'\n'+
                    r'\hline' + '\n' +
                    r'& ' + namestring + r' & Source \\' +'\n'+
                    r'\multicolumn{' + str(n_chunk + 3) + r'}{l}{\textbf{Other identifiers}:} \\' + '\n' +
                    r'& \tess Input Catalog' + tic_id_str + r' & \\' + '\n' +
                    r'& TYCHO-2' + tycho_id_str + r' & \\'  + '\n' +
                    r'& 2MASS' + twomass_id_str + r' & \\' + '\n' +
                    r'& Gaia DR3' + gaia_id_str + r' & \\' + '\n' +
                    r'\hline' + '\n' +
                    r'\multicolumn{' + str(n_chunk + 3) + r'}{l}{\textbf{Astrometric Parameters}:} \\' + '\n')
            else:
                fout.write(r'\begin{tabular}{ll' + colstring + '}'+'\n'+
                    r'\hline' + '\n' +
                    r'& ' + namestring + r'\\' +'\n'+
                    r'\multicolumn{' + str(n_chunk + 2) + r'}{l}{\textbf{Other identifiers}:} \\' + '\n' +
                    r'& \tess Input Catalog' + tic_id_str + r'\\' + '\n' +
                    r'& TYCHO-2' + tycho_id_str + r'\\'  + '\n' +
                    r'& 2MASS' + twomass_id_str + r'\\' + '\n' +
                    r'\hline' + '\n' +
                    r'\multicolumn{' + str(n_chunk + 2) + r'}{l}{\textbf{Astrometric Parameters}:} \\' + '\n')

            write(_row_slice(ra_arr, sl), fout)
            write(_row_slice(dec_arr, sl), fout)
            write(_row_slice(pmra_arr, sl), fout)
            write(_row_slice(pmdec_arr, sl), fout)
            write(_row_slice(parallax_arr, sl), fout)
            write(_row_slice(vsini_arr, sl), fout)
            fout.write(r'\multicolumn{' + str(n_chunk + extra_cols) + r'}{l}{\textbf{Photometric Parameters}:} \\' + '\n')
            write(_row_slice(gaia_g_arr, sl), fout)
            write(_row_slice(gaia_bp_arr, sl), fout)
            write(_row_slice(gaia_rp_arr, sl), fout)
            write(_row_slice(tmag_arr, sl), fout)
            write(_row_slice(j_2mass_arr, sl), fout)
            write(_row_slice(h_2mass_arr, sl), fout)
            write(_row_slice(k_2mass_arr, sl), fout)
            write(_row_slice(wise1_arr, sl), fout)
            write(_row_slice(wise2_arr, sl), fout)
            write(_row_slice(wise3_arr, sl), fout)
            if wise4_targets > 0:
                write(_row_slice(wise4_arr, sl), fout)

            fout.write(r'\hline' + '\n' +
                       r'\end{tabular}' + '\n' +
                       r'} % end resizebox' + '\n')
            if is_last:
                # the minipage notes block is only used in the final table
                fout.write(notes)
            fout.write(r'\end{table*}')
            if not is_last and not MNRAS:
                # keep every piece except the last on the same table number, which \contcaption does itself
                fout.write('\n' + r'\addtocounter{table}{-1}')

    for k, fname in enumerate(filenames):
        sl = slice(k * chunk_size, min((k + 1) * chunk_size, n_targets))
        _write_chunk(fname, sl, is_first=(k == 0), is_last=(k == n_chunks - 1))
