import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
from astroquery.vizier import Vizier
from astropy.coordinates import Angle
from grab_tres_vsini import grab_tres_vsini
import re
import ast
from urllib.request import urlopen
import warnings

def remove_sci_notation(x):
    '''
    Removes scientific notation from a number and returns it as a string.
    
    Parameters
    -----------
    x: the number to remove scientific notation from
    '''
    return np.format_float_positional(x, trim='-')

def round_sig_figs(x, num_sig_figs):
    '''
    Rounds a number to a specified number of significant figures.

    Parameters
    -----------
    x: the number to round
    num_sig_figs: the number of significant figures to round to
    '''
    return '{:g}'.format(float('{:.{p}g}'.format(x, p=num_sig_figs)))

def robust_decimal_errors(val, up_err, low_err):
    '''
    Ensures that the uncertainties have at least as many positions behind the decimal as the value itself.

    Parameters
    -----------
    val: the value to round
    up_err: the upper error on the value
    low_err: the lower error on the value
    '''
    val_str = str(remove_sci_notation(val))
    up_err_str = str(remove_sci_notation(up_err))
    low_err_str = str(remove_sci_notation(low_err))

    if '.' in val_str:
        val_decimal_places = len(val_str.split('.')[1])
    else:
        val_decimal_places = 0

    if '.' in low_err_str:
        low_err_decimal_places = len(low_err_str.split('.')[1])
    else:
        low_err_decimal_places = 0

    if '.' in up_err_str:
        up_err_decimal_places = len(up_err_str.split('.')[1])
    else:
        up_err_decimal_places = 0
    
    if val_decimal_places > 0:
        if low_err_decimal_places < val_decimal_places:
            if '.' not in low_err_str:
                low_err_str += '.'
            low_err_str = low_err_str + '0' * (val_decimal_places - low_err_decimal_places)
        if up_err_decimal_places < val_decimal_places:
            if '.' not in up_err_str:
                up_err_str += '.'
            up_err_str = up_err_str + '0' * (val_decimal_places - up_err_decimal_places)
    
    if low_err_decimal_places > 0:
        if low_err_decimal_places < up_err_decimal_places: # pad zeros to the lower error if it has fewer decimal places than the upper error
            if '.' not in low_err_str:
                low_err_str += '.'
            low_err_str = low_err_str + '0' * (up_err_decimal_places - low_err_decimal_places)
    if up_err_decimal_places > 0:
        if up_err_decimal_places < low_err_decimal_places:
            if '.' not in up_err_str:
                up_err_str += '.'
            up_err_str = up_err_str + '0' * (low_err_decimal_places - up_err_decimal_places)

    # Final check to remove trailing zeros and ensure same number of decimal places in errors
    if '.' in low_err_str:
        low_err_decimal_places = len(low_err_str.split('.')[1])
    if '.' in up_err_str:
        up_err_decimal_places = len(up_err_str.split('.')[1])
    if low_err_decimal_places != up_err_decimal_places:
        while (low_err_decimal_places > up_err_decimal_places) and (low_err_str[-1] == '0'):
            low_err_str = low_err_str[:-1]
            low_err_decimal_places -= 1
        while (up_err_decimal_places > low_err_decimal_places) and (up_err_str[-1] == '0'):
            up_err_str = up_err_str[:-1]
            up_err_decimal_places -= 1

    if val_decimal_places == 1 and low_err_decimal_places == 0 and up_err_decimal_places == 0:
        val = int(val) # to remove trailing zeros
        val_str = str(val)

    return val_str, up_err_str, low_err_str

def format_value_with_errors(val, up_err, low_err):
    '''
    Format a value and its asymmetric uncertainties so:
    - scientific notation is removed
    - the uncertainties have the same number of decimal places as the value

    Parameters
    -----------
    val: the value to format
    up_err: the upper uncertainty on the value
    low_err: the lower uncertainty on the value
    '''
    try:
        val_f = float(val)
    except Exception:
        return '---', None, None

    if not np.isfinite(val_f):
        return '---', None, None

    # try to coerce errors to floats; if they are not finite, treat as missing
    try:
        up_f = float(up_err)
    except Exception:
        up_f = np.nan
    try:
        low_f = float(low_err)
    except Exception:
        low_f = np.nan

    if (not np.isfinite(up_f)) or (not np.isfinite(low_f)):
        # No reliable errors: just return the value without sci notation
        val_str = remove_sci_notation(val_f)
        return val_str, None, None

    # Use robust_decimal_errors to line up decimal places and remove sci notation
    val_str, up_str, low_str = robust_decimal_errors(val_f, up_f, low_f)

    return val_str, up_str, low_str

def grab_medians(path, file_prefix, bimodal=False):
    '''
    Collects median values from EXOFASTv2 output files at the defined path.

    Parameters
    -----------
    path: path to the fit files. Ex: '/Users/jack/Research/pipelines/system_figure_pipeline/data/'
    file_prefix: prefix used in EXOFASTv2 output files. Ex: 'toi1855'
    bimodal: boolean to decide whether or not to use the files that are output by splitpdf.pro
    '''

    median_names= ['parname', 'median_value', 'upper_error', 'lower_error', 'scinot']
    if bimodal == False:
        medians = pd.read_csv(path + file_prefix + '.median.csv', names=median_names, header=None, skiprows=1)
    else:
        medians = pd.read_csv(path + file_prefix + '.csv', names=median_names, header=None, skiprows=1)

    medians_corrected = medians.copy()
    for i in range(len(medians_corrected)):
        median_corrections = median_scinot_corrections(medians, medians_corrected.parname[i])
        median_corrections = [float(x) for x in median_corrections]
        medians_corrected.loc[i, 'median_value'] = median_corrections[0]
        medians_corrected.loc[i, 'upper_error'] = median_corrections[1]
        medians_corrected.loc[i, 'lower_error'] = median_corrections[2]
    return medians_corrected

def median_scinot_corrections(median, parname):
    '''
    Multiplies parameters in the EXOFASTv2 median table by the scientific notation exponent.

    median: pandas DataFrame for median table
    param: input parameter name
    '''

    scinot = median.scinot[median.parname==parname].iloc[0]

    if type(scinot) == str:
        exp_search = re.findall(r'\\times 10\^{(.*)}', scinot)
        exponent = int(exp_search[0])
    else:
        exponent = 0

    param = median.median_value[median.parname==parname].iloc[0]
    param_corrected = param * 10**exponent

    uperr = median.upper_error[median.parname==parname].iloc[0]
    uperr_corrected = uperr * 10**exponent

    lowerr = median.lower_error[median.parname==parname].iloc[0]
    lowerr_corrected = lowerr * 10**exponent
    return robust_decimal_errors(param_corrected, uperr_corrected, lowerr_corrected)

def grab_priors(file_prefix, path):
    '''
    Collects prior values from a prior file labeled 'toi####.priors.final'

    Parameters
    -----------
    file_prefix: prefix used in EXOFASTv2 output files. Ex: 'toi1855'
    path: path to the fit files. Ex: '/Users/jack/Research/pipelines/system_figure_pipeline/data/'
    '''

    columns = ['variable', 'meanvalue', 'stdev', 'low_bound', 'up_bound', 'starting_value'] # these column names only make sense for gaussian column names
    priors = pd.read_csv(path + file_prefix + '.priors.final', sep=r'\s+', skiprows=1, header=None, comment='#', names=columns, on_bad_lines='skip')

    for i in range(len(priors)):
        # find linked parameters (meanvalue holds the name of another variable)
        # and replace them with that variable's value
        val = priors.meanvalue[i]
        if not isinstance(val, str):
            continue
        # try the exact name first, then fall back to stripping the '_0' suffix
        for candidate in (val, val.replace('_0', '')):
            match = priors.loc[priors.variable == candidate, 'meanvalue']
            if len(match):
                priors.loc[i, 'meanvalue'] = match.iloc[0]
                break
    priors['meanvalue'] = priors['meanvalue'].astype(float) # to ensure that all mean values are floats
    return priors

def make_median_string(medians, param, array):
    '''
    Turns median values into strings containing LaTeX, ready to be entered into the median table.

    Parameters
    -----------
    medians: Pandas DataFrame containing the median values obtained using the grab_medians function.
    param: the parameter to generate a string for
    array: the array corresponding to the table row that the parameter should be appended to

    Exactly one cell string is appended per call so that the row can be sliced by target index
    (e.g. when splitting a long target list across multiple tables).
    '''
    param = param+'_0'

    if medians.parname.isin([param]).any() == True:
        val = medians.median_value[medians.parname == param].iloc[0]
        uperr = medians.upper_error[medians.parname == param].iloc[0]
        loerr = medians.lower_error[medians.parname == param].iloc[0]

        val_str, up_str, low_str = format_value_with_errors(val, uperr, loerr)

        if val_str == '---':
            array.append('& ---')
            return

        if (up_str is None) or (low_str is None):
            # no reliable errors
            array.append('& $' + val_str + '$ ')
            return

        if float(uperr) == float(loerr):
            errstring = r' \pm ' + up_str
        else:
            errstring = r'^{+' + up_str + '}_{-' + low_str + '}'

        array.append('& $' + val_str + errstring + '$ ')
    else:
        array.append('& ---')

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

def write(param_arr,file):
    for ii in param_arr:
        file.write(ii)
    file.write(r'\\'+'\n')

def lit_table(target_list, path, file_prefix=None, outputpath='.', vsini_type='gaia', vsini_external=None, tres_username=None, tres_password=None,
              add_source_column=False, grab_mags_from_sedfile=True, max_targets_per_table=5):
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
            warnings.warn(f'Multiple sources in Vizier Gaia query of TIC {TIC_IDs[i]}. Selecting closest source.')
        data_gaia = data_gaia[0]

        twomass_columns=['_r', 'Jmag', 'e_Jmag', 'Hmag', 'e_Hmag', 'Kmag', 'e_Kmag']
        v2mass = Vizier(columns=twomass_columns, catalog='II/246/out')
        data_2MASS = v2mass.query_region('TIC ' + str(TIC_IDs[i]), radius=Angle(6, "arcsec"))
        data_2MASS = data_2MASS[0] # Vizier returns a list of tables, but we only want the first one
        data_2MASS.sort('_r') # sort by distance from the target star
        if len(data_2MASS) > 1:
            warnings.warn(f'Multiple sources in Vizier 2MASS query of TIC {TIC_IDs[i]}. Selecting closest source.')
        data_2MASS = data_2MASS[0]
        

        WISE_columns=['_r', 'W1mag', 'e_W1mag', 'W2mag', 'e_W2mag', 'W3mag', 'e_W3mag', 'W4mag', 'e_W4mag']
        vwise = Vizier(columns=WISE_columns, catalog='II/311/wise')
        data_WISE = vwise.query_region('TIC ' + str(TIC_IDs[i]), radius=Angle(6, "arcsec"))
        data_WISE = data_WISE[0] # Vizier returns a list of tables, but we only want the first one
        data_WISE.sort('_r') # sort by distance from the target star
        if len(data_WISE) > 1:
            warnings.warn(f'Multiple sources in Vizier WISE query of TIC {TIC_IDs[i]}. Selecting closest source.')
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
            wise4count = 0 # WISE4 magnitudes are often not reported or used for any targets. This variable keeps track of the WISE4 mags in fits

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
                wise4count += 1
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
        if wise4count > 0:
            gen_lit_str(wise4_arr, round(float(wise4), 3), round(float(wise4_err), 3))
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
        if wise4count > 0:
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
        colstring = 'cc' + 'c' * n_chunk
        namestring = ''.join(' & ' + t for t in chunk_targets)
        tic_id_str = ''.join(tic_id_list[sl])
        tycho_id_str = ''.join(tycho_id_list[sl])
        twomass_id_str = ''.join(twomass_id_list[sl])
        gaia_id_str = ''.join(gaia_id_list[sl])

        caption = (r'\caption{Measured Properties from Literature}' if is_first
                   else r'\caption{\textit{(Continued)}}')

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
                    r'& \tess Input Catalog' + tic_id_str + r'\\' + '\n' +
                    r'& TYCHO-2' + tycho_id_str + r' & \\'  + '\n' +
                    r'& 2MASS' + twomass_id_str + r' & \\' + '\n' +
                    r'& Gaia DR3' + gaia_id_str + r' & \\' + '\n' +
                    r'\hline' + '\n' +
                    r'\multicolumn{' + str(n_chunk + 3) + r'}{l}{\textbf{Astrometric Parameters}:} \\' + '\n')
            else:
                fout.write(r'\begin{tabular}{l l' + colstring + '}'+'\n'+
                    r'\hline' + '\n' +
                    r'& ' + namestring + r'\\' +'\n'+
                    r'\multicolumn{' + str(n_chunk + 2) + r'}{l}{\textbf{Other identifiers}:} \\' + '\n' +
                    r'& \tess Input Catalog' + tic_id_str + r' & \\' + '\n' +
                    r'& TYCHO-2' + tycho_id_str + r' & \\'  + '\n' +
                    r'& 2MASS' + twomass_id_str + r' & \\' + '\n' +
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
            if add_source_column == True:
                if wise4count > 0:
                    write(_row_slice(wise4_arr, sl), fout)
            else:
                write(_row_slice(wise4_arr, sl), fout)

            fout.write(r'\hline' + '\n' +
                       r'\end{tabular}' + '\n' +
                       r'} % end resizebox' + '\n')
            if is_last:
                # the minipage notes block is only used in the final table
                fout.write(notes)
            fout.write(r'\end{table*}')
            if not is_last:
                # keep every piece except the last on the same table number
                fout.write('\n' + r'\addtocounter{table}{-1}')

    for k, fname in enumerate(filenames):
        sl = slice(k * chunk_size, min((k + 1) * chunk_size, n_targets))
        _write_chunk(fname, sl, is_first=(k == 0), is_last=(k == n_chunks - 1))



def med_table(target_list, path, file_prefix_list, outputpath='.', bimodal=False, parameters=None,
              probabilities=None, max_targets_per_table=5):
    '''
    Generates a median table given the path to EXOFASTv2 output files.

    Parameters
    -----------
    target_list: an array of strings containing the names of each target. Names should match those in EXOFASTv2 output files.
    path: path of EXOFASTv2 output files
    file_prefix_list: list of prefixes used in EXOFASTv2 file generation
    outputpath: the folder in which the table should be generated. Current working directory by default
    bimodal: boolean to decide whether or not to use the files that are output by splitpdf.pro. When True the output
        files are named bimodal_median_table.tex (rather than median_table.tex), the low- and high-mass solutions of a
        given target are kept together in the same table when the target list is split (see max_targets_per_table), and
        the header is written three rows deep: the system name spanning its solution columns, then "Low-mass solution"
        or "High-mass solution" for each column, then the probabilities if they were given. The system name and the
        solution label are read off the target name ('TOI-4138 (Low Mass)'), falling back to the file prefix
        ('257060897.lowmass') when the target name carries no low-/high-mass marker.
    parameters: optional list of parameter names to include in the table. Defaults to the set of parameters
        used in Schulte+ 2025
    probabilities: optional list of the probability of each solution, ordered to match file_prefix_list. Each entry is
        written as "X\\% probability" on its own header line, directly beneath the solution the column holds. Entries
        may be numbers understood as percentages (32 or 32.5), pre-formatted strings ('32%', '~50\\% probability'), or
        None for a blank cell.
    max_targets_per_table: the maximum number of targets shown in a single table. If len(target_list) exceeds this value, the
        targets are split across multiple median_table .tex files. Every table after the first is captioned "\\textit{(Continued)}"
        and every table except the last gets "\\addtocounter{table}{-1}" so that all pieces share one table number. The
        \\begin{flushleft} notes block is only written in the last table. Set to None (or 0) to force a single table.
    '''

    def _normalize_param_name(name):
        return re.sub(r'[^a-z0-9]', '', str(name).lower())

    def _format_probability(prob):
        '''Turns an entry of `probabilities` into a header cell such as '32\\% probability'.'''
        if prob is None:
            return ''
        if isinstance(prob, str):
            text = prob.strip()
            if not text:
                return ''
        else:
            if isinstance(prob, float) and np.isnan(prob):
                return ''
            text = f'{float(prob):g}'
        text = re.sub(r'(?<!\\)%', r'\\%', text) # escape bare percent signs so LaTeX does not comment out the line
        if 'probability' not in text.lower():
            if not text.endswith(r'\%'):
                text += r'\%'
            text += ' probability'
        return text

    if probabilities is not None and len(probabilities) != len(file_prefix_list):
        raise ValueError(f'probabilities has {len(probabilities)} entries but file_prefix_list has {len(file_prefix_list)}; '
                         'the two must line up index for index.')
    probability_cells = [_format_probability(p) for p in probabilities] if probabilities is not None else None

    default_parameters = [
        'mstar', 'rstar', 'lstar', 'rhostar', 'logg', 'teff', 'feh', 'initfeh', 'age', 'eep', 'Av',
        'distance', 'Period', 'rp', 'mp', 'tc', 't0', 'a', 'ideg', 'e', 'omegadeg', 'teq', 'tcirc',
        'k', 'slope', 'p', 'ar', 'depth_TESS', 'tau', 't14', 'b', 'rhop', 'loggp', 'q', 'dr'
    ]

    if parameters is None:
        parameters = default_parameters
    selected_parameters = {_normalize_param_name(param) for param in parameters}

    # Initializing strings with LaTeX for each parameter
    mstars=[r'$M_*$ & Mass (\msun) ']
    rstars=[r'$R_*$ & Radius (\rsun) ']
    rstarseds=[r'$R_{*,SED}$ & Radius$ (\rsun) ']
    lstars=[r'$L_*$ & Luminosity (\lsun) ']
    fbols=[r'$F_{\rm Bol}$ & Bolometric flux (cgs) ']
    rhostars=[r'$\rho_*$ & Density (cgs) ']
    loggs=[r'$\log{g}$ & Surface gravity (cgs) ']
    teffs=[r'$T_{\rm eff}$ & Effective temperature (K) ']
    teffseds=[r'$T_{\rm eff,SED}$ & Effective temperature (K) ']
    fehs=[r'$[{\rm Fe/H}]$ & Metallicity (dex) ']
    initfehs=[r'$[{\rm Fe/H}]_{0}$ & Initial metallicity (dex) ']
    ages=[r'Age & Age (Gyr) ']
    eeps=[r'EEP & Equivalent evolutionary phase ']
    logmstars=[r'$\log{M_*}$ & Mass ($\log{\msun}$) ']
    avs=[r'$A_V$ & V-band extinction (mag) ']
    errscales=[r'$\sigma_{SED}$ & SED photometry error scaling ']
    plaxes=[r'$\varpi$ & Parallax (mas) ']
    dists=[r'$d$ & Distance (pc) ']

    periods=[r'$P$ & Period (days) ']
    rps=[r'$R_{\rm P}$ & Radius (\rj) ']
    mps=[r'$M_{\rm P}$ & Mass (\mj) ']
    mpsuns=[r'$M_{\rm P}$ & Mass (\msun) ']
    tcs=[r'$T_C$ & Time of conjunction (\bjdtdb) ']
    tts=[r'$T_T$ & Time of minimum projected separation (\bjdtdb) ']
    t0s=[r'$T_0$ & Optimal conjunction time (\bjdtdb) ']
    semimajors=[r'$a$ & Semi-major axis (AU) ']
    idegs=[r'$i$ & Inclination (Degrees) ']
    eccs=[r'$e$ & Eccentricity ']
    odegs=[r'$\omega_*$ & Argument of periastron (Degrees) ']
    teqs=[r'$\teq$ & Equilibrium temperature (K) ']
    tcircs=[r'$\tau_{\rm circ}$ & Tidal circularization timescale (Gyr) ']
    ks=[r'$K$ & RV semi-amplitude (m/s) ']
    slopes=[r'$\dot{\gamma}$ & RV slope (m/s/day) ']
    ps=[r'$R_{\rm P}/R_*$ & Radius of planet in stellar radii  ']
    ars=[r'$a/R_*$ & Semi-major axis in stellar radii  ']
    deltas=[r'$\delta$ & $\left(R_{\rm P}/R_*\right)^2$ ']
    Tdepths=[r'Depth & \tess flux decrement at mid-transit ']
    taus=[r'$\tau$ & Ingress/egress transit duration (days) ']
    t14s=[r'$T_{14}$ & Total transit duration (days) ']
    tfwhms=[r'$T_{FWHM}$ & FWHM transit duration (days) ']
    bs=[r'$b$ & Transit impact parameter ']
    cosis=[r'$\cos i$ & Cosine of inclination ']
    bss=[r'$b_S$ & Eclipse impact parameter ']
    tauss=[r'$\tau_S$ & Ingress/egress eclipse duration (days) ']
    ts14s=[r'$T_{S,14}$ & Total eclipse duration (days) ']
    tfwhmss=[r'$T_{S,FWHM}$ & FWHM eclipse duration (days) ']
    depth25s=[r'$\delta_{S,2.5\mu m}$ & Blackbody eclipse depth at 2.5$\mu$m (ppm) ']
    depth50s=[r'$\delta_{S,5.0\mu m}$ & Blackbody eclipse depth at 5.0$\mu$m (ppm) ']
    depth75s=[r'$\delta_{S,7.5\mu m}$ & Blackbody eclipse depth at 7.5$\mu$m (ppm) ']
    rhops=[r'$\rho_{\rm P}$ & Density (cgs) ']
    loggps=[r'$\log{g_{\rm P}}$ & Surface gravity (cgs) ']
    safronovs=[r'$\Theta$ & Safronov number ']
    faves=[r'$\fave$ & Incident flux (\fluxcgs) ']
    tps=[r'$T_{\rm P}$ & Time of periastron (\bjdtdb) ']
    tss=[r'$T_{\rm S}$ & Time of eclipse (\bjdtdb) ']
    tas=[r'$T_{\rm A}$ & Time of ascending node (\bjdtdb) ']
    tds=[r'$T_{\rm D}$ & Time of descending node (\bjdtdb) ']
    vcves=[r'$V_c/V_e$ ']
    ecosws=[r'$e\cos{\omega_*}$ & ']
    esinws=[r'$e\sin{\omega_*}$ & ']
    msinis=[r'$M_{\rm P}\sin{i}$ & Minimum mass (\mj) ']
    qs=[r'$M_{\rm P}/M_*$ & Mass ratio  ']
    drs=[r'$d/R_*$ & Separation at mid-transit  ']
    pts=[r'$P_{\rm T}$ & A priori non-grazing transit prob  ']
    ptgs=[r'$P_{\rm T,G}$ & A priori transit prob  ']
    pss=[r'$P_{\rm S}$ & A priori non-grazing eclipse prob ']
    psgs=[r'$P_{\rm S,G}$ & A priori eclipse prob  ']

    stellar_parameter_specs = [
        ('mstar', 'mstar', mstars),
        ('rstar', 'rstar', rstars),
        ('rstarsed', 'rstarsed', rstarseds),
        ('lstar', 'lstar', lstars),
        ('fbol', 'fbol', fbols),
        ('rhostar', 'rhostar', rhostars),
        ('logg', 'logg', loggs),
        ('teff', 'teff', teffs),
        ('teffsed', 'teffsed', teffseds),
        ('feh', 'feh', fehs),
        ('initfeh', 'initfeh', initfehs),
        ('age', 'age', ages),
        ('eep', 'eep', eeps),
        ('logmstar', 'logmstar', logmstars),
        ('Av', 'Av', avs),
        ('errscale', 'errscale', errscales),
        ('parallax', 'parallax', plaxes),
        ('distance', 'distance', dists),
    ]

    planetary_parameter_specs = [
        ('Period', 'Period', periods),
        ('rp', 'rp', rps),
        ('mp', 'mp', mps),
        ('mpsun', 'mpsun', mpsuns),
        ('tc', 'tc', tcs),
        ('tt', 'tt', tts),
        ('t0', 't0', t0s),
        ('a', 'a', semimajors),
        ('ideg', 'ideg', idegs),
        ('e', 'e', eccs),
        ('omegadeg', 'omegadeg', odegs),
        ('teq', 'teq', teqs),
        ('tcirc', 'tcirc', tcircs),
        ('k', 'k', ks),
        ('slope', 'slope', slopes),
        ('p', 'p', ps),
        ('ar', 'ar', ars),
        ('delta', 'delta', deltas),
        ('depth_TESS', 'depth_TESS', Tdepths),
        ('tau', 'tau', taus),
        ('t14', 't14', t14s),
        ('tfwhm', 'tfwhm', tfwhms),
        ('b', 'b', bs),
        ('cosi', 'cosi', cosis),
        ('bs', 'bs', bss),
        ('taus', 'taus', tauss),
        ('t14s', 't14s', ts14s),
        ('tfwhms', 'tfwhms', tfwhmss),
        ('eclipsedepth25', 'eclipsedepth25', depth25s),
        ('eclipsedepth50', 'eclipsedepth50', depth50s),
        ('eclipsedepth75', 'eclipsedepth75', depth75s),
        ('rhop', 'rhop', rhops),
        ('loggp', 'loggp', loggps),
        ('safronov', 'safronov', safronovs),
        ('fave', 'fave', faves),
        ('tp', 'tp', tps),
        ('ts', 'ts', tss),
        ('ta', 'ta', tas),
        ('td', 'td', tds),
        ('vcve', 'vcve', vcves),
        ('ecosw', 'ecosw', ecosws),
        ('esinw', 'esinw', esinws),
        ('msini', 'msini', msinis),
        ('q', 'q', qs),
        ('dr', 'dr', drs),
        ('pt', 'pt', pts),
        ('ptg', 'ptg', ptgs),
        ('ps', 'ps', pss),
        ('psg', 'psg', psgs),
    ]

    parameter_specs = stellar_parameter_specs + planetary_parameter_specs

    # Setting up to save the table as a .tex file

    if os.path.exists(outputpath) == False:
        os.mkdir(outputpath)

    # Deciding how to split the targets across tables
    n_targets = len(target_list)

    # Columns are first collected into groups that the split is not allowed to break up. Normally every
    # target is its own group, but in bimodal mode the low- and high-mass solutions of a target belong to
    # the same group so that they always land in the same piece of the table.
    # matches the low-/high-mass markers used in target names and in the file prefixes written by
    # splitpdf.pro, e.g. ' (Low-mass solution)', ' (High Mass)', '.lowmass', '_highmass'
    solution_marker = re.compile(r'[\s._-]*[(\[]?\s*(low|high)[\s._-]*mass\b(?:[\s._-]*(?:solutions?|solns?))?\s*[)\]]?',
                                 re.IGNORECASE)

    def _strip_solution_marker(text):
        '''Removes any low-/high-mass marker from `text`, tidying the whitespace the removal leaves behind.'''
        return re.sub(r'\s+', ' ', solution_marker.sub(' ', text)).strip()

    def _solution_group_key(index):
        '''Returns a key identifying the target that column `index` belongs to.

        The low-/high-mass markers written by splitpdf.pro are stripped from the file prefix, so that
        e.g. '257060897.lowmass' and '257060897.highmass' share a key. If the prefix carries no marker
        the target name is tried instead ('TOI-4138 (Low Mass)' -> 'toi-4138'), and a column with no
        marker at all is left in a group of its own.
        '''
        for candidate in (str(file_prefix_list[index]), str(target_list[index])):
            if solution_marker.search(candidate):
                stripped = _strip_solution_marker(candidate)
                if stripped:
                    return stripped.lower()
        return f'__ungrouped_{index}__'

    def _system_name(index):
        '''Returns the system name heading column `index`: its target name with the low-/high-mass marker
        removed ('TOI-4138 (Low Mass)' -> 'TOI-4138'). A target name carrying no marker is used as it is.
        '''
        name = str(target_list[index])
        return _strip_solution_marker(name) or name

    def _solution_label(index):
        '''Returns 'Low-mass solution' or 'High-mass solution' for column `index`, or '' if it is neither.'''
        for candidate in (str(target_list[index]), str(file_prefix_list[index])):
            match = solution_marker.search(candidate)
            if match:
                return match.group(1).capitalize() + '-mass solution'
        return ''

    groups = []
    if bimodal:
        group_index = {}
        for ii in range(n_targets):
            key = _solution_group_key(ii)
            if key in group_index:
                groups[group_index[key]].append(ii)
            else:
                group_index[key] = len(groups)
                groups.append([ii])
    else:
        groups = [[ii] for ii in range(n_targets)]

    # Packing the groups into tables. A group is never split, so a group that is on its own larger than
    # max_targets_per_table gets a table to itself and overruns the limit.
    if not max_targets_per_table or max_targets_per_table >= n_targets:
        chunks = [list(groups)] # everything fits in one table
    else:
        limit = int(max_targets_per_table)
        chunks = []
        current = []
        for group in groups:
            if current and sum(len(g) for g in current) + len(group) > limit:
                chunks.append(current)
                current = []
            current.append(group)
        if current:
            chunks.append(current)
    n_chunks = len(chunks)

    # Deriving the output filename(s). The first table keeps the classic 'median_table.tex' name
    # ('bimodal_median_table.tex' for a bimodal table), bumping a numeric suffix if it already exists;
    # continuation tables append '_2', '_3', ...
    base_name = 'bimodal_median_table' if bimodal else 'median_table'
    first_name = f'{base_name}.tex'
    suffix = 2
    while os.path.exists(f'{outputpath}/{first_name}'):
        first_name = f'{base_name}_{suffix}.tex'
        suffix += 1
    stem = first_name[:-len('.tex')]
    filenames = [first_name] + [f'{stem}_{k}.tex' for k in range(2, n_chunks + 1)]
    print('Saving this table as ' + ', '.join(filenames) + '...')

    for ii in range(len(target_list)):
        medians = grab_medians(path=path, file_prefix=file_prefix_list[ii], bimodal=bimodal)

        for param_name, median_key, labels in parameter_specs:
            if _normalize_param_name(param_name) in selected_parameters or _normalize_param_name(median_key) in selected_parameters:
                make_median_string(medians, median_key, labels)

    # Collecting priors to put at the top of the table (one ' & ...' entry per target so they can be sliced)
    parallax_prior = [] # initializing lists
    metallicity_prior = []
    extinction_prior = []
    dilution_prior = []
    show_priors = (bimodal == False) # the prior files are only read for a normal (non-split) fit
    if show_priors:
        dilute_bool = np.zeros_like(target_list) # to keep track of which targets were fit for dilution

        def prior_value(priortable, column, variable):
            # EXOFASTv2 prior files are inconsistent about capitalization (e.g. 'Av' vs 'av'),
            # so match the variable name case-insensitively; return NaN if it is absent.
            rows = priortable[priortable.variable.str.lower() == variable.lower()]
            return rows[column].iloc[0] if len(rows) else np.nan

        for ii in range(len(target_list)):
            priortable = grab_priors(file_prefix_list[ii], path)
            parallax_prior_mean = prior_value(priortable, 'meanvalue', 'parallax')
            parallax_prior_stdev = prior_value(priortable, 'stdev', 'parallax')
            parallax_prior.append(r'& $\mathcal{G}$[' + round_sig_figs(parallax_prior_mean, 5) + r', ' + round_sig_figs(parallax_prior_stdev, 5) + r'] ')
            metallicity_prior_mean = prior_value(priortable, 'meanvalue', 'feh')
            metallicity_prior_stdev = prior_value(priortable, 'stdev', 'feh')
            metallicity_prior.append(r'& $\mathcal{G}$[' + round_sig_figs(metallicity_prior_mean, 5) + r', ' + round_sig_figs(metallicity_prior_stdev, 5) + r'] ')
            extinction_prior_upperbound = prior_value(priortable, 'up_bound', 'Av')
            extinction_prior.append(r'& $\mathcal{U}$[0, ' + round_sig_figs(extinction_prior_upperbound, 5) + r'] ')

            for x in priortable.variable: # finding the dilution term
                if 'dilute' in x:
                    # use the first dilution term (dilute_0); later dilute_N rows are
                    # linked to it and carry a zeroed stdev
                    dilution_prior_mean = priortable.meanvalue[priortable.variable == x].iloc[0]
                    dilution_prior_stdev = priortable.stdev[priortable.variable == x].iloc[0]

                    dilute_bool[ii] = 1
                    break
            if dilute_bool[ii]:
                dilution_prior.append(r'& $\mathcal{G}$[' + round_sig_figs(dilution_prior_mean, 5) + r', ' + remove_sci_notation(float(round_sig_figs(dilution_prior_stdev, 5))) + r'] ')
                # above line should be cleaned up in a future version. Maybe make a new function that removes scientific notation and sets sig figs for all numbers
            else:
                dilution_prior.append(r'& --- ')

    # Generating the preamble

    preamble = (r'\providecommand{\bjdtdb}{\ensuremath{\rm {BJD_{TDB}}}}'+'\n'+
    r'\providecommand{\feh}{\ensuremath{\left[{\rm Fe}/{\rm H}\right]}}'+'\n'+
    r'\providecommand{\teff}{\ensuremath{T_{\rm eff}}}'+'\n'+
    r'\providecommand{\teq}{\ensuremath{T_{\rm eq}}}'+'\n'+
    r'\providecommand{\ecosw}{\ensuremath{e\cos{\omega_*}}}'+'\n'+
    r'\providecommand{\esinw}{\ensuremath{e\sin{\omega_*}}}'+'\n'+
    r'\providecommand{\msun}{\ensuremath{\,M_\Sun}}'+'\n'+
    r'\providecommand{\rsun}{\ensuremath{\,R_\Sun}}'+'\n'+
    r'\providecommand{\lsun}{\ensuremath{\,L_\Sun}}'+'\n'+
    r'\providecommand{\mj}{\ensuremath{\,M_{\rm J}}}'+'\n'+
    r'\providecommand{\rj}{\ensuremath{\,R_{\rm J}}}'+'\n'+
    r'\providecommand{\me}{\ensuremath{\,M_{\rm E}}}'+'\n'+
    r'\providecommand{\re}{\ensuremath{\,R_{\rm E}}}'+'\n'+
    r'\providecommand{\fave}{\langle F \rangle}'+'\n'+
    r'\providecommand{\fluxcgs}{10$^9$ erg s$^{-1}$ cm$^{-2}$}'+'\n'+
    r'\providecommand{\tess}{\textit{TESS}\xspace}'+'\n')

    def _row_slice(arr, idx):
        '''Return a table row for the targets at indices `idx`: the label cell(s) followed by those targets' cells.'''
        cells = arr[1:]
        return arr[:1] + [cells[i] for i in idx]

    def _write_chunk(fname, chunk_groups, is_first, is_last):
        idx = [i for group in chunk_groups for i in group]
        n_chunk = len(idx)
        colstring = 'lc' + 'c' * n_chunk
        if bimodal:
            # one system name per group, spanning that system's solution columns
            namestring = ''.join(r' & \multicolumn{' + str(len(group)) + r'}{c}{' + _system_name(group[0]) + r'}'
                                 for group in chunk_groups)
        else:
            namestring = ''.join(' & ' + str(target_list[i]) for i in idx)
        title = (r'Median Values and 68\% Confidence Intervals for Solutions which are Bimodal in Mass' if bimodal
                 else r'Median Values and 68\% Confidence Intervals for Fitted Stellar and Planetary Parameters')
        caption = r'\caption{' + title + '}' if is_first else r'\caption{\textit{(Continued)}}'

        with open(f'{outputpath}/{fname}', 'w') as fout:
            fout.write(preamble)
            fout.write(r'\begin{table*}' + '\n' +
                       r'\centering' + '\n' +
                       caption + '\n')
            if is_first:
                fout.write(r'\label{tab:' + ('bimodal' if bimodal else 'median') + '}' + '\n')
            fout.write(r'\scriptsize' + '\n' +
                       r'\begin{tabular}{ll' + colstring + '}'+'\n'+
                       r'\hline' + '\n' +
                       r'& ' + namestring + r'\\' +'\n')
            if bimodal:
                # which solution each column holds, under the system name it belongs to
                fout.write(r'& ' + ''.join(' & ' + _solution_label(i) for i in idx) + r'\\' + '\n')
            if probability_cells is not None:
                # the probability of each solution goes on its own header line, under the solution name
                fout.write(r'& ' + ''.join(' & ' + probability_cells[i] for i in idx) + r'\\' + '\n')
            fout.write(r'\hline' + '\n')
            if show_priors:
                # the priors are only read for a normal fit, so a bimodal table has no priors block
                fout.write(r'\multicolumn{' + str(n_chunk + 2) + r'}{l}{\textbf{Priors}:} \\' + '\n' +
                           r'$\pi$ & Gaia Parallax (mas)' + ''.join(parallax_prior[i] for i in idx) + r'\\' + '\n' +
                           r'$[{\rm Fe/H}]$ & Metallicity (dex)' + ''.join(metallicity_prior[i] for i in idx) + r'\\' + '\n' +
                           r'$A_V$ & V-band extinction (mag)' + ''.join(extinction_prior[i] for i in idx) + r'\\' + '\n' +
                           r'$D_T$ & Dilution in \tess' + ''.join(dilution_prior[i] for i in idx) + r'\\' + '\n' +
                           r'\hline' + '\n')
            fout.write(r'\multicolumn{' + str(n_chunk + 2) + r'}{l}{\textbf{Stellar Parameters}:} \\' + '\n')

            for param_name, median_key, labels in stellar_parameter_specs:
                if _normalize_param_name(param_name) in selected_parameters:
                    write(_row_slice(labels, idx), fout)

            fout.write(r'\multicolumn{' + str(n_chunk + 2) + r'}{l}{\textbf{Planetary Parameters}:} \\' + '\n')
            for param_name, median_key, labels in planetary_parameter_specs:
                if _normalize_param_name(param_name) in selected_parameters:
                    write(_row_slice(labels, idx), fout)

            # conclude with the closing rules; the flushleft notes block is only used in the final table
            fout.write(r'\hline' + '\n' +
                       r'\end{tabular}' + '\n')
            if is_last and show_priors:
                fout.write(r'\begin{flushleft}' + '\n' +
                           r'\textbf{Notes:} The priors for each system are labeled as $\mathcal{G}$[mean, standard deviation] if they are Gaussian priors and $\mathcal{U}$[lower limit, upper limit] if they are uniform priors.' + '\n' +
                           r'\end{flushleft}' + '\n')
            fout.write(r'\end{table*}')
            if not is_last:
                # keep every piece except the last on the same table number
                fout.write('\n' + r'\addtocounter{table}{-1}')

    for k, fname in enumerate(filenames):
        _write_chunk(fname, chunks[k], is_first=(k == 0), is_last=(k == n_chunks - 1))

def _extract_grid_rows(html, heading_text):
    """Extract the JavaScript row-data array for the grid beneath a named section heading.

    Parameters
    ----------
    html : str
        The full HTML page content.
    heading_text : str
        The visible heading text to locate, such as "Time Series Observations".

    Returns
    -------
    list[dict]
        A list of row dictionaries parsed from the embedded JavaScript array.
    """
    heading_match = re.search(
        r'<div[^>]*class=["\']grid_header["\'][^>]*>\s*' + re.escape(heading_text) + r'\b',
        html,
        flags=re.IGNORECASE | re.DOTALL,
    )
    if not heading_match:
        raise RuntimeError(f"Could not find section: {heading_text}")

    # Search only the content that appears after the matching heading.
    after_heading = html[heading_match.end():]
    row_match = re.search(
        r'var\s+(rowData\d+)\s*=\s*(\[[\s\S]*?\]);\s*//\s*Grid options',
        after_heading,
        flags=re.IGNORECASE | re.DOTALL,
    )
    if not row_match:
        raise RuntimeError(f"Could not find row data after section: {heading_text}")

    # Convert the JavaScript literals to Python literals before parsing.
    data = (
        row_match.group(2)
        .replace("true", "True")
        .replace("false", "False")
        .replace("null", "None")
    )

    return ast.literal_eval(data)


def get_followup_table(tic_id):
    """Fetch the follow-up observations for a given TIC ID and return a cleaned DataFrame.

    Parameters
    ----------
    tic_id : str
        TESS Input Catalog identifier, with or without the leading "TIC " prefix.

    Returns
    -------
    pandas.DataFrame
        A cleaned table with telescope, date, camera, filter, and size metadata.
    """
    if tic_id.startswith('TIC '):
        tic_id = tic_id.replace('TIC ', '')
    url = "https://exofop.ipac.caltech.edu/tess/target.php?id=" + tic_id
    with urlopen(url, timeout=20) as response:
        html = response.read().decode('utf-8', 'ignore')

    rows = _extract_grid_rows(html, 'Time Series Observations')

    df = pd.DataFrame(rows)
    df_short = df[['tstel', 'tsdate', 'tscam', 'tsfilt', 'tspix', 'tspsf', 'tspar']].copy()
    df_short.columns = ['Telescope', 'Date', 'Camera', 'Filter', r'Pix. Scale ($\arcsec$/pix)', r'PSF FWHM ($\arcsec$)', r'Aper. Rad. ($\arcsec$)']

    # Extract telescope size (m) and remove it from the Telescope column.
    df_short['Tel. Size (m)'] = df_short['Telescope'].str.extract(r'\((\d*\.?\d*)\s*m\)', expand=False).astype(float)
    df_short['Telescope'] = df_short['Telescope'].str.replace(r'\s*\(\d*\.?\d*\s*m\)', '', regex=True).str.strip()

    df_short = df_short.sort_values(by='Date', ascending=True).reset_index(drop=True)
    df_short = df_short[['Telescope', 'Tel. Size (m)', 'Date', 'Camera', 
                         'Filter', r'Pix. Scale ($\arcsec$/pix)', r'PSF FWHM ($\arcsec$)', r'Aper. Rad. ($\arcsec$)']]

    # reformat date to Year Mon Day format
    df_short['Date'] = pd.to_datetime(df_short['Date'], format='%Y-%m-%d')
    df_short['Date'] = df_short['Date'].dt.strftime('%Y %b %d')
    
    for filter in df_short['Filter']:
        if filter is not None:
            # Escape special LaTeX characters in the filter names
            escaped_filter = filter.replace('#', r'\#').replace('_', r'\_')
            df_short.loc[df_short['Filter'] == filter, 'Filter'] = escaped_filter
    
    # Truncate trailing zeros and convert floats to strings for LaTeX formatting
    float_columns = ['Tel. Size (m)', r'Pix. Scale ($\arcsec$/pix)', r'PSF FWHM ($\arcsec$)', r'Aper. Rad. ($\arcsec$)']
    for col in float_columns:
        df_short[col] = df_short[col].apply(lambda x: ('{:.3f}'.format(x)).rstrip('0').rstrip('.') if pd.notnull(x) else x)

    # Replace NaNs with '---'
    df_short = df_short.fillna('---')

    return df_short


def generate_master_followup_table(tic_list, toi_list):
    """Build a combined follow-up table for many TIC/TOI pairs.

    Parameters
    ----------
    tic_list : list[str]
        TIC identifiers, with or without the "TIC " prefix.
    toi_list : list[str]
        TOI identifiers, with or without the "TOI-" prefix.

    Returns
    -------
    pandas.DataFrame
        A concatenated table with TIC and TOI columns added to each row group.
    """
    master_df = pd.DataFrame()

    for tic_id, toi_id in zip(tic_list, toi_list):
        if not tic_id.startswith('TIC '):
            tic_id = 'TIC ' + tic_id

        if toi_id.startswith('TOI-'):
            toi_id = toi_id.replace('TOI-', '')
        try:
            df = get_followup_table(tic_id)
        except Exception as e:
            print(f"Error occurred while fetching follow-up table for TOI-{toi_id}: {e}")
            continue

        # Add TIC and TOI ids to the first row for this target.
        df['TIC ID'] = tic_id.replace('TIC ', '')
        df.loc[1:, 'TIC ID'] = ''
        df['TOI Number'] = toi_id
        df.loc[1:, 'TOI Number'] = ''

        master_df = pd.concat([master_df, df], ignore_index=True)
        master_df = master_df[['TIC ID', 'TOI Number', 'Telescope', 'Tel. Size (m)', 'Date', 'Camera',
                                'Filter', r'Pix. Scale ($\arcsec$/pix)', r'PSF FWHM ($\arcsec$)', r'Aper. Rad. ($\arcsec$)']]
    return master_df

def convert_table_to_latex_and_save(df, filename):
    """Convert a DataFrame to LaTeX format and save it to a file.

    Parameters
    ----------
    df : pandas.DataFrame
        The DataFrame to convert.
    filename : str
        The name of the file to save the LaTeX output.
    """
    latex_str = df.to_latex(index=False, escape=False, caption='Summary of Follow-up Observations', label='tab:followup', longtable=True)

    # required for longtable to work in a two-column document
    latex_str = '\\onecolumn\n' + latex_str + '\n\\twocolumn'
    
    with open(filename, 'w') as f:
        f.write(latex_str)

def generate_followup_table(tic_list, toi_list, output_filename):
    """Generate the follow-up table and save it as a LaTeX file.

    Parameters
    ----------
    tic_list : list[str]
        List of TIC identifiers.
    toi_list : list[str]
        List of TOI identifiers.
    output_filename : str
        The filename for the output LaTeX file.
    """
    master_df = generate_master_followup_table(tic_list, toi_list)
    convert_table_to_latex_and_save(master_df, output_filename)

def get_hri_table(tic_id):
    """Fetch the high-resolution imaging observations for a given TIC ID and return a cleaned DataFrame.

    Parameters
    ----------
    tic_id : str
        TESS Input Catalog identifier, with or without the leading "TIC " prefix.

    Returns
    -------
    pandas.DataFrame
        A cleaned table with telescope, date, camera, filter, and size metadata.
    """
    if tic_id.startswith('TIC '):
        tic_id = tic_id.replace('TIC ', '')
    url = "https://exofop.ipac.caltech.edu/tess/target.php?id=" + tic_id
    with urlopen(url, timeout=20) as response:
        html = response.read().decode('utf-8', 'ignore')

    rows = _extract_grid_rows(html, 'Imaging Observations')

    df = pd.DataFrame(rows)
    df_short = df[['itel', 'idate', 'iinst', 'itype', 'ifilt', 'ipix', 'ipsf', 'icont']].copy()
    df_short.columns = ['Telescope', 'Date', 'Instrument', 'Imaging Type', 'Filter', r'Pix. Scale ($\arcsec$/pix)', r'PSF FWHM ($\arcsec$)', 'Contrast']

    df_short = df_short.sort_values(by='Date', ascending=True).reset_index(drop=True)

    for contrast in df_short['Contrast']:
        if contrast is not None:
            # Replace 'delta' with the LaTeX delta symbol
            new_contrast = contrast.replace('delta', r'$\Delta$')
            df_short.loc[df_short['Contrast'] == contrast, 'Contrast'] = new_contrast

    for filter in df_short['Filter']:
        if filter is not None:
            # Escape special LaTeX characters in the filter names
            escaped_filter = filter.replace('#', r'\#').replace('_', r'\_')
            df_short.loc[df_short['Filter'] == filter, 'Filter'] = escaped_filter

    # reformat date to Year Mon Day format
    df_short['Date'] = pd.to_datetime(df_short['Date'], format='%Y-%m-%d')
    df_short['Date'] = df_short['Date'].dt.strftime('%Y %b %d')

    # Truncate trailing zeros and convert floats to strings for LaTeX formatting
    float_columns = [r'Pix. Scale ($\arcsec$/pix)', r'PSF FWHM ($\arcsec$)']
    for col in float_columns:
        df_short[col] = df_short[col].apply(lambda x: ('{:.3f}'.format(x)).rstrip('0').rstrip('.') if pd.notnull(x) else x)

    # Replace NaNs and empty strings with '---'
    df_short = df_short.fillna('---')
    df_short = df_short.replace('', '---')

    return df_short


def generate_master_table(tic_list, toi_list):
    """Build a combined HRI table for many TIC/TOI pairs.

    Parameters
    ----------
    tic_list : list[str]
        TIC identifiers, with or without the "TIC " prefix.
    toi_list : list[str]
        TOI identifiers, with or without the "TOI-" prefix.

    Returns
    -------
    pandas.DataFrame
        A concatenated table with TIC and TOI columns added to each row group.
    """
    master_df = pd.DataFrame()

    for tic_id, toi_id in zip(tic_list, toi_list):
        if not tic_id.startswith('TIC '):
            tic_id = 'TIC ' + tic_id

        if toi_id.startswith('TOI-'):
            toi_id = toi_id.replace('TOI-', '')
        try:
            df = get_hri_table(tic_id)
        except Exception as e:
            print(f"Error occurred while fetching HRI table for TOI-{toi_id}: {e}")
            continue

        # Add TIC and TOI ids to the first row for this target.
        df['TIC ID'] = tic_id.replace('TIC ', '')
        df.loc[1:, 'TIC ID'] = ''
        df['TOI Number'] = toi_id
        df.loc[1:, 'TOI Number'] = ''

        master_df = pd.concat([master_df, df], ignore_index=True)
        master_df = master_df[['TIC ID', 'TOI Number', 'Telescope', 'Date', 'Instrument', 'Imaging Type', 'Filter', 
                               r'Pix. Scale ($\arcsec$/pix)', r'PSF FWHM ($\arcsec$)', 'Contrast']]
    return master_df

def generate_hri_table(tic_list, toi_list, output_filename):
    """Generate the HRI table and save it as a LaTeX file.

    Parameters
    ----------
    tic_list : list[str]
        List of TIC identifiers.
    toi_list : list[str]
        List of TOI identifiers.
    output_filename : str
        The filename for the output LaTeX file.
    """
    master_df = generate_master_table(tic_list, toi_list)
    convert_table_to_latex_and_save(master_df, output_filename)

# EXOFASTv2 SED bandnames that map onto the Gaia magnitude rows of the secondary star table
SECONDARY_GAIA_BANDS = {
    'G': ['Gaia_G_EDR3', 'Gaia'],
    'BP': ['Gaia_BP_EDR3', 'GaiaBP'],
    'RP': ['Gaia_RP_EDR3', 'GaiaRP'],
}

# EXOFASTv2 SED bandnames that map onto the AO/speckle contrast rows of the secondary star table
SECONDARY_CONTRAST_BANDS = {
    'I': ['I', 'Ic', 'Icont'],
    'J': ['J', 'J2M', 'Jcont'],
    'H': ['H', 'H2M', 'Hcont'],
    'K': ['K', 'K2M', 'Ks', 'Kcont', 'Kshort'],
}


def parse_sed_stars(sed_path):
    """Parse an EXOFASTv2 SED file into the magnitudes belonging to each star.

    The final column of an SED file is the star index. An index of ``1`` means the magnitude
    belongs to star 1 alone, ``0,1`` means it is the blended magnitude of stars 0 and 1, and
    ``1-0`` means it is a magnitude difference (a contrast) between stars 1 and 0.

    Parameters
    ----------
    sed_path : str
        Path to the EXOFASTv2 .sed file.

    Returns
    -------
    dict
        ``{star_index: {'mags': {bandname: (value, error)}, 'deltamags': {bandname: (value, error)}}}``
        for every star index that appears in the file.
    """
    columns = ['bandname', 'magnitude', 'used_errors', 'catalog_errors', 'star_index']
    sedtable = pd.read_csv(sed_path, sep=r'\s+', skiprows=1, header=None, names=columns,
                           comment='#', dtype=str)

    stars = {}
    for _, row in sedtable.iterrows():
        index = row['star_index']
        if not isinstance(index, str):
            continue # magnitudes without a star index belong to a single star fit
        index = index.strip()
        band = str(row['bandname']).strip()
        try:
            value = float(row['magnitude'])
            error = float(row['used_errors'])
        except (TypeError, ValueError):
            continue

        if '-' in index: # a contrast, e.g. '1-0'
            star = int(index.split('-')[0])
            key = 'deltamags'
        elif ',' in index: # a blended magnitude, e.g. '0,1'. Not attributable to one star
            continue
        else:
            star = int(index)
            key = 'mags'

        stars.setdefault(star, {'mags': {}, 'deltamags': {}})[key][band] = (value, error)

    return stars


def get_stellar_companions(tic_id):
    """Fetch the stellar companions detected in high-resolution imaging from ExoFOP.

    Parameters
    ----------
    tic_id : str or int
        TESS Input Catalog identifier, with or without the leading "TIC " prefix.

    Returns
    -------
    list[dict]
        One entry per ExoFOP companion detection, with the angular separation and its
        uncertainty (arcsec), position angle (degrees), filter, magnitude difference, and
        observation date. An empty list is returned if the target has no companions listed.
    """
    tic_id = str(tic_id)
    if tic_id.startswith('TIC '):
        tic_id = tic_id.replace('TIC ', '')
    url = "https://exofop.ipac.caltech.edu/tess/target.php?id=" + tic_id
    with urlopen(url, timeout=20) as response:
        html = response.read().decode('utf-8', 'ignore')

    try:
        rows = _extract_grid_rows(html, 'Stellar Companions')
    except RuntimeError:
        return []

    def _to_float(value):
        try:
            return float(value)
        except (TypeError, ValueError):
            return None

    detections = []
    for row in rows:
        sep = _to_float(row.get('sep'))
        if sep is None:
            continue # a detection without a separation is useless here
        detections.append({
            'sep': sep,
            'sep_err': _to_float(row.get('sep_e')),
            'pa': _to_float(row.get('pa')),
            # ExoFOP filter strings look like 'Ks: 2.15 (0.32) microns'
            'band': str(row.get('filt') or '').split(':')[0].strip(),
            'dmag': _to_float(row.get('dmag')),
            'date': row.get('odate'),
        })
    return detections


def _group_companion_detections(detections, sep_tol=0.1, pa_tol=5.0):
    """Group ExoFOP companion detections that refer to the same physical companion.

    ExoFOP lists one row per filter (and per observation) for each companion, so rows are
    grouped by their angular separation and position angle.

    Parameters
    ----------
    detections : list[dict]
        Detections from get_stellar_companions.
    sep_tol : float
        Separation tolerance, in arcsec, within which two detections are considered the same star.
    pa_tol : float
        Position angle tolerance, in degrees, within which two detections are considered the same star.

    Returns
    -------
    list[list[dict]]
        The detections, grouped by companion.
    """
    groups = []
    for detection in sorted(detections, key=lambda d: d['sep']):
        for group in groups:
            sep_ref = np.mean([d['sep'] for d in group])
            if abs(detection['sep'] - sep_ref) > sep_tol:
                continue
            pa_values = [d['pa'] for d in group if d['pa'] is not None]
            if (detection['pa'] is not None) and pa_values:
                pa_ref = np.mean(pa_values)
                # wrap the difference into [-180, 180] so companions near a PA of 0 still group
                pa_diff = abs((detection['pa'] - pa_ref + 180) % 360 - 180)
                if pa_diff > pa_tol:
                    continue
            group.append(detection)
            break
        else:
            groups.append([detection])
    return groups


def _match_companion_group(groups, deltamags, dmag_tol=0.05):
    """Identify which ExoFOP companion corresponds to a secondary star in an SED file.

    The contrasts in an SED file are taken from the same imaging observations that ExoFOP
    reports, so a secondary star is matched to the companion whose magnitude differences
    agree with the contrasts in the SED file.

    Parameters
    ----------
    groups : list[list[dict]]
        Companion detections grouped by _group_companion_detections.
    deltamags : dict
        ``{bandname: (value, error)}`` contrasts for the secondary star, from the SED file.
    dmag_tol : float
        Tolerance, in magnitudes, for calling an SED contrast and an ExoFOP contrast the same measurement.

    Returns
    -------
    tuple
        ``(group, matched_detections, ambiguous)``. The group is None if no companion could be
        matched, and ambiguous is True if more than one companion matched equally well.
    """
    best_group = None
    best_matches = []
    best_score = 0
    ambiguous = False

    for group in groups:
        matches = []
        for value, _ in deltamags.values():
            for detection in group:
                if (detection['dmag'] is not None) and (abs(detection['dmag'] - value) <= dmag_tol):
                    matches.append(detection)
                    break
        if len(matches) > best_score:
            best_group, best_matches, best_score = group, matches, len(matches)
            ambiguous = False
        elif (len(matches) == best_score) and (best_score > 0):
            ambiguous = True

    if best_score == 0:
        # nothing to match against (or no contrasts in the SED file): a lone companion is unambiguous
        if len(groups) == 1:
            return groups[0], groups[0], False
        return None, [], len(groups) > 1

    return best_group, best_matches, ambiguous


def _companion_separation(detections):
    """Combine the separations of several detections of one companion into a value and uncertainty.

    Parameters
    ----------
    detections : list[dict]
        The detections of a single companion.

    Returns
    -------
    tuple
        The mean angular separation (arcsec) and its uncertainty, which combines the reported
        measurement uncertainties with the scatter between detections. The uncertainty is None
        if neither is available.
    """
    separations = np.array([d['sep'] for d in detections], dtype=float)
    sep = float(np.mean(separations))

    scatter = float(np.std(separations)) if len(separations) > 1 else 0.0
    reported = [d['sep_err'] for d in detections if d['sep_err'] is not None]
    measurement = float(np.mean(reported)) if reported else 0.0

    sep_err = np.sqrt(scatter**2 + measurement**2)
    return sep, (sep_err if sep_err > 0 else None)


def _round_with_errors(value, up_err, low_err, num_sig_figs=2):
    """Round a value and its uncertainties so the uncertainties carry a set number of significant figures.

    The value and both uncertainties are written to the same number of decimal places, which is set
    by the smaller of the two uncertainties.

    Parameters
    ----------
    value : float
        The value to round.
    up_err : float
        The upper uncertainty on the value.
    low_err : float
        The lower uncertainty on the value.
    num_sig_figs : int
        The number of significant figures to keep in the uncertainties.

    Returns
    -------
    tuple[str, str, str]
        The rounded value and uncertainties, as strings.
    """
    def decimal_places(err):
        # the number of decimal places left once the uncertainty is rounded to num_sig_figs
        # significant figures. Trailing zeros are dropped along the way, so an uncertainty
        # reported to one significant figure does not gain a digit it never had.
        if (err is None) or (not np.isfinite(err)) or (err == 0):
            return None
        err_str = remove_sci_notation(float(round_sig_figs(err, num_sig_figs)))
        return len(err_str.split('.')[1]) if '.' in err_str else 0

    places = [p for p in (decimal_places(up_err), decimal_places(low_err)) if p is not None]
    decimals = max(places) if places else 3 # fall back to the 3 decimals used for magnitudes elsewhere

    def as_string(x):
        return '{:.{p}f}'.format(np.round(x, decimals), p=max(decimals, 0))

    return as_string(value), as_string(up_err), as_string(low_err)


def gen_secondary_str(array, value, up_err=None, low_err=None):
    """Append a value and its uncertainties to a row of the secondary star table.

    Parameters
    ----------
    array : list
        The list of strings corresponding to the table row.
    value : float or None
        The value to append. None (or a non-finite value) writes the '---' filler.
    up_err : float or None
        The upper uncertainty on the value. If None, only the value is written.
    low_err : float or None
        The lower uncertainty on the value. Defaults to the upper uncertainty.
    """
    if (value is None) or (not np.isfinite(float(value))):
        array.append('& --- ')
        return

    if low_err is None:
        low_err = up_err
    if up_err is None:
        up_err = low_err

    if (up_err is None) or (not np.isfinite(float(up_err))) or (not np.isfinite(float(low_err))):
        # no reliable errors, so the value is written on its own
        array.append('& $' + remove_sci_notation(np.round(float(value), 3)) + '$ ')
        return

    val_str, up_str, low_str = _round_with_errors(float(value), float(up_err), float(low_err))
    if up_str == low_str:
        array.append('& $' + val_str + r' \pm ' + up_str + '$ ')
    else:
        array.append('& $' + val_str + r'^{+' + up_str + '}_{-' + low_str + '}$ ')


def _is_unbound(star_type):
    """Whether a classification describes a star that is not bound to the target star.

    An unbound star, such as an unrelated background star, does not share the distance of the
    system, so its angular separation cannot be turned into a projected separation.

    Parameters
    ----------
    star_type : str or None
        A classification from the star_types argument of secondary_stars_table(), such as
        'Bound companion' or 'Background star'. None counts as unclassified, and so as bound.
    """
    if star_type is None:
        return False
    return any(word in str(star_type).lower() for word in ('background', 'unbound'))


def secondary_stars_table(target_list, path, file_prefix, tic_list=None, host_list=None, outputpath='.',
                          distance_source='median', distances_external=None,
                          angular_separations=None, star_types=None,
                          sep_tol=0.1, pa_tol=5.0, dmag_tol=0.05, max_stars_per_table=5):
    '''
    Generates a table of the secondary stars in a set of EXOFASTv2 fits, with one column per
    secondary star. Angular separations are taken from the stellar companions reported on ExoFOP,
    projected separations are calculated from those separations and the distance to the system, and
    the magnitudes and AO/speckle contrasts are collected from the EXOFASTv2 SED files.

    Secondary stars are identified by the star index in the final column of the SED file: an index
    of '1' marks a magnitude belonging to star 1, and an index of '1-0' marks a contrast between
    star 1 and the target star. Each secondary star is matched to an ExoFOP companion by comparing
    its contrasts against the magnitude differences reported there, since both come from the same
    imaging observations.

    Every input list holds one entry per secondary star, so a host with two secondary stars is
    listed twice. Secondary stars that share a file prefix are taken in the order they appear in
    the SED file, meaning the first is star 1, the second is star 2, and so on.

    Parameters
    -----------
    target_list: an array of strings containing the name of each secondary star. Ex: ['TOI-3988 B', 'TIC 123456789'],
        where a bound companion is nominally named after its host and an unbound background star is named by its TIC ID
    path: path of EXOFASTv2 output files
    file_prefix: an array of strings containing the file prefix used in the EXOFASTv2 fit of each secondary star's host
    tic_list: optional array of the TIC IDs of the hosts, used to look up their companions on ExoFOP. If None, the file
        prefixes are used, since fits are usually named after the TIC ID of the host
    host_list: optional array of the names of the hosts. If given, they are written in a 'Planet Host' row beneath the star
        names and the classification row
    outputpath: the folder in which the table should be generated. Current working directory by default
    distance_source: accepts 'median' to use the distance from the EXOFASTv2 fit, 'gaia' to invert the Gaia DR3 parallax,
        or 'external' to provide distances directly
    distances_external: N by 2 array containing the distance in pc (first column) and its uncertainty (second column) for
        each secondary star, used when distance_source='external'
    angular_separations: optional dict of manual angular separations, in arcsec, used in place of the ExoFOP values. Keys are
        the names in target_list, and values may be a separation or a (separation, uncertainty) pair
    star_types: optional dict of classifications describing what each secondary star is, such as 'Bound companion' or
        'Background star'. Keys are the names in target_list. The classification row is only written if at least one
        secondary star has been classified. A classification naming the star as a background or unbound star suppresses
        its projected separation, since such a star does not share the distance of the system
    sep_tol: separation tolerance, in arcsec, within which two ExoFOP detections are treated as the same companion
    pa_tol: position angle tolerance, in degrees, within which two ExoFOP detections are treated as the same companion
    dmag_tol: tolerance, in magnitudes, for matching an SED contrast to an ExoFOP magnitude difference
    max_stars_per_table: the maximum number of secondary stars shown in a single table. If the number of secondary stars exceeds
        this value, the stars are split across multiple secondary_stars_table .tex files. Every table after the first is captioned
        "\\textit{(Continued)}" and every table except the last gets "\\addtocounter{table}{-1}" so that all pieces share one table
        number. The \\begin{minipage} notes block is only written in the last table. Set to None (or 0) to force a single table.
    '''

    # Setting up to save the table as a .tex file

    if os.path.exists(outputpath) == False:
        os.mkdir(outputpath)

    # Falling back on the file prefixes for the TIC IDs, since fits are usually named after the TIC ID of the host

    if tic_list is None:
        if all(str(prefix).strip().isdigit() for prefix in file_prefix):
            tic_list = [str(prefix).strip() for prefix in file_prefix]
        else:
            raise ValueError('tic_list is needed when the file prefixes are not the TIC IDs of the hosts.')

    # initializing rows
    sep_arr = [r'$\rho$ & Angular separation ($\arcsec$) ']
    projsep_arr = [r'$\rho_{\rm proj}$ & Projected separation (AU) ']
    gaia_g_arr = [r'${\rm G}$ & Gaia $G$ mag. ']
    gaia_bp_arr = [r'$G_{\rm BP}$ & Gaia $G_{\rm BP}$ mag. ']
    gaia_rp_arr = [r'$G_{\rm RP}$ & Gaia $G_{\rm RP}$ mag. ']
    contrast_arrs = {
        'I': [r'$\Delta I$ & $I$-band contrast (mag) '],
        'J': [r'$\Delta J$ & $J$-band contrast (mag) '],
        'H': [r'$\Delta H$ & $H$-band contrast (mag) '],
        'K': [r'$\Delta K$ & $K$-band contrast (mag) '],
    }

    star_labels = [] # the column headers, one per secondary star
    host_labels = [] # the host of each secondary star, so the columns can be traced back to a system
    type_labels = [] # what each secondary star is, if it has been classified

    # a host with more than one secondary star appears more than once, so its files are only read,
    # and its ExoFOP page only queried, the first time it comes up
    sed_cache = {}
    median_cache = {}
    parallax_cache = {}
    companion_cache = {}

    for i in range(len(target_list)):
        label = target_list[i]
        prefix = file_prefix[i]
        tic = str(tic_list[i])

        # collecting the magnitudes of every star in the fit
        if prefix not in sed_cache:
            sed_cache[prefix] = parse_sed_stars(path + prefix + '.sed')
        stars = sed_cache[prefix]
        secondary_indices = sorted(index for index in stars if index > 0)

        # secondary stars that share a file prefix are taken in the order they are listed in target_list
        occurrence = list(file_prefix[:i]).count(prefix)
        if len(secondary_indices) == 0:
            warnings.warn(f'There are no secondary stars in the SED file of {prefix}. Skipping {label}.')
            continue
        if occurrence >= len(secondary_indices):
            warnings.warn(f'More secondary stars are listed for {prefix} than its SED file holds. Skipping {label}.')
            continue
        star = secondary_indices[occurrence]
        mags = stars[star]['mags']
        deltamags = stars[star]['deltamags']

        star_labels.append(label)
        host_labels.append(str(host_list[i]) if host_list is not None else '---')

        # what the secondary star is, e.g. a bound companion or an unrelated background star
        star_type = star_types.get(label) if star_types else None
        type_labels.append(str(star_type) if star_type is not None else '---')
        unbound = _is_unbound(star_type)

        # collecting the distance, which turns the angular separation into a projected separation
        distance = distance_up_err = distance_low_err = None
        if unbound:
            pass # an unbound star lies at its own distance, so it is given no projected separation
        elif distance_source == 'median':
            if prefix not in median_cache:
                median_cache[prefix] = grab_medians(path=path, file_prefix=prefix)
            medians = median_cache[prefix]
            if medians.parname.isin(['distance_0']).any():
                distance = float(medians.median_value[medians.parname == 'distance_0'].iloc[0])
                distance_up_err = float(medians.upper_error[medians.parname == 'distance_0'].iloc[0])
                distance_low_err = float(medians.lower_error[medians.parname == 'distance_0'].iloc[0])
            else:
                warnings.warn(f'There is no distance in the median table of {prefix}. The projected separation of {label} will be omitted.')
        elif distance_source == 'gaia':
            if tic not in parallax_cache:
                vgaia = Vizier(columns=['_r', 'Plx', 'e_Plx'], catalog='I/355/gaiadr3')
                data_gaia = vgaia.query_region('TIC ' + tic, radius=Angle(6, "arcsec"))
                data_gaia = data_gaia[0]
                data_gaia.sort('_r') # sort by distance from the target star
                parallax_cache[tic] = (float(data_gaia['Plx'][0]), float(data_gaia['e_Plx'][0]))
            parallax, parallax_err = parallax_cache[tic]
            distance = 1000 / parallax
            distance_up_err = distance_low_err = distance * parallax_err / parallax
        elif distance_source == 'external':
            distance = float(distances_external[i][0])
            distance_up_err = distance_low_err = float(distances_external[i][1])

        # collecting the companions that ExoFOP reports for the host
        if tic not in companion_cache:
            try:
                detections = get_stellar_companions(tic)
            except Exception as e:
                warnings.warn(f'Could not fetch the ExoFOP companions of TIC {tic}: {e}')
                detections = []
            companion_cache[tic] = _group_companion_detections(detections, sep_tol=sep_tol, pa_tol=pa_tol)
        groups = companion_cache[tic]

        # the angular separation, either provided by hand or matched to an ExoFOP companion
        sep = sep_err = None
        manual = angular_separations.get(label) if angular_separations else None
        if manual is not None:
            if np.ndim(manual) == 0:
                sep, sep_err = float(manual), None
            else:
                sep, sep_err = float(manual[0]), float(manual[1])
        elif len(groups) > 0:
            group, matches, ambiguous = _match_companion_group(groups, deltamags, dmag_tol=dmag_tol)
            if group is None:
                warnings.warn(f'Could not match {label} to an ExoFOP companion. Its separation will be omitted.')
            elif ambiguous:
                warnings.warn(f'{label} matches more than one ExoFOP companion equally well. Its separation will be omitted.')
            else:
                sep, sep_err = _companion_separation(matches)
        else:
            warnings.warn(f'ExoFOP lists no companions for TIC {tic}. The separation of {label} will be omitted.')

        gen_secondary_str(sep_arr, sep, sep_err)

        # the projected separation, in AU, from the angular separation and the distance
        if unbound or (sep is None) or (distance is None):
            gen_secondary_str(projsep_arr, None)
        else:
            projsep = sep * distance
            # add the fractional uncertainties of the separation and the distance in quadrature
            frac_sep = (sep_err / sep) if sep_err else 0.0
            projsep_up_err = projsep * np.sqrt(frac_sep**2 + (distance_up_err / distance)**2)
            projsep_low_err = projsep * np.sqrt(frac_sep**2 + (distance_low_err / distance)**2)
            gen_secondary_str(projsep_arr, projsep, projsep_up_err, projsep_low_err)

        # the Gaia magnitudes of the secondary star itself
        for band, array in (('G', gaia_g_arr), ('BP', gaia_bp_arr), ('RP', gaia_rp_arr)):
            value = error = None
            for bandname in SECONDARY_GAIA_BANDS[band]:
                if bandname in mags:
                    value, error = mags[bandname]
                    break
            gen_secondary_str(array, value, error)

        # the AO/speckle contrasts between the secondary star and the target star
        for band, array in contrast_arrs.items():
            value = error = None
            for bandname in SECONDARY_CONTRAST_BANDS[band]:
                if bandname in deltamags:
                    value, error = deltamags[bandname]
                    break
            gen_secondary_str(array, value, error)

    n_stars = len(star_labels)
    if n_stars == 0:
        warnings.warn('No secondary stars were found for any target. No table was generated.')
        return

    hosts_named = any(host != '---' for host in host_labels)
    classified = any(star_type != '---' for star_type in type_labels)

    # only keep the contrast rows that at least one secondary star has a measurement in
    contrast_rows = [array for array in contrast_arrs.values()
                     if any(cell != '& --- ' for cell in array[1:])]

    # Generating the preamble

    preamble = (r'\providecommand{\tess}{\textit{TESS}\xspace}' + '\n')

    notes = (r'\vspace{2mm}' + '\n' +
             r'\begin{minipage}{\textwidth}' + '\n' +
             r'\footnotesize' + '\n' +
             r'\textbf{Notes:}' + '\n' +
             r'Angular separations and contrasts are reported on ExoFOP. Where a companion was detected in ' +
             r'several filters, the separations are averaged and their scatter is included in the uncertainty.\\' + '\n' +
             r'\end{minipage}' + '\n')

    # Deciding how to split the secondary stars across tables
    if not max_stars_per_table or max_stars_per_table >= n_stars:
        chunk_size = n_stars
    else:
        chunk_size = int(max_stars_per_table)
    n_chunks = max(1, int(np.ceil(n_stars / chunk_size)))

    # Deriving the output filename(s). The first table keeps the classic 'secondary_stars_table.tex'
    # name (bumping a numeric suffix if it already exists); continuation tables append '_2', '_3', ...
    first_name = 'secondary_stars_table.tex'
    suffix = 2
    while os.path.exists(f'{outputpath}/{first_name}'):
        first_name = f'secondary_stars_table_{suffix}.tex'
        suffix += 1
    stem = first_name[:-len('.tex')]
    filenames = [first_name] + [f'{stem}_{k}.tex' for k in range(2, n_chunks + 1)]
    print('Saving this table as ' + ', '.join(filenames) + '...')

    def _row_slice(arr, sl):
        '''Return a table row for the secondary stars in slice `sl`: the label cell(s) followed by the sliced star cells.'''
        return arr[:1] + arr[1:][sl]

    def _write_chunk(fname, sl, is_first, is_last):
        chunk_labels = star_labels[sl]
        n_chunk = len(chunk_labels)
        colstring = 'c' * n_chunk # one column per secondary star, after the two label columns
        namestring = ''.join(' & ' + label for label in chunk_labels)
        hoststring = ''.join(' & ' + host for host in host_labels[sl])
        typestring = ''.join(' & ' + star_type for star_type in type_labels[sl])

        caption = (r'\caption{Observed Properties the Secondary Stars}' if is_first
                   else r'\caption{\textit{(Continued)}}')

        with open(f'{outputpath}/{fname}', 'w') as fout:
            fout.write(preamble)
            fout.write(r'\begin{table*}' + '\n' +
                       r'\centering' + '\n' +
                       caption + '\n')
            if is_first:
                fout.write(r'\label{tab:secondary}' + '\n')
            fout.write(r'\scriptsize' + '\n')

            fout.write(r'\begin{tabular}{l l' + colstring + '}' + '\n' +
                       r'\hline' + '\n' +
                       r'& ' + namestring + r'\\' + '\n')
            # the classification is written first, so that the planet host is not read as one
            if classified:
                fout.write(r'& Classification' + typestring + r'\\' + '\n')
            if hosts_named:
                fout.write(r'& Planet Host' + hoststring + r'\\' + '\n')
            fout.write(r'\hline' + '\n')

            write(_row_slice(sep_arr, sl), fout)
            write(_row_slice(projsep_arr, sl), fout)
            write(_row_slice(gaia_g_arr, sl), fout)
            write(_row_slice(gaia_bp_arr, sl), fout)
            write(_row_slice(gaia_rp_arr, sl), fout)
            for array in contrast_rows:
                write(_row_slice(array, sl), fout)

            fout.write(r'\hline' + '\n' +
                       r'\end{tabular}' + '\n')
            if is_last:
                # the minipage notes block is only used in the final table
                fout.write(notes)
            fout.write(r'\end{table*}')
            if not is_last:
                # keep every piece except the last on the same table number
                fout.write('\n' + r'\addtocounter{table}{-1}')

    for k, fname in enumerate(filenames):
        sl = slice(k * chunk_size, min((k + 1) * chunk_size, n_stars))
        _write_chunk(fname, sl, is_first=(k == 0), is_last=(k == n_chunks - 1))
