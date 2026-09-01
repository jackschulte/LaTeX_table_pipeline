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
    priors = pd.read_csv(path + file_prefix + '.priors.final', sep=r'\s+', skiprows=1, header=None, comment='#', names=columns)

    for i in range(len(priors)):
        # find linked parameters and replace them with the first instance of the parameter
        if type(priors.meanvalue[i]) == str:
            if (priors.meanvalue[i] in priors.variable.values) or (priors.meanvalue[i].replace('_0', '') in priors.variable.values):
                priors.loc[i, 'meanvalue'] = priors.meanvalue[priors.variable == priors.meanvalue[i].replace('_0', '')].iloc[0]
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

        array.append('& $' + val_str)
        array.append(errstring + '$ ')
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
        if ('B' or 'C') in toi[-1]:
            ticid_A = TOI_df.loc[toi_id]['TIC ID']
            # Query the TIC v8.2 for the TIC identifier of the nearest sources
            tic_columns = ['_r', 'TIC']
            vtic = Vizier(columns=tic_columns, catalog='IV/39/tic82')
            data_tic = vtic.query_region('TIC ' + str(ticid_A), radius=Angle(6, "arcsec"))
        else:
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



def med_table(target_list, path, file_prefix_list, outputpath='.', bimodal=False, multistar=False, parameters=None):
    '''
    Generates a median table given the path to EXOFASTv2 output files.

    Parameters
    -----------
    target_list: an array of strings containing the names of each target. Names should match those in EXOFASTv2 output files.
    path: path of EXOFASTv2 output files
    file_prefix_list: list of prefixes used in EXOFASTv2 file generation
    outputpath: the folder in which the table should be generated. Current working directory by default
    parameters: optional list of parameter names to include in the table. Defaults to the set of parameters
        used in Schulte+ 2025
    '''

    def _normalize_param_name(name):
        return re.sub(r'[^a-z0-9]', '', str(name).lower())

    default_parameters = [
        'mstar', 'rstar', 'lstar', 'rhostar', 'logg', 'teff', 'feh', 'initfeh', 'age', 'eep', 'Av',
        'distance', 'Period', 'rp', 'mp', 'tc', 't0', 'a', 'ideg', 'e', 'omegadeg', 'teq', 'tcirc',
        'k', 'p', 'ar', 'depth_TESS', 'tau', 't14', 'b', 'rhop', 'loggp', 'q', 'dr'
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

    newfile = 'median_table.tex'

    # if this file exists, come up with a new name
    i = 2
    while os.path.exists(f'{outputpath}/{newfile}'):
        newfile = 'median_table_' + str(i) + '.tex'
        i += 1
    print(f'Saving this table as {newfile}...')

    for ii in range(len(target_list)):
        medians = grab_medians(path=path, file_prefix=file_prefix_list[ii], bimodal=bimodal)

        for param_name, median_key, labels in parameter_specs:
            if _normalize_param_name(param_name) in selected_parameters or _normalize_param_name(median_key) in selected_parameters:
                make_median_string(medians, median_key, labels)

    colstring = 'lc'
    namestring = ''
    
    for ii in range(len(target_list)):
        colstring+='c'
        namestring += (' & ' + target_list[ii])

    # Collecting priors to put at the top of the table
    parallax_prior = '' # initializing strings
    metallicity_prior = ''
    extinction_prior = ''
    dilution_prior = ''
    if bimodal == False:
        dilute_bool = np.zeros_like(target_list) # to keep track of which targets were fit for dilution
        for ii in range(len(target_list)):
            priortable = grab_priors(file_prefix_list[ii], path)
            parallax_prior_mean = priortable.meanvalue[priortable.variable == 'parallax'].iloc[0]
            parallax_prior_stdev = priortable.stdev[priortable.variable == 'parallax'].iloc[0]
            parallax_prior += (r'& $\mathcal{G}$[' + round_sig_figs(parallax_prior_mean, 5) + r', ' + round_sig_figs(parallax_prior_stdev, 5) + r'] ')
            metallicity_prior_mean = priortable.meanvalue[priortable.variable == 'feh'].iloc[0]
            metallicity_prior_stdev = priortable.stdev[priortable.variable == 'feh'].iloc[0]
            metallicity_prior += (r'& $\mathcal{G}$[' + round_sig_figs(metallicity_prior_mean, 5) + r', ' + round_sig_figs(metallicity_prior_stdev, 5) + r'] ')
            extinction_prior_upperbound = priortable.up_bound[priortable.variable == 'Av'].iloc[0]
            extinction_prior += (r'& $\mathcal{U}$[0, ' + round_sig_figs(extinction_prior_upperbound, 5) + r'] ')

            for x in priortable.variable: # finding the dilution term
                match = re.findall('dilute', x)
                if len(match) > 0:
                    dilute_colname = (match[0])
                    dilution_prior_mean = priortable.meanvalue[priortable.variable == dilute_colname].iloc[0]
                    dilution_prior_stdev = priortable.stdev[priortable.variable == dilute_colname].iloc[0]

                    dilute_bool[ii] = 1
            if dilute_bool[ii]:
                dilution_prior += (r'& $\mathcal{G}$[' + round_sig_figs(dilution_prior_mean, 5) + r', ' + remove_sci_notation(float(round_sig_figs(dilution_prior_stdev, 5))) + r'] ')
                # above line should be cleaned up in a future version. Maybe make a new function that removes scientific notation and sets sig figs for all numbers
            else:
                dilution_prior += (r'& --- ')

    # Generating the preamble
    
    with open(f'{outputpath}/{newfile}', 'w') as fout: 
        fout.write(r'\providecommand{\bjdtdb}{\ensuremath{\rm {BJD_{TDB}}}}'+'\n'+
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
    r'\providecommand{\tess}{\textit{TESS}\xspace}'+'\n'+
    r'\begin{table*}'+'\n'+
    r'\centering'+'\n'+
    r'\caption{Median Values and 68\% Confidence Intervals for Fitted Stellar and Planetary Parameters}'+'\n'+
    r'\label{tab:median}'+'\n'+
    r'\scriptsize'+'\n'+
    r'\begin{tabular}{ll' + colstring + '}'+'\n'+
    r'\hline' + '\n' +
    r'& ' + namestring + r'\\' +'\n'+
    r'\hline' + '\n' +
    r'\multicolumn{' + str(len(target_list) + 2) + r'}{l}{\textbf{Priors}:} \\' + '\n' +
    r'$\pi$ & Gaia Parallax (mas)' + parallax_prior + r'\\' + '\n' +
    r'$[{\rm Fe/H}]$ & Metallicity (dex)' + metallicity_prior + r'\\' + '\n' +
    r'$A_V$ & V-band extinction (mag)' + extinction_prior + r'\\' + '\n' + 
    r'$D_T$ & Dilution in \tess' + dilution_prior + r'\\' + '\n' +
    r'\hline' + '\n' +               
    r'\multicolumn{' + str(len(target_list) + 2) + r'}{l}{\textbf{Stellar Parameters}:} \\' + '\n' )
    #r'\smallskip\\\multicolumn{2}{l}{Stellar Parameters:}&\smallskip\\'+'\n')

            
        for param_name, median_key, labels in stellar_parameter_specs:
            if _normalize_param_name(param_name) in selected_parameters:
                write(labels, fout)

        fout.write(r'\multicolumn{' + str(len(target_list) + 2) + r'}{l}{\textbf{Planetary Parameters}:} \\' + '\n')
        for param_name, median_key, labels in planetary_parameter_specs:
            if _normalize_param_name(param_name) in selected_parameters:
                write(labels, fout)
        
        # conclude with \enddata at the bottom of the input .tex file
        fout.write(r'\hline' + '\n' + 
                   r'\end{tabular}' + '\n' +
                   r'\begin{flushleft}' + '\n' +
                   r'\textbf{Notes:} The priors for each system are labeled as $\mathcal{G}$[mean, standard deviation] if they are Gaussian priors and $\mathcal{U}$[lower limit, upper limit] if they are uniform priors.' + '\n' +
                   r'\end{flushleft}' + '\n' +
                   r'\end{table*}')

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