"""Formatting and ExoFOP scraping helpers shared by the LaTeX table generators."""

import numpy as np
import pandas as pd
import re
import ast
import logging
import time
from http.client import HTTPException
from urllib.request import urlopen

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

def write(param_arr,file):
    for ii in param_arr:
        file.write(ii)
    file.write(r'\\'+'\n')

EXOFOP_TARGET_URL = 'https://exofop.ipac.caltech.edu/tess/target.php?id='


class ExoFOPFetchError(Exception):
    """A target's ExoFOP page could not be fetched in full.

    Deliberately not a RuntimeError, which is what _extract_grid_rows raises for a section a
    complete page does not have, so that callers treating a missing section as "nothing
    listed" do not mistake a failed fetch for one.
    """


def fetch_exofop_page(tic_id, attempts=3, delay=5):
    """Fetch a target's ExoFOP page, retrying when it comes back cut off.

    ExoFOP now and then answers with a page that stops just after its navigation bar: it has
    a success status but none of the target's data, and which targets it happens to varies
    from one request to the next. A complete page always ends in "</html>", so a page that
    does not is requested again, as is one whose request fails outright, up to `attempts`
    times with `delay` seconds between tries. Each retry is logged as a warning.

    Parameters
    ----------
    tic_id : str or int
        TESS Input Catalog identifier, with or without the leading "TIC " prefix.
    attempts : int
        How many times to request the page before giving up.
    delay : float
        Seconds to wait before each retry.

    Returns
    -------
    str
        The complete HTML of the page.

    Raises
    ------
    ExoFOPFetchError
        If no attempt returns the complete page.
    """
    tic_id = str(tic_id)
    if tic_id.startswith('TIC '):
        tic_id = tic_id.replace('TIC ', '')

    problem = None
    for attempt in range(1, attempts + 1):
        if attempt > 1:
            logging.warning(f'Retrying the ExoFOP page of TIC {tic_id} ({problem}); attempt {attempt} of {attempts}.')
            time.sleep(delay)
        try:
            with urlopen(EXOFOP_TARGET_URL + tic_id, timeout=20) as response:
                html = response.read().decode('utf-8', 'ignore')
        except (OSError, HTTPException) as e:
            problem = f'the request failed: {e}'
            continue
        if html.rstrip().lower().endswith('</html>'):
            return html
        problem = f'the page was cut off after {len(html)} characters'

    raise ExoFOPFetchError(f'Could not fetch the full ExoFOP page of TIC {tic_id} in {attempts} attempts; '
                           f'last, {problem}.')


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


# ExoFOP filter names mapped onto the notation used in the literature. Every standard filter
# is set in math mode so that it is italicised. Keys are matched case-sensitively, since
# lower-case single letters are Sloan bands while the upper-case ones are Johnson-Cousins
# (e.g. "i" is Sloan i-prime but "I" is Cousins I).
FOLLOWUP_FILTER_NAMES = {
    # Sloan, written with the conventional prime
    'up': r"$u'$", 'gp': r"$g'$", 'rp': r"$r'$", 'ip': r"$i'$", 'zp': r"$z'$",
    'u': r"$u'$", 'g': r"$g'$", 'r': r"$r'$", 'i': r"$i'$", 'z': r"$z'$",
    'SDSS-u': r"$u'$", 'SDSS-g': r"$g'$", 'SDSS-r': r"$r'$", 'SDSS-i': r"$i'$", 'SDSS-z': r"$z'$",
    'Sloan-u': r"$u'$", 'Sloan-g': r"$g'$", 'Sloan-r': r"$r'$", 'Sloan-i': r"$i'$", 'Sloan-z': r"$z'$",
    # EXOFASTv2's names for the Sloan bands, as its lightcurve filenames give them
    'Sloanu': r"$u'$", 'Sloang': r"$g'$", 'Sloanr': r"$r'$", 'Sloani': r"$i'$", 'Sloanz': r"$z'$",
    # z-short
    'zs': r'$z_s$', 'z_s': r'$z_s$', 'z-s': r'$z_s$', 'zshort': r'$z_s$', 'z-short': r'$z_s$',
    # Cousins
    'Rc': r'$R_c$', 'rc': r'$R_c$', 'Ic': r'$I_c$', 'ic': r'$I_c$',
    # Johnson-Cousins and near-infrared, italicised but otherwise unchanged
    'U': r'$U$', 'B': r'$B$', 'V': r'$V$', 'R': r'$R$', 'I': r'$I$',
    'J': r'$J$', 'H': r'$H$', 'K': r'$K$', 'Ks': r'$K_s$', 'K_s': r'$K_s$',
}

# separators that join several filters into one ExoFOP entry, e.g. "gp-ip" or "g, r, i, z_s"
_FILTER_SEPARATOR = re.compile(r'(\s*[,+/&-]\s*)')


def format_filter_name(filter_name):
    """Convert an ExoFOP filter name into the notation used in the literature.

    Names listed in FOLLOWUP_FILTER_NAMES are translated directly ("ip" becomes "$i'$"),
    every standard filter being set in math mode so that it comes out italicised. Any name
    mentioning a clear filter is collapsed to "Clear", unitalicised, since the bandpasses
    ExoFOP quotes alongside it ("clear: 650 (500) nm") differ from observation to
    observation. An entry naming several filters at once is translated only when every
    piece of it is a known filter, so compound names such as "gp-ip" are converted while
    descriptive ones such as "g-narrow" are left alone. An entry that names a standard
    filter and then the instrument and bandpass behind a colon, as "R-ASTEP+: 850 (138) nm"
    does, keeps only the filter.

    Parameters
    ----------
    filter_name : str
        The filter name as it appears on ExoFOP.

    Returns
    -------
    str
        The converted name, or the original name with LaTeX special characters escaped
        when no conversion is known.
    """
    name = filter_name.strip()

    # Clear filters are quoted with a variety of bandpasses; they are all just "Clear".
    if 'clear' in name.lower():
        return 'Clear'

    if name in FOLLOWUP_FILTER_NAMES:
        return FOLLOWUP_FILTER_NAMES[name]

    # Split on the separators, keeping them, so the original spacing is preserved.
    pieces = _FILTER_SEPARATOR.split(name)
    filters = pieces[::2]
    if len(filters) > 1 and all(piece in FOLLOWUP_FILTER_NAMES for piece in filters):
        pieces[::2] = [FOLLOWUP_FILTER_NAMES[piece] for piece in filters]
        return ''.join(pieces)

    # An entry such as "R-ASTEP+: 850 (138) nm" names a standard filter and then the
    # instrument and bandpass it was measured with. Keep the filter and drop the rest, but
    # only when a colon marks that description off, so that narrow-band names like
    # "g-narrow" are not mistaken for the standard band they begin with.
    if ':' in name:
        band = name.split(':')[0].strip()
        if band in FOLLOWUP_FILTER_NAMES:
            return FOLLOWUP_FILTER_NAMES[band]
        # drop a trailing instrument name, as in "R-ASTEP+"
        band = band.split('-')[0].strip()
        if band in FOLLOWUP_FILTER_NAMES:
            return FOLLOWUP_FILTER_NAMES[band]

    # Unrecognised name: escape the LaTeX special characters and leave it as it is.
    return filter_name.replace('#', r'\#').replace('_', r'\_')

def convert_table_to_latex_and_save(df, filename, caption='Summary of Follow-up Observations', label='tab:followup',
                                    fontsize=None):
    """Convert a DataFrame to LaTeX format and save it to a file.

    Parameters
    ----------
    df : pandas.DataFrame
        The DataFrame to convert.
    filename : str
        The name of the file to save the LaTeX output.
    caption : str
        The caption for the LaTeX table.
    label : str
        The label for the LaTeX table.
    fontsize : str, optional
        A LaTeX font size command, e.g. "\\scriptsize", written directly after \\onecolumn to
        shrink the table. The size stays in effect until the document resets it.
    """
    latex_str = df.to_latex(index=False, escape=False, caption=caption, label=label, longtable=True)

    # required for longtable to work in a two-column document
    preamble = '\\onecolumn\n'
    if fontsize is not None:
        preamble += fontsize + '\n'
    latex_str = preamble + latex_str + '\n\\twocolumn'
    
    with open(filename, 'w') as f:
        f.write(latex_str)
