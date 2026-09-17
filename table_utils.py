"""Formatting and ExoFOP scraping helpers shared by the LaTeX table generators."""

import numpy as np
import pandas as pd
import re
import ast
import logging
import time
from decimal import Decimal, InvalidOperation, ROUND_HALF_UP, localcontext
from http.client import HTTPException
from urllib.request import urlopen

# every number is quantized inside this context, which is wide enough for the longest value a
# table holds (a time of conjunction written to eight decimal places) without the default
# 28-digit precision truncating it
_DECIMAL_PRECISION = 60


def to_decimal(x):
    '''
    Turns a number into a Decimal without losing the precision it was written with.

    A string keeps exactly the digits it carries, so a median EXOFASTv2 wrote as '0.04530' stays five
    decimal places wide; float('0.04530') prints as 0.0453 and would cost the median its last
    significant figure. A float is converted through repr(), the shortest decimal string that round
    trips, rather than through Decimal(float), which would expose the whole tail of the binary
    representation (0.1 becoming 0.1000000000000000055511151231257827).

    Parameters
    -----------
    x: the number to convert, as a string, float, int, or Decimal

    Raises
    -------
    InvalidOperation
        If x does not spell a number.
    '''
    if isinstance(x, Decimal):
        return x
    if isinstance(x, str):
        return Decimal(x.strip())
    if isinstance(x, (int, np.integer)):
        return Decimal(int(x))
    value = float(x)
    if not np.isfinite(value):
        raise InvalidOperation(f'{x!r} is not a finite number.')
    return Decimal(repr(value))

def is_finite_number(x):
    '''
    Whether a value spells a finite number, and so can be written into a table.

    Catches the None and masked entries a catalogue query returns for a missing measurement as well as
    the NaN and infinity a fit can produce, so that callers can write the '---' filler instead.

    Parameters
    -----------
    x: the value to test, of any type
    '''
    if x is None or isinstance(x, np.ma.core.MaskedConstant):
        return False
    try:
        to_decimal(x)
    except (InvalidOperation, ValueError, TypeError):
        return False
    return True

def decimal_places(x):
    '''
    Counts the digits a number carries after the decimal point when it is written positionally.

    This is the precision of a number as written rather than as valued, so it is only exact for a
    number passed as a string: '1.50' has two decimal places while the float 1.50 has one, since the
    trailing zero is gone by the time Python has parsed it.

    Parameters
    -----------
    x: the number to measure, as a string, float, int, or Decimal
    '''
    return max(-to_decimal(x).as_tuple().exponent, 0)

def format_to_decimals(x, decimals):
    '''
    Writes a number positionally with exactly the given number of digits after the decimal point.

    Trailing zeros are kept, so a value rounded to three decimal places always shows three of them,
    and half-way cases round away from zero (0.0625 to three decimals is 0.063) rather than to the
    nearest even digit as numpy and Python's own round() do.

    Parameters
    -----------
    x: the number to write, as a string, float, int, or Decimal
    decimals: the number of digits to write after the decimal point
    '''
    decimals = max(int(decimals), 0)
    with localcontext() as context:
        context.prec = _DECIMAL_PRECISION
        quantized = to_decimal(x).quantize(Decimal(1).scaleb(-decimals), rounding=ROUND_HALF_UP)
    if not quantized:
        quantized = quantized.copy_abs() # so that a small negative value is written as 0.00, not -0.00
    return format(quantized, 'f')

def format_sig_figs(x, num_sig_figs):
    '''
    Writes a number positionally to a fixed number of significant figures, keeping trailing zeros.

    Unlike round_sig_figs, which reports the shortest string that carries the rounded value, this keeps
    the zeros that state the precision: 6.2 to three significant figures is '6.20'. Use it where the
    result sets the precision that an uncertainty then has to match.

    Parameters
    -----------
    x: the number to write
    num_sig_figs: the number of significant figures to keep
    '''
    value = to_decimal(x)
    if not value:
        return format_to_decimals(value, num_sig_figs - 1)
    with localcontext() as context:
        context.prec = _DECIMAL_PRECISION
        exponent = value.adjusted() - num_sig_figs + 1
        rounded = value.quantize(Decimal(1).scaleb(exponent), rounding=ROUND_HALF_UP)
    return format(rounded, 'f')

def remove_sci_notation(x):
    '''
    Removes scientific notation from a number and returns it as a string.

    A number passed as a string keeps exactly the digits it carries, so that a value read out of a file
    is written with the precision it was recorded with. A float is written with the shortest digits
    that round trip, which drops any trailing zero the number no longer remembers.

    Parameters
    -----------
    x: the number to remove scientific notation from
    '''
    if isinstance(x, str):
        try:
            return format(to_decimal(x), 'f')
        except InvalidOperation:
            return x
    return np.format_float_positional(x, trim='-')

def round_sig_figs(x, num_sig_figs):
    '''
    Rounds a number to a specified number of significant figures.

    The result is always positional: '{:g}' alone switches to scientific notation for a number below
    1e-5 or above the requested number of figures, which would put an 'e-06' into the table.

    Trailing zeros are dropped, so 6.2 to three significant figures comes back as '6.2'. Where the
    precision has to be stated rather than merely carried, as for a value an uncertainty is written
    alongside, use format_sig_figs instead.

    Parameters
    -----------
    x: the number to round
    num_sig_figs: the number of significant figures to round to
    '''
    return remove_sci_notation(float('{:.{p}g}'.format(float(x), p=num_sig_figs)))

def _decimals_showing(decimals, *values):
    '''
    Widens a number of decimal places until every non-zero value given still shows a digit.

    A value that is small enough to round away to zero at the precision asked for says nothing at all:
    an uncertainty printed as '0.0' claims a measurement was exact. Where that would happen the
    precision is widened to the place of the value's leading significant figure, which is the fewest
    decimal places that keep it, and the widening is logged so that the mismatch is visible.

    Parameters
    -----------
    decimals: the number of decimal places asked for
    values: the values that have to survive being written at that precision
    '''
    widened = decimals
    for value in values:
        number = to_decimal(value)
        if not number:
            continue # an uncertainty of exactly zero is meant to be written as zero
        if number.copy_abs() < Decimal(5).scaleb(-widened - 1): # it would round to zero
            widened = -number.adjusted()
    if widened != decimals:
        logging.warning(f'{values} cannot be written at {decimals} decimal places without an '
                        f'uncertainty rounding away to zero, so {widened} are used instead. The value '
                        'was probably rounded more coarsely than its uncertainty.')
    return widened

def robust_decimal_errors(val, up_err, low_err, decimals=None):
    '''
    Writes a value and its uncertainties to the same number of decimal places.

    The value sets the precision: both uncertainties are rounded to the number of decimal places the
    value is written with, so that the three numbers line up and no trailing zero is dropped from any
    of them. Pass the numbers as strings, as they are read out of an EXOFASTv2 median file, to keep
    the precision EXOFASTv2 chose for them; passing floats loses every trailing zero before this
    function ever sees the number.

    If the value was written coarsely enough that an uncertainty would round away to zero, the
    precision is widened until the uncertainty keeps its leading digit and the value is padded to
    match, since an uncertainty of zero states something that was never measured.

    Parameters
    -----------
    val: the value to round
    up_err: the upper error on the value
    low_err: the lower error on the value
    decimals: the number of decimal places to write all three with. Taken from the precision of the
        value when it is None, which is what a median read as a string already carries.
    '''
    if decimals is None:
        decimals = decimal_places(val)
    decimals = _decimals_showing(decimals, up_err, low_err)
    return (format_to_decimals(val, decimals),
            format_to_decimals(up_err, decimals),
            format_to_decimals(low_err, decimals))

def format_value_and_error(val, err, decimals=None):
    '''
    Writes a value and its symmetric uncertainty to the same number of decimal places.

    The symmetric case of robust_decimal_errors, for a catalogue value quoted with a single error.

    Parameters
    -----------
    val: the value to round
    err: the error on the value
    decimals: the number of decimal places to write both with. Taken from the precision of the value
        when it is None.
    '''
    val_str, err_str, _ = robust_decimal_errors(val, err, err, decimals=decimals)
    return val_str, err_str

def grab_medians(path, file_prefix, bimodal=False):
    '''
    Collects median values from EXOFASTv2 output files at the defined path.

    The median value and the two uncertainties are returned as strings holding exactly the digits
    EXOFASTv2 wrote, with any scientific-notation exponent already folded in. Reading them as numbers
    would throw away the precision the file states: a median of '0.04530' would come back as 0.0453
    and the table would lose the last significant figure EXOFASTv2 chose to report. Callers that want
    to compute with a value rather than print it should call float() on it themselves.

    Parameters
    -----------
    path: path to the fit files. Ex: '/Users/jack/Research/pipelines/system_figure_pipeline/data/'
    file_prefix: prefix used in EXOFASTv2 output files. Ex: 'toi1855'
    bimodal: boolean to decide whether or not to use the files that are output by splitpdf.pro
    '''

    median_names= ['parname', 'median_value', 'upper_error', 'lower_error', 'scinot']
    suffix = '.csv' if bimodal else '.median.csv'
    medians = pd.read_csv(path + file_prefix + suffix, names=median_names, header=None, skiprows=1,
                          dtype=str)
    for column in median_names:
        medians[column] = medians[column].str.strip()

    medians_corrected = medians.copy()
    for i in medians.index:
        exponent = _scinot_exponent(medians.scinot[i])
        for column in ('median_value', 'upper_error', 'lower_error'):
            medians_corrected.loc[i, column] = _apply_exponent(medians[column][i], exponent)
    return medians_corrected

def _scinot_exponent(scinot):
    '''
    Reads the exponent out of the scientific-notation column of an EXOFASTv2 median file.

    The column holds a LaTeX fragment such as '\\times 10^{-10}' for the parameters EXOFASTv2 factors a
    power of ten out of, and nothing at all for the rest.

    Parameters
    -----------
    scinot: the contents of the column, which is NaN where the parameter has no exponent
    '''
    if not isinstance(scinot, str):
        return 0
    exp_search = re.findall(r'\\times 10\^{(.*)}', scinot)
    return int(exp_search[0]) if exp_search else 0

def _apply_exponent(value, exponent):
    '''
    Multiplies a number written as a string by a power of ten, keeping every digit it was written with.

    Shifting the decimal point rather than multiplying by a float keeps the significant figures exact,
    so a bolometric flux of '5.40' with an exponent of -10 becomes '0.000000000540' and keeps the
    trailing zero that states its precision.

    Parameters
    -----------
    value: the number to scale, as a string
    exponent: the power of ten to apply
    '''
    try:
        number = to_decimal(value)
    except (InvalidOperation, ValueError, TypeError):
        return value # not a number, so leave it for the caller to notice
    if exponent == 0:
        return format(number, 'f')
    with localcontext() as context:
        context.prec = _DECIMAL_PRECISION
        return format(number.scaleb(exponent), 'f')

def median_scinot_corrections(median, parname):
    '''
    Applies the scientific-notation exponent of a parameter in the EXOFASTv2 median table to its value
    and uncertainties, returning all three as strings with their digits intact.

    Parameters
    -----------
    median: pandas DataFrame for median table
    parname: input parameter name
    '''

    row = median[median.parname == parname].iloc[0]
    exponent = _scinot_exponent(row.scinot)
    return tuple(_apply_exponent(row[column], exponent)
                 for column in ('median_value', 'upper_error', 'lower_error'))


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
    # near-infrared narrow bands, as high-resolution imagers use them
    'Jcont': r'$Jcont$', 'Hcont': r'$Hcont$', 'Kcont': r'$Kcont$',
    'Brgamma': r'$Br\gamma$', 'Br-gamma': r'$Br\gamma$',
}

# units a filter's wavelength can be given in, written upright
FILTER_UNITS = {'nm': 'nm', 'microns': r'\textmu m', 'micron': r'\textmu m', 'um': r'\textmu m'}

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
    does, keeps only the filter. One named only by its central wavelength, as a speckle
    camera's "562: 562 (54) nm" is, becomes "562 nm", with neither the number nor the unit
    italicised.

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

        # a filter named by its central wavelength keeps the wavelength and the unit it is in
        wavelength = re.fullmatch(r'\d+(?:\.\d+)?', name.split(':')[0].strip())
        unit = re.search(r'\)\s*([A-Za-z]+)\s*$', name)
        if wavelength and unit and unit.group(1).lower() in FILTER_UNITS:
            return f'{wavelength.group()} {FILTER_UNITS[unit.group(1).lower()]}'

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
