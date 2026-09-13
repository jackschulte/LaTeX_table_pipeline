"""Generates the table of ground-based photometric follow-up observations."""

import os
import pandas as pd
import re
import glob
import subprocess
import logging
from table_utils import (FOLLOWUP_FILTER_NAMES, _FILTER_SEPARATOR, _extract_grid_rows,
                         convert_table_to_latex_and_save, fetch_exofop_page, format_filter_name)

# EXOFASTv2 names its transit files "n<YYYYMMDD>.<filter>.<telescope>.<target>.dat",
# e.g. "n20211119.Sloani.KeplerCam.TOI-3788.dat"
LIGHTCURVE_DATE = re.compile(r'^n(\d{8})$')

# An aperture written with a decimal point, as in "GdP-0.4m", is split in two by the full
# stops that separate the fields; this matches the half that got left behind.
LIGHTCURVE_APERTURE = re.compile(r'^\d+m$')


def parse_lightcurve_filename(filename):
    """Read the date, filter, and telescope out of an EXOFASTv2 lightcurve filename.

    The fields are separated by full stops, but a telescope name may contain one itself when
    its aperture is part of the name ("GdP-0.4m"), so a field that is only an aperture is
    joined back onto the telescope before it. Whatever follows the telescope is the target
    and, for the TESS lightcurves, the sector and cadence they were binned to; none of it
    is read.

    Parameters
    ----------
    filename : str
        A filename such as "n20211119.Sloani.KeplerCam.TOI-3788.dat".

    Returns
    -------
    dict or None
        The Date, Filter, and Telescope, or None if the name is not a lightcurve.
    """
    fields = os.path.splitext(os.path.basename(filename))[0].split('.')
    if len(fields) < 3:
        return None

    date = LIGHTCURVE_DATE.match(fields[0])
    if date is None:
        return None

    telescope = fields[2]
    for field in fields[3:]:
        if not (telescope[-1:].isdigit() and LIGHTCURVE_APERTURE.match(field)):
            break
        telescope += '.' + field

    return {'Date': pd.to_datetime(date.group(1), format='%Y%m%d'),
            'Filter': fields[1],
            'Telescope': telescope}


def read_lightcurve_files(target_folder_names, file_pattern='n2*.dat',
                          hpcc_path='jschulte@rsync.hpcc.msu.edu:/mnt/research/Exoplanet_Lab/jack/Global_Fits/'):
    """Compile the date, filter, and telescope of every transit lightcurve fit for each object.

    The files themselves are never copied over; only their names are read, which is all that
    is needed to know which observations went into a fit. When hpcc_path names a remote
    machine the listing is done with a single ssh call, so the connection is only made once
    however many objects are asked for.

    Parameters
    ----------
    target_folder_names : str or list[str]
        The folder holding each object's fit files, relative to hpcc_path, e.g.
        "meep3/toi3788". A single string is accepted for one object.
    file_pattern : str
        The glob matching the lightcurve files inside each folder. The default matches the
        "n<YYYYMMDD>" every lightcurve name opens with, without picking up the other files
        a fit folder holds.
    hpcc_path : str
        The path the folders sit under. A "user@host:" prefix reads the names over ssh; a
        path without one is read from the local filesystem instead.

    Returns
    -------
    dict[str, pandas.DataFrame]
        One table per target folder, with the Date, Filter, Telescope, and Filename of each
        lightcurve. A folder with no lightcurve files gets an empty table.
    """
    if isinstance(target_folder_names, str):
        target_folder_names = [target_folder_names]

    if ':' in hpcc_path:
        host, base_path = hpcc_path.split(':', 1)
    else:
        host, base_path = None, hpcc_path
    base_path = base_path.rstrip('/')

    patterns = [f"{base_path}/{folder.strip('/')}/{file_pattern}" for folder in target_folder_names]

    if host is None:
        filepaths = [path for pattern in patterns for path in sorted(glob.glob(pattern))]
    else:
        # ls expands the globs itself, and prints each path as it was given, so one call
        # covers every object and the paths still say which folder they came from.
        result = subprocess.run(['ssh', host, 'ls -1 ' + ' '.join(patterns)],
                                capture_output=True, text=True)
        filepaths = [line.strip() for line in result.stdout.splitlines() if line.strip()]
        if not filepaths and result.stderr.strip():
            print(f'No lightcurve files could be read: {result.stderr.strip()}')

    lightcurves = {folder: [] for folder in target_folder_names}
    for filepath in filepaths:
        lightcurve = parse_lightcurve_filename(filepath)
        if lightcurve is None:
            continue

        # Attribute the file to the folder its path passes through.
        for folder in target_folder_names:
            if f"/{folder.strip('/')}/" in filepath:
                lightcurve['Filename'] = os.path.basename(filepath)
                lightcurves[folder].append(lightcurve)
                break

    return {folder: pd.DataFrame(rows, columns=['Date', 'Filter', 'Telescope', 'Filename'])
            for folder, rows in lightcurves.items()}


# A telescope size written into the telescope's name, which the size column already gives:
# "1m0" and "0m4p" as LCO writes them, "0.4m", "0-4m", "036m", or "12in", with or without a
# separator. Only a number with a unit is taken, so model numbers such as "CDK20" or "T100" stay.
TELESCOPE_SIZE_IN_NAME = re.compile(r'[-_\s]?(?:\d+m\d+p?|\d*\.\d+m|\d+-\d+m|\d+m|\d+(?:\.\d+)?\s?(?:cm|in|inch))$')


def _trim_telescope_size(name):
    """Remove a size written into a telescope's name, e.g. "1m0" from "LCO-McD-1m0"."""
    if not isinstance(name, str):
        return name
    trimmed = TELESCOPE_SIZE_IN_NAME.sub('', name).rstrip(' -_')
    return trimmed or name


# the columns of one target's follow-up table, before its TIC and TOI are added
FOLLOWUP_COLUMNS = ['Telescope', 'Tel. Size (m)', 'Date', 'Camera', 'Filter', r'Pix. Scale ($\arcsec$/pix)',
                    r'PSF FWHM ($\arcsec$)', r'Aper. Rad. ($\arcsec$)']


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
    html = fetch_exofop_page(tic_id)

    rows = _extract_grid_rows(html, 'Time Series Observations')

    df = pd.DataFrame(rows)
    df_short = df[['tstel', 'tsdate', 'tscam', 'tsfilt', 'tspix', 'tspsf', 'tspar']].copy()
    df_short.columns = ['Telescope', 'Date', 'Camera', 'Filter', r'Pix. Scale ($\arcsec$/pix)', r'PSF FWHM ($\arcsec$)', r'Aper. Rad. ($\arcsec$)']

    # Extract telescope size (m) and remove it from the Telescope column.
    df_short['Tel. Size (m)'] = df_short['Telescope'].str.extract(r'\((\d*\.?\d*)\s*m\)', expand=False).astype(float)
    df_short['Telescope'] = df_short['Telescope'].str.replace(r'\s*\(\d*\.?\d*\s*m\)', '', regex=True).str.strip()
    # Some names repeat the size inside them too ("LCO-McD-1m0"), which is trimmed as well.
    df_short['Telescope'] = df_short['Telescope'].apply(_trim_telescope_size)

    df_short = df_short.sort_values(by='Date', ascending=True).reset_index(drop=True)
    df_short = df_short[FOLLOWUP_COLUMNS]

    # reformat date to Year Mon Day format
    df_short['Date'] = pd.to_datetime(df_short['Date'], format='%Y-%m-%d')
    df_short['Date'] = df_short['Date'].dt.strftime('%Y %b %d')
    
    # Convert the filter names to the notation used in the literature, escaping the LaTeX
    # special characters in any name that has no known conversion.
    df_short['Filter'] = df_short['Filter'].apply(
        lambda x: format_filter_name(x) if isinstance(x, str) else x)
    
    # Truncate trailing zeros and convert floats to strings for LaTeX formatting
    float_columns = ['Tel. Size (m)', r'Pix. Scale ($\arcsec$/pix)', r'PSF FWHM ($\arcsec$)', r'Aper. Rad. ($\arcsec$)']
    for col in float_columns:
        df_short[col] = df_short[col].apply(lambda x: ('{:.3f}'.format(x)).rstrip('0').rstrip('.') if pd.notnull(x) else x)

    # Replace NaNs with '---'
    df_short = df_short.fillna('---')

    return df_short


# ExoFOP names a few telescopes by their instrument where the fit files name the site the
# instrument sits at, a difference no comparison of the names themselves can bridge. The key
# is the name ExoFOP uses and the value the name the filenames use; add pairings as they
# turn up among the lightcurves reported as having no ExoFOP match.
TELESCOPE_ALIASES = {'CDK20': 'El_Sauce'}


def _canonical_name(name):
    """Reduce a telescope or filter name to the letters and digits it is built from.

    The LaTeX a filter name is written with ("$i'$"), the prefixes EXOFASTv2 spells out
    ("Sloani"), and the punctuation observatories are quoted with ("LCOGT-CTIO") all get
    stripped, so that the names can be compared across the two sources.
    """
    name = re.sub(r"[$\\{}'_]", '', str(name))
    name = re.sub(r'[^a-z0-9]', '', name.lower())
    # EXOFASTv2 writes the Sloan bands as "Sloani", the table as "i"
    if name.startswith('sloan'):
        name = name[len('sloan'):]
    return name


def _name_tokens(name):
    """The words a telescope or filter name is built from, e.g. "lco", "hal", "0m35".

    Words of one or two characters are left out, being too short to tell anything apart.
    """
    name = re.sub(r"[$\\{}']", '', str(name))
    tokens = (_canonical_name(token) for token in re.split(r'[^A-Za-z0-9]+', name))
    return {token for token in tokens if len(token) > 2}


def _names_agree(table_name, file_name):
    """Whether a telescope or filter in the table is the one named in a lightcurve filename.

    The two sources name the same telescope differently: ExoFOP gives the observatory and
    the instrument in full, while a filename abbreviates and often keeps only one of them.
    They are taken to agree when either whole name contains the other ("LCOGT" matches
    "LCO-SAAO-1m0"), or when any word of one matches a word of the other, which is what
    pairs "Unistellar eVscope2" with "eVscope_TO" and "Celestron C11" with "Herges-C11_Exo".
    Names of one or two characters have to match outright instead, so that the Sloan "i" is
    not read as the Cousins "Ic".
    """
    table_full, file_full = _canonical_name(table_name), _canonical_name(file_name)
    if not table_full or not file_full:
        return False
    if len(table_full) <= 2 or len(file_full) <= 2:
        return table_full == file_full
    if table_full in file_full or file_full in table_full:
        return True

    # the separators are gone from the canonical names, so the words are taken from the
    # names as they were written
    return any(table_token in file_token or file_token in table_token
               for table_token in _name_tokens(table_name)
               for file_token in _name_tokens(file_name))


def _telescopes_agree(table_name, file_name):
    """Whether the table's telescope is the one a lightcurve filename names.

    The names are compared as they are written, and then again through TELESCOPE_ALIASES,
    which carries the telescopes the two sources call by unrelated names.
    """
    if _names_agree(table_name, file_name):
        return True
    alias = TELESCOPE_ALIASES.get(str(table_name).strip())
    return alias is not None and _names_agree(alias, file_name)


# ExoFOP and the lightcurve filenames can date the same night a day apart, in either
# direction, so an observation is looked for this many days either side of a file's date.
DATE_TOLERANCE_DAYS = 1


def _band(name):
    """The band a filter name refers to, without its subscript: "$z_s$" and "Sloanz" are both "z"."""
    name = _canonical_name(name)
    return name[0] if len(name) == 2 and name[1] in 'cs' else name


def _filter_components(name):
    """The filters an ExoFOP filter entry names, e.g. "$g'$" and "$i'$" of "$g'$-$i'$".

    Only an entry made up entirely of recognised filters is split; anything else, such as
    "g-narrow,NaD", is a single filter whatever punctuation it contains.
    """
    pieces = _FILTER_SEPARATOR.split(str(name))[::2]
    known = set(FOLLOWUP_FILTER_NAMES.values())
    return pieces if len(pieces) > 1 and all(piece in known for piece in pieces) else [str(name)]


def _filters_agree(table_filter, file_filter):
    """Whether an ExoFOP filter entry includes the band a lightcurve filename names."""
    band = _band(file_filter)
    return bool(band) and any(_band(component) == band for component in _filter_components(table_filter))


def _observed_filter(table_filter, file_filter):
    """The filter to show in the row of one lightcurve from an ExoFOP observation.

    ExoFOP describes a multi-band observation in one entry ("$g'$, $r'$, $i'$, $z_s$"), but
    each band has a lightcurve and so a row of its own, which shows only its own band. The
    band is written as ExoFOP writes it, that being the more specific of the two ("$z_s$"
    where the filename can only say "Sloanz").
    """
    band = _band(file_filter)
    matching = [component for component in _filter_components(table_filter) if _band(component) == band]
    return matching[0] if len(matching) == 1 else table_filter


def build_lightcurve_rows(exofop, lightcurves, toi_id=''):
    """Write one row for each of a target's ground-based lightcurves, from ExoFOP where possible.

    The lightcurve files decide which rows the table has. Every non-TESS file gets exactly
    one, and ExoFOP observations that no file accounts for are left out. Each file is paired
    with an ExoFOP observation dated within DATE_TOLERANCE_DAYS of it whose telescope or
    filter agrees, preferring one on the same date, then one whose telescope agrees, then one
    whose filter does, and last one fewer files have been paired with already, so that two
    identical entries on one night go to two different files. A pairing across dates is named
    in a warning, since the row then gives ExoFOP's date rather than the filename's. A file with
    no such observation still gets a row, giving the date, telescope (less any size in its
    name) and filter as its filename does and "---" for everything else, and is named in a
    warning too.

    Parameters
    ----------
    exofop : pandas.DataFrame
        The target's follow-up table from get_followup_table, with its dates still written
        as "%Y %b %d". It may be empty.
    lightcurves : pandas.DataFrame
        The target's lightcurve files, as returned by read_lightcurve_files.
    toi_id : str
        The TOI the rows belong to, used only to label the warning.

    Returns
    -------
    pandas.DataFrame
        One row per non-TESS lightcurve, in the order the lightcurves were taken.
    """
    # TESS lightcurves are fit alongside the ground-based ones but never have a row here
    ground_based = lightcurves[(lightcurves['Telescope'] != 'TESS') & (lightcurves['Filter'] != 'TESS')]
    ground_based = ground_based.sort_values(['Date', 'Filename'])

    dates = pd.to_datetime(exofop['Date'], format='%Y %b %d', errors='coerce')
    pairings = {index: 0 for index in exofop.index}

    rows, redated, unmatched, near_miss = [], [], [], False
    for _, lightcurve in ground_based.iterrows():
        candidates, nearby = [], []
        for index, date in dates.items():
            if pd.isna(date):
                continue
            offset = abs((date - lightcurve['Date']).days)
            if offset > DATE_TOLERANCE_DAYS:
                continue
            telescope_agrees = _telescopes_agree(exofop.at[index, 'Telescope'], lightcurve['Telescope'])
            filter_agrees = _filters_agree(exofop.at[index, 'Filter'], lightcurve['Filter'])
            if telescope_agrees or filter_agrees:
                candidates.append(((offset, not telescope_agrees, not filter_agrees, pairings[index]), index))
            else:
                nearby.append(f"{exofop.at[index, 'Date']} {exofop.at[index, 'Telescope']} "
                              f"({exofop.at[index, 'Filter']})")

        if candidates:
            rank, index = min(candidates)
            pairings[index] += 1
            row = exofop.loc[index].to_dict()
            row['Filter'] = _observed_filter(row['Filter'], lightcurve['Filter'])
            if rank[0]:
                redated.append(f"{lightcurve['Filename']} vs {row['Telescope']} on {row['Date']}")
        else:
            row = dict.fromkeys(FOLLOWUP_COLUMNS, '---')
            row['Date'] = lightcurve['Date'].strftime('%Y %b %d')
            # trimmed before escaping, as the size can follow an underscore
            row['Telescope'] = _trim_telescope_size(lightcurve['Telescope']).replace('#', r'\#').replace('_', r'\_')
            row['Filter'] = format_filter_name(lightcurve['Filter'])
            unmatched.append(lightcurve['Filename']
                             + (f" (nearby on ExoFOP: {', '.join(nearby)})" if nearby else ''))
            near_miss |= bool(nearby)
        rows.append(row)

    if redated:
        logging.warning(f'TOI-{toi_id}: {len(redated)} of its {len(ground_based)} non-TESS lightcurve '
                        f'files are paired with an ExoFOP observation on a different date, and their '
                        f'rows give the ExoFOP date: {"; ".join(redated)}')
    if unmatched:
        logging.warning(f'TOI-{toi_id}: {len(unmatched)} of its {len(ground_based)} non-TESS lightcurve '
                        f'files have no matching ExoFOP observation, so their rows give only the '
                        f'date, telescope and filter: {"; ".join(unmatched)}.'
                        + (' If a nearby ExoFOP observation is the same one under another telescope '
                           'name, pair the names up in TELESCOPE_ALIASES.' if near_miss else ''))

    return pd.DataFrame(rows, columns=FOLLOWUP_COLUMNS)


def generate_master_followup_table(tic_list, toi_list, fit_observations_only=False,
                                  target_folder_names=None, **lightcurve_kwargs):
    """Build a combined follow-up table for many TIC/TOI pairs.

    Parameters
    ----------
    tic_list : list[str]
        TIC identifiers, with or without the "TIC " prefix.
    toi_list : list[str]
        TOI identifiers, with or without the "TOI-" prefix.
    fit_observations_only : bool
        Whether to list the observations that were actually fit rather than everything on
        ExoFOP. The lightcurve files then decide the rows, one for each non-TESS file, and
        ExoFOP fills in their details; see build_lightcurve_rows. Requires
        target_folder_names.
    target_folder_names : list[str], optional
        The fit folder of each target, in the same order as tic_list, e.g. "meep3/toi3788".
    **lightcurve_kwargs
        Passed on to read_lightcurve_files, e.g. hpcc_path or file_pattern.

    Returns
    -------
    pandas.DataFrame
        A concatenated table with TIC and TOI columns added to each row group.
    """
    lightcurves = {}
    if fit_observations_only:
        if target_folder_names is None:
            raise ValueError('target_folder_names is needed to tell which observations were fit')
        if len(target_folder_names) != len(tic_list):
            raise ValueError('target_folder_names must have one folder per target')
        # one call, so the ssh connection is made once for every target together
        lightcurves = read_lightcurve_files(target_folder_names, **lightcurve_kwargs)

    master_df = pd.DataFrame()

    for target_index, (tic_id, toi_id) in enumerate(zip(tic_list, toi_list)):
        if not tic_id.startswith('TIC '):
            tic_id = 'TIC ' + tic_id

        if toi_id.startswith('TOI-'):
            toi_id = toi_id.replace('TOI-', '')
        try:
            df = get_followup_table(tic_id)
        except Exception as e:
            logging.warning(f"Error occurred while fetching follow-up table for TOI-{toi_id}: {e}")
            if not fit_observations_only:
                continue
            # the lightcurves still decide the rows, which just go without ExoFOP's details
            df = pd.DataFrame(columns=FOLLOWUP_COLUMNS)

        # Build the rows from the lightcurves while they still belong to one target.
        if fit_observations_only:
            target_lightcurves = lightcurves[target_folder_names[target_index]]
            if target_lightcurves.empty:
                logging.warning(f'No lightcurve files found for TOI-{toi_id}; keeping all of its '
                                f'ExoFOP observations instead.')
            else:
                df = build_lightcurve_rows(df, target_lightcurves, toi_id=toi_id)
            if df.empty:
                logging.warning(f'TOI-{toi_id} has no observations to list; it is left out of the table.')
                continue

        # Add TIC and TOI ids to the first row for this target.
        df['TIC ID'] = tic_id.replace('TIC ', '')
        df.loc[1:, 'TIC ID'] = ''
        df['TOI'] = toi_id
        df.loc[1:, 'TOI'] = ''

        master_df = pd.concat([master_df, df], ignore_index=True)
        master_df = master_df[['TIC ID', 'TOI'] + FOLLOWUP_COLUMNS]
    return master_df

def generate_followup_table(tic_list, toi_list, output_filename, print_summary=True,
                            fit_observations_only=False, target_folder_names=None,
                            **lightcurve_kwargs):
    """Generate the follow-up table and save it as a LaTeX file.

    Parameters
    ----------
    tic_list : list[str]
        List of TIC identifiers.
    toi_list : list[str]
        List of TOI identifiers.
    output_filename : str
        The filename for the output LaTeX file.
    print_summary : bool
        Whether to print the number of follow-up observations and the number of unique
        telescopes that went into the table. Enabled by default.
    fit_observations_only : bool
        Whether to list the observations that were actually fit rather than everything on
        ExoFOP. The lightcurve filenames decide the rows: each non-TESS lightcurve gets one,
        filled in from its ExoFOP observation, or with only its date, telescope and filter
        when ExoFOP has none. Requires target_folder_names.
    target_folder_names : list[str], optional
        The fit folder of each target, in the same order as tic_list, e.g. "meep3/toi3788".
    **lightcurve_kwargs
        Passed on to read_lightcurve_files, e.g. hpcc_path or file_pattern.
    """
    master_df = generate_master_followup_table(tic_list, toi_list,
                                               fit_observations_only=fit_observations_only,
                                               target_folder_names=target_folder_names,
                                               **lightcurve_kwargs)

    if print_summary:
        n_observations = len(master_df)
        n_telescopes = master_df['Telescope'].nunique()
        print(f'Number of follow-up observations: {n_observations}')
        print(f'Number of unique telescopes: {n_telescopes}')

    convert_table_to_latex_and_save(master_df, output_filename, fontsize=r'\scriptsize')
