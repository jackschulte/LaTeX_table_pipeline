"""Generates the table of observed properties of the secondary stars."""

import numpy as np
import os
import pandas as pd
from astroquery.vizier import Vizier
from astropy.coordinates import Angle
import logging
from table_utils import _extract_grid_rows, fetch_exofop_page, grab_medians, remove_sci_notation, round_sig_figs, write

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

    Raises
    ------
    ExoFOPFetchError
        If the page cannot be fetched in full, so that a failed fetch is not mistaken for a
        target without companions.
    """
    html = fetch_exofop_page(tic_id)

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
                          sep_tol=0.1, pa_tol=5.0, dmag_tol=0.05, max_stars_per_table=5, MNRAS=False):
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
    MNRAS: set to True to use the MNRAS class's \\contcaption for the continuation tables, in place of the
        "\\caption{\\textit{(Continued)}}" and "\\addtocounter{table}{-1}" pair written otherwise. \\contcaption
        supplies the continuation wording and holds the table number itself, so neither is needed alongside it.
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
    gaia_g_arr = [r'$G$ & Gaia $G$ mag. ']
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
            logging.warning(f'There are no secondary stars in the SED file of {prefix}. Skipping {label}.')
            continue
        if occurrence >= len(secondary_indices):
            logging.warning(f'More secondary stars are listed for {prefix} than its SED file holds. Skipping {label}.')
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
                logging.warning(f'There is no distance in the median table of {prefix}. The projected separation of {label} will be omitted.')
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
                logging.warning(f'Could not fetch the ExoFOP companions of TIC {tic}: {e}')
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
                logging.warning(f'Could not match {label} to an ExoFOP companion. Its separation will be omitted.')
            elif ambiguous:
                logging.warning(f'{label} matches more than one ExoFOP companion equally well. Its separation will be omitted.')
            else:
                sep, sep_err = _companion_separation(matches)
        else:
            logging.warning(f'ExoFOP lists no companions for TIC {tic}. The separation of {label} will be omitted.')

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
        logging.warning('No secondary stars were found for any target. No table was generated.')
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

        if is_first:
            caption = r'\caption{Observed Properties of Secondary Stars}'
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
                fout.write(r'\label{tab:secondarylit}' + '\n')
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
            if not is_last and not MNRAS:
                # keep every piece except the last on the same table number, which \contcaption does itself
                fout.write('\n' + r'\addtocounter{table}{-1}')

    for k, fname in enumerate(filenames):
        sl = slice(k * chunk_size, min((k + 1) * chunk_size, n_stars))
        _write_chunk(fname, sl, is_first=(k == 0), is_last=(k == n_chunks - 1))
