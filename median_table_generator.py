"""Generates the tables of median EXOFASTv2 fit parameters."""

import logging
import numpy as np
import os
import pandas as pd
import re
from table_utils import (grab_medians, is_finite_number, remove_sci_notation, robust_decimal_errors,
                         round_sig_figs, write)

def format_value_with_errors(val, up_err, low_err):
    '''
    Format a value and its asymmetric uncertainties so:
    - scientific notation is removed
    - the uncertainties are written to the same number of decimal places as the value

    The value governs the precision: both uncertainties are rounded to the number of decimal places the
    value is written with, so that the three numbers line up and none of them loses a trailing zero. The
    three arrive as the strings grab_medians read out of the EXOFASTv2 median file, which carry exactly
    the digits EXOFASTv2 chose to report; turning them into floats first would drop every trailing zero
    and with it the precision the fit reported.

    Parameters
    -----------
    val: the value to format
    up_err: the upper uncertainty on the value
    low_err: the lower uncertainty on the value

    Returns
    -------
    tuple[str, str | None, str | None]
        The value and its two uncertainties as strings. A value that is missing or not a number comes
        back as the '---' filler with no uncertainties, and a value whose uncertainties are missing
        comes back on its own, both signalled by a None in place of each uncertainty.
    '''
    if not is_finite_number(val):
        return '---', None, None

    if not (is_finite_number(up_err) and is_finite_number(low_err)):
        # No reliable errors: just return the value without sci notation
        return remove_sci_notation(val), None, None

    return robust_decimal_errors(val, up_err, low_err)

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
            match = priors.loc[priors.variable == candidate]
            if len(match):
                # a variable can appear several times: once under '# intializing linked parameters'
                # holding only a starting value, and again where its prior is actually declared.
                # Link to the row carrying the prior when there is one.
                declared = match[_defines_gaussian_prior(match)]
                source = declared if len(declared) else match
                priors.loc[i, 'meanvalue'] = source.meanvalue.iloc[0]
                break
    # to ensure that all mean values are floats. A link whose target is absent (e.g. the target row is
    # commented out in the prior file) is left as NaN rather than raising and losing the whole table.
    priors['meanvalue'] = pd.to_numeric(priors['meanvalue'], errors='coerce')
    return priors

def _defines_gaussian_prior(priortable):
    '''
    Flags the rows of a prior table that actually declare a Gaussian prior, as opposed to the rows that
    merely repeat a variable.

    EXOFASTv2 writes a variable's name more than once in a '.priors.final' file. Under the
    '# intializing linked parameters' header it lists starting values alone (no standard deviation), and a
    parameter linked to another is written as 'dilute_1 dilute_0 0' with a zeroed standard deviation. Only a
    row with a finite, positive standard deviation states a Gaussian prior, so this is what distinguishes
    the prior from a starting value that happens to be listed first.

    Parameters
    -----------
    priortable: Pandas DataFrame of priors obtained using the grab_priors function

    Returns a boolean Series aligned with priortable.
    '''
    stdev = pd.to_numeric(priortable['stdev'], errors='coerce')
    return np.isfinite(stdev) & (stdev > 0)

def _defines_uniform_prior(priortable):
    '''
    Flags the rows of a prior table that declare a uniform prior, which EXOFASTv2 writes with a standard
    deviation of -1 followed by the lower and upper bounds. As with the Gaussian case this separates the
    declared prior from rows that only repeat the variable's starting value.

    Parameters
    -----------
    priortable: Pandas DataFrame of priors obtained using the grab_priors function

    Returns a boolean Series aligned with priortable.
    '''
    # the bounds are tested for presence rather than finiteness so that a prior bounded at infinity is
    # still recognized as declared
    return priortable['low_bound'].notna() & priortable['up_bound'].notna()

def make_median_string(medians, param, array, star_index=0):
    '''
    Turns median values into strings containing LaTeX, ready to be entered into the median table.

    Parameters
    -----------
    medians: Pandas DataFrame containing the median values obtained using the grab_medians function.
    param: the parameter to generate a string for
    array: the array corresponding to the table row that the parameter should be appended to
    star_index: the index of the star the parameter belongs to, matching the suffix EXOFASTv2 writes in the
        median file. 0 is the target star (teff_0), 1 is the first secondary star (teff_1), and so on.

    Exactly one cell string is appended per call so that the row can be sliced by target index
    (e.g. when splitting a long target list across multiple tables).
    '''
    param = f'{param}_{star_index}'

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

        if up_str == low_str: # compared as written, so that errors that round alike share one \pm
            errstring = r' \pm ' + up_str
        else:
            errstring = r'^{+' + up_str + '}_{-' + low_str + '}'

        array.append('& $' + val_str + errstring + '$ ')
    else:
        array.append('& ---')


def med_table(target_list, path, file_prefix_list, outputpath='.', bimodal=False, parameters=None,
              probabilities=None, max_targets_per_table=5, secondary_stars=False, host_list=None,
              star_types=None, MNRAS=False):
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
    MNRAS: set to True to use the MNRAS class's \\contcaption for the continuation tables, in place of the
        "\\caption{\\textit{(Continued)}}" and "\\addtocounter{table}{-1}" pair written otherwise. \\contcaption
        supplies the continuation wording and holds the table number itself, so neither is needed alongside it.
    secondary_stars: set to True to build a table of the median stellar parameters of the secondary stars instead of the target
        stars. The values are read from the same median files, off the parameters EXOFASTv2 writes with a star index above 0
        (teff_1, feh_1, ...), so target_list holds the name of each secondary star while file_prefix_list holds the fit of its
        host, as in secondary_stars_table(). Only stellar parameters are written, since the planetary parameters of the fit
        belong to the target star, and the priors block is replaced by the classification and planet host rows described under
        star_types and host_list. Hosts whose fit has no secondary star are dropped from the table, and their names are
        printed. The output is named secondary_median_table.tex. A host with more than one secondary star is listed once per
        star, and the repeats take the star indices in the order they appear, so the first is star 1 and the second is star 2.
        Passing an integer instead of True reads that star index for every column (secondary_stars=2 reads teff_2, feh_2, ...).
        The symbols of the initial metallicity, age, V-band extinction and distance rows are marked with a superscript
        asterisk, and a note explaining that a bound companion takes those four from the primary star is written beneath the
        last table. The note is left out if none of the four are in the table.
    host_list: optional array of the names of the hosts, one per entry of target_list. If given, they are written in a 'Planet
        Host' row beneath the star names. Only used when secondary_stars is set.
    star_types: optional dict of classifications describing what each secondary star is, such as 'Bound companion' or
        'Background star'. Keys are the names in target_list. The classification row is only written if at least one secondary
        star has been classified. Only used when secondary_stars is set.
    '''

    # `secondary_stars` picks which star of the fit the table is built from: False keeps the target star,
    # True gives each column the next secondary star of its host, and an integer fixes the star index
    secondary = bool(secondary_stars)
    fixed_star_index = int(secondary_stars) if (secondary and secondary_stars is not True) else None
    if fixed_star_index is not None and fixed_star_index < 1:
        raise ValueError(f'secondary_stars must be True, False, or a star index of 1 or more, not {secondary_stars!r}.')

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
    if host_list is not None and len(host_list) != len(file_prefix_list):
        raise ValueError(f'host_list has {len(host_list)} entries but file_prefix_list has {len(file_prefix_list)}; '
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

    # a bound companion takes these four from the primary star rather than being fit for them, so in a
    # secondary-star table their symbols carry a superscript asterisk that the note beneath the table
    # explains. The marked symbol is written out for each one, since where the asterisk goes depends on the
    # sub- and superscripts the symbol already carries.
    fixed_to_primary = ('initfeh', 'age', 'Av', 'distance')
    if secondary:
        for labels, marked_symbol in ((initfehs, r'$[{\rm Fe/H}]_{0}^{*}$'),
                                      (ages, r'Age$^{*}$'),
                                      (avs, r'$A_V^{*}$'),
                                      (dists, r'$d^{*}$')):
            _, separator, description = labels[0].partition(' & ')
            labels[0] = marked_symbol + separator + description
    show_secondary_note = secondary and any(_normalize_param_name(name) in selected_parameters
                                            for name in fixed_to_primary)

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

    # the planetary parameters of a fit belong to the target star, so a secondary-star table holds
    # nothing but the stellar parameters
    parameter_specs = stellar_parameter_specs if secondary else stellar_parameter_specs + planetary_parameter_specs

    # Setting up to save the table as a .tex file

    if os.path.exists(outputpath) == False:
        os.mkdir(outputpath)

    # Reading the median files up front, so that a secondary-star table can work out which star each column
    # holds, and drop the hosts that have no secondary star, before the columns are laid out
    median_tables = [grab_medians(path=path, file_prefix=prefix, bimodal=bimodal) for prefix in file_prefix_list]
    star_indices = [0] * len(file_prefix_list) # the star of its fit that each column is read from

    if secondary:
        # a star suffix on its own does not mark a secondary star, since EXOFASTv2 also indexes the transit,
        # telescope and limb darkening parameters (dilute_1, gamma_1 and u1_1 all turn up in a fit of a single
        # star), so only the suffixes carried by one of the stellar parameters are counted
        stellar_keys = {median_key for _, median_key, _ in stellar_parameter_specs}

        def _secondary_indices(medians):
            '''Returns the sorted star indices of the secondary stars held in a median table.'''
            found = set()
            for parname in medians.parname.astype(str):
                key, _, index = parname.rpartition('_')
                if key in stellar_keys and index.isdigit() and int(index) > 0:
                    found.add(int(index))
            return sorted(found)

        kept = []
        dropped = []
        star_indices = []
        for ii, prefix in enumerate(file_prefix_list):
            available = _secondary_indices(median_tables[ii])
            if fixed_star_index is not None:
                index = fixed_star_index if fixed_star_index in available else None
            else:
                # a host with more than one secondary star is listed once per star, and the repeats take the
                # star indices in the order they appear
                occurrence = list(file_prefix_list[:ii]).count(prefix)
                index = available[occurrence] if occurrence < len(available) else None
            if index is None:
                dropped.append(str(target_list[ii]))
            else:
                kept.append(ii)
                star_indices.append(index)

        missing = f'star {fixed_star_index}' if fixed_star_index is not None else 'secondary star'
        if dropped:
            print(f'No {missing} was found in the fit of ' + ', '.join(dropped) +
                  '; leaving them out of the table...')
        if not kept:
            raise ValueError(f'None of the given fits have a {missing}, so there is no table to write. Median '
                             'files hold the secondary stars as teff_1, feh_1, and so on.')
        target_list = [target_list[ii] for ii in kept]
        file_prefix_list = [file_prefix_list[ii] for ii in kept]
        median_tables = [median_tables[ii] for ii in kept]
        if host_list is not None:
            host_list = [host_list[ii] for ii in kept]
        if probability_cells is not None:
            probability_cells = [probability_cells[ii] for ii in kept]

    # The rows that identify each secondary star, as in secondary_stars_table(): what the star is, and which
    # planet host it belongs to. Neither row is written unless at least one column has something to put in it.
    host_labels = []
    type_labels = []
    if secondary:
        for ii, name in enumerate(target_list):
            host_labels.append(str(host_list[ii]) if host_list is not None else '---')
            star_type = star_types.get(str(name)) if star_types else None
            type_labels.append(str(star_type) if star_type is not None else '---')
    hosts_named = any(host != '---' for host in host_labels)
    classified = any(star_type != '---' for star_type in type_labels)

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
    # ('bimodal_median_table.tex' for a bimodal table, 'secondary_median_table.tex' for a secondary-star
    # table), bumping a numeric suffix if it already exists; continuation tables append '_2', '_3', ...
    base_name = ('bimodal_' if bimodal else '') + ('secondary_' if secondary else '') + 'median_table'
    first_name = f'{base_name}.tex'
    suffix = 2
    while os.path.exists(f'{outputpath}/{first_name}'):
        first_name = f'{base_name}_{suffix}.tex'
        suffix += 1
    stem = first_name[:-len('.tex')]
    filenames = [first_name] + [f'{stem}_{k}.tex' for k in range(2, n_chunks + 1)]
    print('Saving this table as ' + ', '.join(filenames) + '...')

    for ii in range(len(target_list)):
        medians = median_tables[ii]

        for param_name, median_key, labels in parameter_specs:
            if _normalize_param_name(param_name) in selected_parameters or _normalize_param_name(median_key) in selected_parameters:
                make_median_string(medians, median_key, labels, star_index=star_indices[ii])

    # Collecting priors to put at the top of the table (one ' & ...' entry per target so they can be sliced)
    parallax_prior = [] # initializing lists
    metallicity_prior = []
    extinction_prior = []
    dilution_prior = []
    # the prior files are only read for a normal (non-split) fit of the target star
    show_priors = (bimodal == False) and (secondary == False)
    if show_priors:
        dilute_bool = np.zeros_like(target_list) # to keep track of which targets were fit for dilution

        def prior_value(priortable, column, variable, kind='gaussian'):
            # EXOFASTv2 prior files are inconsistent about capitalization (e.g. 'Av' vs 'av'),
            # so match the variable name case-insensitively; return NaN if it is absent.
            rows = priortable[priortable.variable.str.lower() == variable.lower()]
            if not len(rows):
                return np.nan
            # a variable is usually listed more than once, so prefer the row that declares the prior over
            # one that only repeats the starting value (see _defines_gaussian_prior)
            declares = _defines_gaussian_prior if kind == 'gaussian' else _defines_uniform_prior
            declared = rows[declares(rows)]
            return (declared if len(declared) else rows)[column].iloc[0]

        def gaussian_prior_cell(mean, stdev, name, prefix):
            '''Writes a Gaussian prior as a table cell, to five significant figures, or the '---' filler.

            A prior whose mean or standard deviation is missing is left blank rather than written as
            'nan', which is what happens when a prior file links a variable to a row that is not there
            (see grab_priors).
            '''
            if not (is_finite_number(mean) and is_finite_number(stdev)):
                logging.warning(f'{prefix} has no usable {name} prior (mean {mean}, standard deviation '
                                f'{stdev}); leaving that cell blank.')
                return r'& --- '
            return r'& $\mathcal{G}$[' + round_sig_figs(mean, 5) + r', ' + round_sig_figs(stdev, 5) + r'] '

        for ii in range(len(target_list)):
            priortable = grab_priors(file_prefix_list[ii], path)
            parallax_prior_mean = prior_value(priortable, 'meanvalue', 'parallax')
            parallax_prior_stdev = prior_value(priortable, 'stdev', 'parallax')
            parallax_prior.append(gaussian_prior_cell(parallax_prior_mean, parallax_prior_stdev,
                                                      'parallax', file_prefix_list[ii]))
            metallicity_prior_mean = prior_value(priortable, 'meanvalue', 'feh')
            metallicity_prior_stdev = prior_value(priortable, 'stdev', 'feh')
            metallicity_prior.append(gaussian_prior_cell(metallicity_prior_mean, metallicity_prior_stdev,
                                                         'metallicity', file_prefix_list[ii]))
            extinction_prior_upperbound = prior_value(priortable, 'up_bound', 'Av', kind='uniform')
            if is_finite_number(extinction_prior_upperbound):
                extinction_prior.append(r'& $\mathcal{U}$[0, ' + round_sig_figs(extinction_prior_upperbound, 5) + r'] ')
            else:
                logging.warning(f'{file_prefix_list[ii]} has no usable V-band extinction prior '
                                f'(upper bound {extinction_prior_upperbound}); leaving that cell blank.')
                extinction_prior.append(r'& --- ')

            # finding the dilution term. A dilution variable turns up several times in a prior file: under
            # '# intializing linked parameters' with only a starting value, once where its prior is
            # declared, and once per linked dilute_N row with a zeroed stdev. Take the row that declares a
            # Gaussian prior rather than whichever row comes first, which would otherwise report the
            # starting value with a NaN uncertainty.
            dilute_rows = priortable[priortable.variable.astype(str).str.contains('dilute')]
            declared_dilution = dilute_rows[_defines_gaussian_prior(dilute_rows)]
            if len(declared_dilution):
                dilution_prior_mean = declared_dilution.meanvalue.iloc[0]
                dilution_prior_stdev = declared_dilution.stdev.iloc[0]
                dilute_bool[ii] = 1
            elif len(dilute_rows):
                # dilution was a fitted parameter but no Gaussian prior was placed on it, so the priors
                # block has nothing to report for this target
                logging.warning(f'{file_prefix_list[ii]} fits dilution but declares no Gaussian dilution prior; '
                                'leaving the dilution prior blank.')
            if dilute_bool[ii]:
                dilution_prior.append(gaussian_prior_cell(dilution_prior_mean, dilution_prior_stdev,
                                                          'dilution', file_prefix_list[ii]))
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
        colstring = 'c' * n_chunk
        if bimodal:
            # one system name per group, spanning that system's solution columns
            namestring = ''.join(r' & \multicolumn{' + str(len(group)) + r'}{c}{' + _system_name(group[0]) + r'}'
                                 for group in chunk_groups)
        else:
            namestring = ''.join(' & ' + str(target_list[i]) for i in idx)
        if bimodal:
            title = r'Median Values and 68\% Confidence Intervals for Solutions which are Bimodal in Mass'
        elif secondary:
            title = r'Median Values and 68\% Confidence Intervals for the Fitted Stellar Parameters of Secondary Stars'
        else:
            title = r'Median Values and 68\% Confidence Intervals for Fitted Stellar and Planetary Parameters'
        if is_first:
            caption = r'\caption{' + title + '}'
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
                label = 'bimodal' if bimodal else ('secondarymedian' if secondary else 'median')
                fout.write(r'\label{tab:' + label + '}' + '\n')
            fout.write(r'\scriptsize' + '\n' +
                       r'\begin{tabular}{ll' + colstring + '}'+'\n'+
                       r'\hline' + '\n' +
                       r'& ' + namestring + r'\\' +'\n')
            if bimodal:
                # which solution each column holds, under the system name it belongs to
                fout.write(r'& ' + ''.join(' & ' + _solution_label(i) for i in idx) + r'\\' + '\n')
            if classified:
                # what each secondary star is, written first so that the planet host is not read as one
                fout.write(r'& Classification' + ''.join(' & ' + type_labels[i] for i in idx) + r'\\' + '\n')
            if hosts_named:
                fout.write(r'& Planet Host' + ''.join(' & ' + host_labels[i] for i in idx) + r'\\' + '\n')
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

            if not secondary:
                # the planetary parameters of the fit belong to the target star, not to a secondary star
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
            if is_last and show_secondary_note:
                fout.write(r'\begin{flushleft}' + '\n' +
                           r'\textbf{Note:} *When the secondary star is a bound companion, the initial metallicity, age, \textit{V}-band extinction, and distance are fixed to those of the primary star.' + '\n' +
                           r'\end{flushleft}' + '\n')
            fout.write(r'\end{table*}')
            if not is_last and not MNRAS:
                # keep every piece except the last on the same table number, which \contcaption does itself
                fout.write('\n' + r'\addtocounter{table}{-1}')

    for k, fname in enumerate(filenames):
        _write_chunk(fname, chunks[k], is_first=(k == 0), is_last=(k == n_chunks - 1))
