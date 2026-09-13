"""Generates the table of high-resolution imaging observations."""

import pandas as pd
from table_utils import _extract_grid_rows, convert_table_to_latex_and_save, fetch_exofop_page, format_filter_name

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
    html = fetch_exofop_page(tic_id)

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

    # Convert the filter names to the notation used in the literature, escaping the LaTeX
    # special characters in any name that has no known conversion.
    df_short['Filter'] = df_short['Filter'].apply(
        lambda x: format_filter_name(x) if isinstance(x, str) else x)

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
        df['TOI'] = toi_id
        df.loc[1:, 'TOI'] = ''

        master_df = pd.concat([master_df, df], ignore_index=True)
        master_df = master_df[['TIC ID', 'TOI', 'Telescope', 'Date', 'Instrument', 'Imaging Type', 'Filter', 
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
    convert_table_to_latex_and_save(master_df, output_filename,
                                    caption='Summary of High-resolution Imaging Observations',
                                    label='tab:hri')
