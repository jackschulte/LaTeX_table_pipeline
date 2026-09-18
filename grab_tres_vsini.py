"""Collects vsini measurements from the TRES/CHIRON reconnaissance spectroscopy site (tess.exoplanets.dk)."""

import numpy as np
import pandas as pd
import re
import logging
import time
import requests
from bs4 import BeautifulSoup
from io import StringIO

TRES_URL = 'http://tess.exoplanets.dk/'


class TresFetchError(Exception):
    """A page of the TRES site could not be fetched, or did not hold what was expected."""


class TresUnreachableError(TresFetchError):
    """The TRES site did not answer a request however many times it was retried."""


class TresSite:
    """A logged-in session on the TRES site, reused for every target.

    The site is an ASP.NET application whose pages are driven by form postbacks, so a target's
    spectra are reached the way a browser reaches them: searching for the target on the main page,
    opening its candidate page, and pressing "Export Textfile", which writes the table of spectra
    to a text file and opens it. Each of those steps is a single plain HTTP request here, with no
    browser behind them.

    The site now and then refuses connections for several minutes at a time, or stops sending a
    page partway through. Every request is therefore retried, up to `attempts` times with a wait
    that grows by `delay` seconds after each failure (seven minutes of waiting in all, by default),
    and a response that has not arrived in full within `timeout` seconds counts as failed. Each
    retry is logged as a warning.

    Parameters
    ----------
    username : str
        Username for the TRES site.
    password : str
        Password for the TRES site.
    attempts : int
        How many times to make each request before giving up.
    delay : float
        Seconds added to the wait before each successive retry.
    timeout : float
        Seconds a response may take to arrive in full.

    Raises
    ------
    TresUnreachableError
        If the site cannot be reached.
    TresFetchError
        If the login is rejected.
    """

    def __init__(self, username, password, attempts=8, delay=15, timeout=60):
        self.session = requests.Session()
        self.attempts = attempts
        self.delay = delay
        self.timeout = timeout

        url, html = self._fetch('GET', 'Login.aspx')
        form = _form_fields(html)
        form.update({'ctl00$cpMainContent$tbUserName': username,
                     'ctl00$cpMainContent$tbPassword': password,
                     'ctl00$cpMainContent$btnLogin': 'Login'})
        url, html = self._fetch('POST', 'Login.aspx', data=form)
        if url.endswith('Login.aspx'):
            raise TresFetchError('The TRES site rejected the login. Check the username and password.')
        self._search_form = _form_fields(html) # the main page, which the first search is posted from

    def _fetch(self, method, url, data=None):
        """Make a request, retrying it when it fails, and return the final URL and the page's text."""
        url = requests.compat.urljoin(TRES_URL, url)
        problem = None
        for attempt in range(1, self.attempts + 1):
            if attempt > 1:
                logging.warning(f'Retrying {url} on the TRES site ({problem}); attempt {attempt} of {self.attempts}.')
                time.sleep(self.delay * (attempt - 1))
            try:
                with self.session.request(method, url, data=data, timeout=(10, 30), stream=True) as response:
                    deadline = time.monotonic() + self.timeout
                    body = b''
                    for chunk in response.iter_content(65536):
                        body += chunk
                        if time.monotonic() > deadline:
                            raise TimeoutError(f'the response took longer than {self.timeout} s')
                    if response.status_code >= 500:
                        problem = f'the server answered with status {response.status_code}'
                        continue
                    if response.status_code >= 400: # a missing page or refused request will not come right on a retry
                        raise TresFetchError(f'{url} answered with status {response.status_code}.')
                    return response.url, body.decode(response.encoding or 'utf-8', 'replace')
            except (requests.RequestException, TimeoutError) as e:
                # urllib3 wraps the reason in a long chain of connection-pool messages
                reason = re.search(r'\[Errno \d+\] [^\'")]+|Read timed out|Connection to \S+ timed out', str(e))
                problem = reason.group(0) if reason else str(e)

        raise TresUnreachableError(f'Could not fetch {url} from the TRES site in {self.attempts} attempts; '
                                   f'last, {problem}.')

    def spectra_table(self, TICID):
        """Return the table of a target's spectra, as given by the "Export Textfile" button of its candidate page.

        Parameters
        ----------
        TICID : str or int
            TESS Input Catalog identifier of the target.

        Returns
        -------
        pandas.DataFrame
            One row per spectrum, with columns including "vsini", "vsini_err", "spec." and "method".

        Raises
        ------
        TresUnreachableError
            If the site cannot be reached.
        TresFetchError
            If the target is not on the site, or its spectra cannot be exported.
        """
        T0name = 'T0' + str(TICID).zfill(9)

        search = dict(self._search_form, **{'__EVENTTARGET': 'ctl00$cpMainContent$lnkbSearch', '__EVENTARGUMENT': '',
                                            'ctl00$cpMainContent$tbSearch': T0name})
        _, html = self._fetch('POST', 'Default.aspx', data=search)
        self._search_form = _form_fields(html) # the results page is the main page again, ready for the next search
        link = BeautifulSoup(html, 'lxml').find('a', href=re.compile('Candidate_Edit'), string=T0name)
        if link is None:
            raise TresFetchError(f'{T0name} is not on the TRES site.')

        candidate_url, html = self._fetch('GET', link['href'])
        export = dict(_form_fields(html), **{'__EVENTTARGET': '', '__EVENTARGUMENT': '',
                                             'ctl00$cpMainContent$Button12': 'Export Textfile'})
        _, html = self._fetch('POST', candidate_url, data=export)
        textfile = re.search(r"window\.open\('([^']+\.txt)'", html)
        if textfile is None:
            raise TresFetchError(f'Exporting the spectra of {T0name} did not give a text file.')

        _, table = self._fetch('GET', textfile.group(1))
        return pd.read_csv(StringIO(table), sep=r'\s+', header=0)


def _form_fields(html):
    """The fields a browser would submit with a page's form, leaving out its buttons."""
    soup = BeautifulSoup(html, 'lxml')
    fields = {}
    for field in soup.select('form input[name]'):
        kind = field.get('type', 'text').lower()
        if kind in ('submit', 'button', 'image', 'reset', 'file'):
            continue
        if kind in ('checkbox', 'radio') and not field.has_attr('checked'):
            continue
        fields[field['name']] = field.get('value', 'on' if kind in ('checkbox', 'radio') else '')
    for field in soup.select('form select[name]'):
        option = field.find('option', selected=True) or field.find('option')
        if option is not None:
            fields[field['name']] = option.get('value', option.get_text())
    for field in soup.select('form textarea[name]'):
        fields[field['name']] = field.get_text()
    return fields


def grab_tres_vsini(username, password, TICIDs):
    '''
    Grabs the vsini measurements of each target and calculates their mean and standard error of the mean.

    Only SPC 2.9 measurements from TRES spectra and measurements from CHIRON spectra are used. The site is
    logged into once for all of the targets. A target whose spectra cannot be collected, such as one not on
    the site, is logged as a warning and given NaN for both values, so that it does not stop the others from
    being collected. If the site itself stops answering, though, the collection is stopped, since every
    remaining target would otherwise wait through its own round of retries only to fail the same way.

    Parameters
    -----------
    username: String including the user's TRES website username
    password: String including the user's TRES website password
    TICIDs: an array of TIC IDs, as strings or integers

    Returns
    -----------
    A list holding a (mean vsini, standard error of the mean) pair for each target, in the order of TICIDs

    Raises
    -----------
    TresUnreachableError: if the site cannot be reached
    TresFetchError: if the login is rejected
    '''

    site = TresSite(username, password)

    vsinis = []
    for TICID in TICIDs:
        try:
            data = site.spectra_table(TICID)
        except TresUnreachableError:
            raise
        except TresFetchError as e:
            logging.warning(f'Could not collect the TRES vsini of TIC {TICID}: {e}')
            vsinis.append((np.nan, np.nan))
            continue

        data = data[(data.method == 'SPC2.9') | (data['spec.'] == 'chiron')]
        vsini = data.vsini

        mean_vsini = np.mean(vsini)
        SEM_vsini = np.std(vsini, ddof=1) / np.sqrt(np.size(vsini)) # standard error of the mean
        vsinis.append((mean_vsini, SEM_vsini))

    return vsinis
