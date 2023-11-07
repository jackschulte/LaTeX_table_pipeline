import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import astropy
import astropy.units as u
import astropy.constants as const
import glob
import pandas as pd
import matplotlib.colors as colors
from astropy.coordinates import SkyCoord
from astroquery.vizier import Vizier
from astropy.coordinates import Angle

def grab_medians(name, path, file_prefix = '.MIST.SED.'):
    '''
    name: name of the system. Ex: 'TOI-1855'
    path: path to the fit files. Ex: '/Users/jack/Research/pipelines/system_figure_pipeline/data/'
    file_prefix: prefix used in EXOFASTv2 output files
    '''

    median_names= ['parname', 'median_value', 'upper_error', 'lower_error']

    medians = pd.read_csv(path + name + file_prefix + 'median.csv', names=median_names, header=None, skiprows=1)
    return medians

def make_median_string(medians, param, array):
    param = param+'_0'

    if medians.parname.isin([param]).any() == True:
    
        uperr = medians.upper_error[medians.parname == param].iloc[0]
        loerr = medians.lower_error[medians.parname == param].iloc[0]

        remove_sci_notation = lambda x: np.format_float_positional(x, trim='-')

        uperr = remove_sci_notation(uperr)
        loerr = remove_sci_notation(loerr)

        if uperr==loerr:
            errstring = r'\pm ' + str(uperr)
        else:
            errstring = r'^{+'+ str(uperr) + '}_{-' + str(loerr) + '}'
        array.append('& $' + str(medians.median_value[medians.parname == param].iloc[0]))
        #if ii<len(paths)-1:
        array.append(errstring + '$ ')
    else:
        array.append('& ---')

def make_and_write_lit_rows(TICIDs, TESS_mags, TESS_mags_err, file):

    # Initializing rows of the table

    ra_arr=[r'$\alpha_{J2000}\ddagger$ & Right Ascension (RA) ']
    dec_arr=[r'$\delta_{J2000}\ddagger$ & Declination (Dec) ']
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
    pmdec_arr=[r'$\mu_{\alpha}$ & Gaia DR3 proper motion in DEC (mas yr$^{-1}$)']
    parallax_arr=[r'$\pi^\dagger$ & Gaia DR3 Parallax (mas) ']
    vbroad_arr=[r'$v\sin{i_\star}$ & Projected rotational velocity (km s$^{-1}$) ']
    # add note that vbroad includes other effects such as macroturbulence, template mismatch, and other instrumental effects


    for i in range(len(TICIDs)):
        # querying for data
        data_gaia = Vizier.query_region('TIC ' + str(TICIDs[i]), radius=Angle(0.001, "deg"), catalog='I/355/gaiadr3')
        data_2MASS = Vizier.query_region('TIC ' + str(TICIDs[i]), radius=Angle(0.001, "deg"), catalog='II/246/out')
        data_WISE = Vizier.query_region('TIC ' + str(TICIDs[i]), radius=Angle(0.001, "deg"), catalog='II/311/wise')
        
        assert len(data_gaia) == 1, f'Multiple Vizier query results for {TICIDs[i]}'
        assert len(data_2MASS) == 1, f'Multiple Vizier query results for {TICIDs[i]}'
        assert len(data_WISE) == 1, f'Multiple Vizier query results for {TICIDs[i]}'

        ra = data_gaia[0]['RA_ICRS'][0]
        ra_angle = Angle(ra, 'deg')
        ra_str = f'{int(ra_angle.hms[0])}:{int(ra_angle.hms[1])}:{round(float(ra_angle.hms[2]), 3)}'
        dec = data_gaia[0]['DE_ICRS'][0]
        dec_angle = Angle(dec, 'deg')
        dec_str = f'{int(dec_angle.dms[0])}:{np.abs(int(dec_angle.dms[1]))}:{np.abs(round(float(dec_angle.dms[2]), 3))}'

        gaia_g = data_gaia[0]['Gmag'][0]
        gaia_g_err = 0.02 # CHANGE THIS
        gaia_bp = data_gaia[0]['BPmag'][0]
        gaia_bp_err = 0.02 # CHANGE THIS
        gaia_rp = data_gaia[0]['RPmag'][0]
        gaia_rp_err = 0.02 # CHANGE THIS

        j_2mass = data_2MASS[0]['Jmag'][0]
        j_2mass_err = data_2MASS[0]['e_Jmag'][0]
        h_2mass = data_2MASS[0]['Hmag'][0]
        h_2mass_err = data_2MASS[0]['e_Hmag'][0]
        k_2mass = data_2MASS[0]['Kmag'][0]
        k_2mass_err = data_2MASS[0]['e_Kmag'][0]

        wise1 = data_WISE[0]['W1mag'][0]
        wise1_err = data_WISE[0]['e_W1mag'][0]
        wise2 = data_WISE[0]['W2mag'][0]
        wise2_err = data_WISE[0]['e_W2mag'][0]
        wise3 = data_WISE[0]['W3mag'][0]
        wise3_err = data_WISE[0]['e_W3mag'][0]
        wise4 = data_WISE[0]['W4mag'][0]
        wise4_err = data_WISE[0]['e_W4mag'][0]

        pmra = data_gaia[0]['pmRA'][0]
        pmra_err = data_gaia[0]['e_pmRA'][0]
        pmdec = data_gaia[0]['pmDE'][0]
        pmdec_err = data_gaia[0]['e_pmDE'][0]
        parallax = data_gaia[0]['Plx'][0]
        parallax_err = data_gaia[0]['e_Plx'][0]
        vbroad = data_gaia[0]['Vbroad'][0]

        gen_value_str(ra_arr, ra_str)
        gen_value_str(dec_arr, dec_str)
        gen_value_str(gaia_g_arr, gaia_g, gaia_g_err)
        gen_value_str(gaia_bp_arr, gaia_bp, gaia_bp_err)
        gen_value_str(gaia_rp_arr, gaia_rp, gaia_rp_err)
        gen_value_str(tmag_arr, TESS_mags[i], TESS_mags_err[i])
        gen_value_str(j_2mass_arr, j_2mass, j_2mass_err)
        gen_value_str(h_2mass_arr, h_2mass, h_2mass_err)
        gen_value_str(k_2mass_arr, k_2mass, k_2mass_err)
        gen_value_str(wise1_arr, wise1, wise1_err)
        gen_value_str(wise2_arr, wise2, wise2_err)
        gen_value_str(wise3_arr, wise3, wise3_err)
        gen_value_str(wise4_arr, wise4, wise4_err)
        gen_value_str(pmra_arr, pmra, pmra_err)
        gen_value_str(pmdec_arr, pmdec, pmdec_err)
        gen_value_str(parallax_arr, parallax, parallax_err)
        gen_value_str(vbroad_arr, vbroad)
        
    
    write(ra_arr, file)
    write(dec_arr, file)
    write(pmra_arr, file)
    write(pmdec_arr, file)
    write(parallax_arr, file)
    write(vbroad_arr, file)
    file.write(r'\multicolumn{' + str(len(TICIDs) + 2) + r'}{l}{\textbf{Photometric Parameters}:} \\' + '\n')
    write(gaia_g_arr, file)
    write(gaia_bp_arr, file)
    write(gaia_rp_arr, file)
    write(tmag_arr, file)
    write(j_2mass_arr, file)
    write(h_2mass_arr, file)
    write(k_2mass_arr, file)
    write(wise1_arr, file)
    write(wise2_arr, file)
    write(wise3_arr, file)
    write(wise4_arr, file)

def gen_value_str(array, value, error=None):
    if (value == None) or (isinstance(value, np.ma.core.MaskedConstant)):
        array.append('& ---')
    elif (error == None) or (isinstance(error, np.ma.core.MaskedConstant)):
        array.append(r'& ' + str(value))
    else:
        array.append(r'& $' + str(value) + r' \pm ' + str(error) + '$ ')

def write(param_arr,file):
    for ii in param_arr:
        file.write(ii)
    file.write(r'\\'+'\n')

def lit_table(target_list, outputpath='.'):
    '''
    target_list: an array of strings containing the names of each target. Nominally, these should be TOI IDs. Ex: ['TOI-1855', 'TOI-2107']
    outputpath: the folder in which the table should be generated. Current working directory by default
    '''

    # Setting up to save the table as a .tex file

    if os.path.exists(outputpath) == False:
        os.mkdir(outputpath)

    newfile = 'lit_table.tex'

    # if this file exists, come up with a new name
    i = 2
    while os.path.exists(f'{outputpath}/{newfile}'):
        newfile = 'lit_table_' + str(i) + '.tex'
        i += 1
    print(f'Saving this table as {newfile}...')

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


    # Generating the preamble

    colstring = 'lc'
    namestring = ''
    
    for ii in range(len(target_list)):
        colstring+='c'
        namestring += (' & \colhead{' + target_list[ii] + '}')
    
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
    r'\tablecolumns{' + str(len(target_list) + 2) + '}'+'\n'+
    r'\tablehead{ & ' + namestring + '}'+'\n'+
    r'\startdata'+'\n'+
    #r'\hline \\' + '\n' + 
    #r'\hline \\' + '\n' + 
    r'\multicolumn{' + str(len(target_list) + 2) + r'}{l}{\textbf{Other identifiers}:} \\' + '\n' +
    r'& \tess Input Catalog & \\' + '\n' +
    r'& TYCHO-2 & \\'  + '\n' +
    r'& 2MASS & \\' + '\n' +
    r'\hline' + '\n' +               
    r'\multicolumn{' + str(len(target_list) + 2) + r'}{l}{\textbf{Astrometric Parameters}:} \\' + '\n' )

        make_and_write_lit_rows(TIC_IDs, TESS_mags, TESS_mags_err, fout)

        fout.write(r'\enddata' + '\n')
    

def med_table(target_list, path, file_prefix = '.MIST.SED.', outputpath='.'):

    '''
    target_list: an array of strings containing the names of each target. Names should match those in EXOFASTv2 output files.
    path: path of EXOFASTv2 output files
    file_prefix: prefix used in EXOFASTv2 file generation
    outputpath: the folder in which the table should be generated. Current working directory by default
    '''

    # set which parameters to output to the .tex file

    # Stellar parameters
    Mstar=True
    Rstar=True
    RstarSED=False
    Lstar=True
    fbol=False
    rhostar=True
    logg=True
    teff=True
    teffsed=False
    feh=True
    initfeh=True
    age=True
    eep=True
    logmstar=False
    Av=True
    errscale=False
    parallax=False # the GAIA parralax is already in the literature table
    distance=True

    # Planetary parameters
    period=True
    rp=True
    mp=True
    mpsun=False #nor sure what the problem is with this
    tc=True
    tt=False
    t0=True
    semimajor=True
    ideg=True
    ecc=True
    omegadeg=True
    teq=True
    tcirc=True
    k=True
    p=True
    ar=True
    delta=False
    depthTESS=True
    tau=True
    t14=True
    tfwhm=False
    b=True
    cosi=False #nor sure what the problem is with this
    bs=True
    taus=False
    ts=False # maybe include this and ts14?
    ts14=False
    tfwhms=False
    depth25=False
    depth50=False
    depth75=False
    rhop=True
    #logp=True
    loggp=True
    safronov=False
    fave=False
    tp=False
    ta=False
    td=False
    vcve=False
    ecosw=False
    esinw=False
    #secosw=True
    #sesinw=True
    msini=False
    q=True
    dr=True
    pt=False
    ptg=False
    ps=False
    psg=False

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

    # Initializing strings with LaTeX for each parameter
    
    mstars=[r'$M_*$ & Mass (\msun) ']
    rstars=[r'$R_*$ & Radius (\rsun) ']
    rstarseds=[r'$R_{*,SED}$ & Radius$ (\rsun) ']
    lstars=[r'$L_*$ & Luminosity (\lsun) ']
    fbols=[r'$F_{Bol}$ & Bolometric flux (cgs) ']
    rhostars=[r'$\rho_*$ & Density (cgs) ']
    loggs=[r'$\log{g}$ & Surface gravity (cgs) ']
    teffs=[r'$T_{\rm eff}$ & Effective temperature (K) ']
    teffseds=[r'$T_{\rm eff,SED}$ & Effective temperature (K) ']
    fehs=[r'$[{\rm Fe/H}]$ & Metallicity (dex) ']
    initfehs=[r'$[{\rm Fe/H}]_{0}$ & Initial metallicity ']
    ages=[r'$Age$ & Age (Gyr) ']
    eeps=[r'$EEP$ & Equal evolutionary phase ']
    logmstars=[r'$\log{M_*}$ & Mass ($\log{\msun}$) ']
    avs=[r'$A_V$ & V-band extinction (mag) ']
    errscales=[r'$\sigma_{SED}$ & SED photometry error scaling ']
    plaxes=[r'$\varpi$ & Parallax (mas) ']
    dists=[r'$d$ & Distance (pc) ']
    
    periods=[r'$P$ & Period (days) ']
    rps=[r'$R_P$ & Radius (\rj) ']
    mps=[r'$M_P$ & Mass (\mj) ']
    mpsuns=[r'$M_P$ & Mass (\msun) ']
    tcs=[r'$T_C$ & Time of conjunction (\bjdtdb) ']
    tts=[r'$T_T$ & Time of minimum projected separation (\bjdtdb) ']
    t0s=[r'$T_0$ & Optimal conjunction time (\bjdtdb) ']
    semimajors=[r'$a$ & Semi-major axis (AU) ']
    idegs=[r'$i$ & Inclination (Degrees) ']
    eccs=[r'$e$ & Eccentricity ']
    odegs=[r'$\omega_*$ & Argument of periastron (Degrees) ']
    teqs=[r'$T_{eq}$ & Equilibrium temperature (K) ']
    tcircs=[r'$\tau_{\rm circ}$ & Tidal circularization timescale (Gyr) ']
    ks=[r'$K$ & RV semi-amplitude (m/s) ']
    ps=[r'$R_P/R_*$ & Radius of planet in stellar radii  ']
    ars=[r'$a/R_*$ & Semi-major axis in stellar radii  ']
    deltas=[r'$\delta$ & $\left(R_P/R_*\right)^2$ ']
    Tdepths=[r'$Depth$ & \tess flux decrement at mid-transit ']
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
    rhops=[r'$\rho_P$ & Density (cgs) ']
    #logps=[r'']
    loggps=[r'$logg_P$ & Surface gravity  ']
    safronovs=[r'$\Theta$ & Safronov number ']
    faves=[r'$\fave$ & Incident flux (\fluxcgs) ']
    tps=[r'$T_P$ & Time of periastron (\bjdtdb) ']
    tss=[r'$T_S$ & Time of eclipse (\bjdtdb) ']
    tas=[r'$T_A$ & Time of ascending node (\bjdtdb) ']
    tds=[r'$T_D$ & Time of descending node (\bjdtdb) ']
    vcves=[r'$V_c/V_e$ ']
    ecosws=[r'$e\cos{\omega_*}$ & ']
    esinws=[r'$e\sin{\omega_*}$ & ']
    #secosws=[r'']
    #sesinws=[r'']
    msinis=[r'$M_P\sin i$ & Minimum mass (\mj) ']
    qs=[r'$M_P/M_*$ & Mass ratio  ']
    drs=[r'$d/R_*$ & Separation at mid-transit  ']
    pts=[r'$P_T$ & A priori non-grazing transit prob  ']
    ptgs=[r'$P_{T,G}$ & A priori transit prob  ']
    pss=[r'$P_S$ & A priori non-grazing eclipse prob ']
    psgs=[r'$P_{S,G}$ & A priori eclipse prob  ']
    
    
    for ii in range(len(target_list)):

        medians = grab_medians(target_list[ii], path, file_prefix=file_prefix)

        if Mstar ==True: 
            make_median_string(medians,'mstar',mstars)
        if Rstar == True:
            make_median_string(medians,'rstar',rstars)
        if RstarSED==True:
            make_median_string(medians,'rstarsed',rstarseds)
        if Lstar == True:
            make_median_string(medians,'lstar',lstars)
        if fbol == True:
            make_median_string(medians,'fbol',fbols)
        if rhostar == True:
            make_median_string(medians,'rhostar',rhostars)
        if logg == True:
            make_median_string(medians,'logg',loggs)
        if teff == True:
            make_median_string(medians,'teff',teffs)
        if teffsed == True:
            make_median_string(medians,'teffsed',teffseds)
        if feh == True:
            make_median_string(medians,'feh',fehs)
        if initfeh == True:
            make_median_string(medians,'initfeh',initfehs)
        if age == True:
            make_median_string(medians,'age',ages)
        if eep == True:
            make_median_string(medians,'eep',eeps)
        if logmstar == True:
            make_median_string(medians,'logmstar',logmstars)
        if Av == True:
            make_median_string(medians,'Av',avs)
        if errscale == True:
            make_median_string(medians,'errscale',errscales)
        if parallax==True:
            make_median_string(medians,'parallax',plaxes)
        if distance==True:
            make_median_string(medians,'distance',dists)
            
        if period == True:
            make_median_string(medians,'Period',periods)
        if rp == True:
            make_median_string(medians,'rp',rps)
        if mp==True:
            make_median_string(medians,'mp',mps)
        if mpsun == True:
            make_median_string(medians,'mpsun',mpsuns)
        if tc == True:
            make_median_string(medians,'tc',tcs)
        if tt == True:
            make_median_string(medians,'tt',tts)
        if t0 == True:
            make_median_string(medians,'t0',t0s)
        if semimajor == True:
            make_median_string(medians,'a',semimajors)
        if ideg == True:
            make_median_string(medians,'ideg',idegs)
        if ecc == True:
            make_median_string(medians,'e',eccs)
        if omegadeg == True:
            make_median_string(medians,'omegadeg',odegs)
        if teq==True:
            make_median_string(medians,'teq',teqs)
        if tcirc == True:
            make_median_string(medians,'tcirc',tcircs)
        if k == True:
            make_median_string(medians,'k',ks)
        if p == True:
            make_median_string(medians,'p',ps)
        if ar==True:
            make_median_string(medians,'ar',ars)          
        if delta ==True:
            make_median_string(medians,'delta',deltas)
        if depthTESS==True:
            make_median_string(medians,'depth_TESS',Tdepths)
        if tau==True:
            make_median_string(medians,'tau',taus)
        if t14==True:
            make_median_string(medians,'t14',t14s)
        if tfwhm==True:
            make_median_string(medians,'tfwhm',tfwhms)
        if b==True:
            make_median_string(medians,'b',bs)
        if cosi==True:
            make_median_string(medians,'cosi',cosis)
        if bs==True:
            make_median_string(medians,'bs',bss)
        if taus==True:
            make_median_string(medians,'taus',tauss)
        if ts14==True:
            make_median_string(medians,'t14s',ts14s)
        if tfwhms==True:
            make_median_string(medians,'tfwhms',tfwhmss)
        if depth25==True:
            make_median_string(medians,'eclipsedepth25',depth25s)
        if depth50==True:
            make_median_string(medians,'eclipsedepth50',depth50s)
        if depth75==True:
            make_median_string(medians,'eclipsedepth75',depth75s)   
        if rhop==True:
            make_median_string(medians,'rhop',rhops)
        #if logp==True:
        #    make_median_string(medians,'logp',logps)
        if loggp==True:
            make_median_string(medians,'loggp',loggps)
        if safronov==True:
            make_median_string(medians,'safronov',safronovs)
        if fave==True:
            make_median_string(medians,'fave',faves)        
        if tp==True:
            make_median_string(medians,'tp',tps)
        if ts==True:
            make_median_string(medians,'ts',tss)
        if ta==True:
            make_median_string(medians,'ta',tas)
        if td==True:
            make_median_string(medians,'td',tds)
        if vcve==True:
            make_median_string(medians,'vcve',vcves)
        if ecosw==True:
            make_median_string(medians,'ecosw',ecosws)
        if esinw==True:
            make_median_string(medians,'esinw',esinws)
        #if secosw==True:
        #    make_median_string(medians,'secosw',secosws)
        #if sesinw==True:
        #    make_median_string(medians,'sesinw',sesinws)
        if msini==True:
            make_median_string(medians,'msini',msinis)
        if q==True:
            make_median_string(medians,'q',qs)
        if dr==True:
            make_median_string(medians,'dr',drs)
        if pt==True:
            make_median_string(medians,'pt',pts)
        if ptg==True:
            make_median_string(medians,'ptg',ptgs)
        if ps==True:
            make_median_string(medians,'ps',pss)
        if psg==True:
            make_median_string(medians,'psg',psgs)

    colstring = 'lc'
    namestring = ''
    
    for ii in range(len(target_list)):
        colstring+='c'
        namestring += (' & \colhead{' + target_list[ii] + '}')

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
    r'\tablecolumns{' + str(len(target_list) + 2) + '}'+'\n'+
    r'\tablehead{& ' + namestring + '}'+'\n'+
    r'\startdata'+'\n'+
    #r'\hline \\' + '\n' + 
    #r'\hline \\' + '\n' + 
    r'\multicolumn{' + str(len(target_list) + 2) + r'}{l}{\textbf{Priors}:} \\' + '\n' +
    r'$\pi$ (Gaussian) & Gaia Parallax (mas) & \\'  + '\n' +
    r'$[{\rm Fe/H}]$ (Gaussian) & Metallicity (dex) & \\' + '\n' +
    r'$A_V$ (upper limit) & V-band extinction (mag) & \\' + '\n' + 
    r'$D_T$ (Gaussian) & Dilution in \tess & \\' + '\n' +
    r'\hline' + '\n' +               
    r'\multicolumn{' + str(len(target_list) + 2) + r'}{l}{\textbf{Stellar Parameters}:} \\' + '\n' )
    #r'\smallskip\\\multicolumn{2}{l}{Stellar Parameters:}&\smallskip\\'+'\n')

            
        if Mstar ==True:
            write(mstars,fout)
        if Rstar==True:
            write(rstars,fout)
        if RstarSED==True:
            write(rstarseds,fout)
        if Lstar == True:
            write(lstars,fout)
        if fbol == True:
            write(fbols,fout)
        if rhostar == True:
            write(rhostars,fout)
        if logg == True:
            write(loggs,fout)
        if teff == True:
            write(teffs,fout)
        if teffsed == True:
            write(teffseds,fout)
        if feh == True:
            write(fehs,fout)
        if initfeh == True:
            write(initfehs,fout)
        if age == True:
            write(ages,fout)
        if eep == True:
            write(eeps,fout)
        if logmstar == True:
            write(logmstars,fout)
        if Av == True:
            write(avs,fout)
        if errscale == True:
            write(errscales,fout)
        if parallax==True:
            write(plaxes,fout)
        if distance==True:
            write(dists,fout)
        
        #fout.write(r'\smallskip\\\multicolumn{2}{l}{Planetary Parameters:}&b\smallskip\\' + '\n')
        fout.write(r'\multicolumn{' + str(len(target_list) + 2) + r'}{l}{\textbf{Planetary Parameters}:} \\' + '\n')
        if period == True:
            write(periods,fout)
        if rp==True:
            write(rps,fout)
        if mp==True:
            write(mps,fout)
        if mpsun == True:
            write(mpsuns,fout)
        if tc == True:
            write(tcs,fout)
        if tt == True:
            write(tts,fout)
        if t0 == True:
            write(t0s,fout)
        if semimajor == True:
            write(semimajors,fout)
        if ideg == True:
            write(idegs,fout)
        if ecc == True:
            write(eccs,fout)
        if omegadeg == True:
            write(odegs,fout)
        if teq==True:
            write(teqs,fout)
        if tcirc == True:
            write(tcircs,fout)
        if k == True:
            write(ks,fout)
        if p == True:
            write(ps,fout)
        if ar==True:
            write(ars,fout)
        if delta ==True:
            write(deltas,fout)
        if depthTESS==True:
            write(Tdepths,fout)
        if tau==True:
            write(taus,fout)
        if t14==True:
            write(t14s,fout)
        if tfwhm==True:
            write(tfwhms,fout)            
        if b==True:
            write(bs,fout)
        if cosi==True:
            write(cosis,fout)
        if bs==True:
            write(bss,fout)
        if taus==True:
            write(tauss,fout)
        if tfwhms==True:
            write(tfwhmss,fout)
        if depth25==True:
            write(depth25s,fout)
        if depth50==True:
            write(depth50s,fout)
        if depth75==True:
            write(depth75s,fout)
        if rhop==True:
            write(rhops,fout)
        #if logp==True:
        #    make_median_string(medians,'logp',logps)
        if loggp==True:
            write(loggps,fout)
        if safronov==True:
            write(safronovs,fout)
        if fave==True:
            write(faves,fout)
        if tp==True:
            write(tps,fout)
        if ts==True:
            write(tss,fout)
        if ts14==True:
            write(ts14s,fout)               
        if ta==True:
            write(tas,fout)
        if td==True:
            write(tds,fout)
        if vcve==True:
            write(vcves,fout)
        if ecosw==True:
            write(ecosws,fout)
        if esinw==True:
            write(esinws,fout)
        #if secosw==True:
        #    make_median_string(medians,'secosw',secosws)
        #if sesinw==True:
        #    make_median_string(medians,'sesinw',sesinws)
        if msini==True:
            write(msinis,fout)
        if q==True:
            write(qs,fout)
        if dr==True:
            write(drs,fout)
        if pt==True:
            write(pts,fout)
        if ptg==True:
            write(ptgs,fout)
        if ps==True:
            write(pss,fout)
        if psg==True:
            write(psgs,fout)
        
        # conclude with \enddata at the bottom of the input .tex file
        fout.write(r'\enddata' + '\n')