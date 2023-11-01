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

def grab_medians(name, path, file_prefix = '.MIST.SED.'):
    '''
    name: name of the system. Ex: 'TOI-1855'
    path: path to the fit files. Ex: '/Users/jack/Research/pipelines/system_figure_pipeline/data/'
    file_prefix: prefix used in EXOFASTv2 output files
    '''

    median_names= ['parname', 'median_value', 'upper_error', 'lower_error']

    medians = pd.read_csv(path + name + file_prefix + 'median.csv', names=median_names, header=None, skiprows=1)
    return medians

def make_string(medians, param, array):
    param = param+'_0'

    if medians.parname.isin([param]).any() == True:
    
        uperr = float(medians.upper_error[medians.parname == param])
        loerr = float(medians.lower_error[medians.parname == param])
        if uperr==loerr:
            errstring = r'\pm ' + str(uperr)
        else:
            errstring = r'^{+'+ str(uperr) + '}_{-' + str(loerr) + '}'
        array.append('& $' + str(float(medians.median_value[medians.parname == param])))
        #if ii<len(paths)-1:
        array.append(errstring + '$ ')
    else:
        array.append('& ---')

def write(param_arr,file):
    for ii in param_arr:
        file.write(ii)
    file.write(r'\\'+'\n')

def med_table(target_list, path, file_prefix = '.MIST.SED.'):

    # set which parameters to output to the .tex file

    # Stellar parameters
    Mstar=True
    Rstar=True
    RstarSED=False
    Lstar=True
    fbol=True
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
    errscale=True
    parallax=True
    distance=True

    # Planetary parameters
    period=True
    rp=True
    mp=True
    mpsun=False #nor sure what the problem is with this
    tc=True
    tt=True
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
    delta=True
    depthTESS=True
    tau=True
    t14=True
    tfwhm=False
    b=True
    cosi=False #nor sure what the problem is with this
    bs=True
    taus=False
    ts14=True
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
    ts=True
    ta=False
    td=False
    vcve=False
    ecosw=True
    esinw=True
    #secosw=True
    #sesinw=True
    msini=False
    q=True
    dr=True
    pt=False
    ptg=False
    ps=False
    psg=False

    newfile = 'median_table.tex'

    i = 1
    while os.path.exists(newfile):
        newfile = 'median_table_' + str(i) + '.tex'
        i += 1

    
    mstars=[r'~~~~$M_*$\dotfill &Mass (\msun)\dotfill ']
    rstars=[r'~~~~$R_*$\dotfill &Radius (\rsun)\dotfill ']
    rstarseds=[r'~~~~$R_{*,SED}$\dotfill &Radius$^{1}$ (\rsun)\dotfill ']
    lstars=[r'~~~~$L_*$\dotfill &Luminosity (\lsun)\dotfill ']
    fbols=[r'~~~~$F_{Bol}$\dotfill &Bolometric Flux (cgs)\dotfill ']
    rhostars=[r'~~~~$\rho_*$\dotfill &Density (cgs)\dotfill ']
    loggs=[r'~~~~$\log{g}$\dotfill &Surface gravity (cgs)\dotfill ']
    teffs=[r'~~~~$T_{\rm eff}$\dotfill &Effective Temperature (K)\dotfill ']
    teffseds=[r'~~~~$T_{\rm eff,SED}$\dotfill &Effective Temperature$^{1}$ (K)\dotfill ']
    fehs=[r'~~~~$[{\rm Fe/H}]$\dotfill &Metallicity (dex)\dotfill ']
    initfehs=[r'~~~~$[{\rm Fe/H}]_{0}$\dotfill &Initial Metallicity$^{2}$ \dotfill ']
    ages=[r'~~~~$Age$\dotfill &Age (Gyr)\dotfill ']
    eeps=[r'~~~~$EEP$\dotfill &Equal Evolutionary Phase$^{3}$ \dotfill ']
    logmstars=[r'~~~~$\log{M_*}$\dotfill &Mass ($\log{\msun}$)\dotfill ']
    avs=[r'~~~~$A_V$\dotfill &V-band extinction (mag)\dotfill ']
    errscales=[r'~~~~$\sigma_{SED}$\dotfill & SED photometry error scaling\dotfill ']
    plaxes=[r'~~~~$\varpi$\dotfill &Parallax (mas)\dotfill ']
    dists=[r'~~~~$d$\dotfill &Distance (pc)\dotfill ']
    
    periods=[r'~~~~$P$\dotfill &Period (days)\dotfill ']
    rps=[r'~~~~$R_P$\dotfill &Radius (\rj)\dotfill ']
    mps=[r'~~~~$M_P$\dotfill &Mass (\mj)\dotfill ']
    mpsuns=[r'~~~~$M_P$\dotfill &Mass (\msun)\dotfill ']
    tcs=[r'~~~~$T_C$\dotfill &Time of conjunction$^{4}$ (\bjdtdb)\dotfill ']
    tts=[r'~~~~$T_T$\dotfill &Time of minimum projected separation$^{5}$ (\bjdtdb)\dotfill ']
    t0s=[r'~~~~$T_0$\dotfill &Optimal conjunction Time$^{6}$ (\bjdtdb)\dotfill ']
    semimajors=[r'~~~~$a$\dotfill &Semi-major axis (AU)\dotfill ']
    idegs=[r'~~~~$i$\dotfill &Inclination (Degrees)\dotfill ']
    eccs=[r'~~~~$e$\dotfill &Eccentricity \dotfill ']
    odegs=[r'~~~~$\omega_*$\dotfill &Argument of Periastron (Degrees)\dotfill ']
    teqs=[r'~~~~$T_{eq}$\dotfill &Equilibrium temperature$^{7}$ (K)\dotfill ']
    tcircs=[r'~~~~$\tau_{\rm circ}$\dotfill &Tidal circularization timescale (Gyr)\dotfill ']
    ks=[r'~~~~$K$\dotfill &RV semi-amplitude (m/s)\dotfill ']
    ps=[r'~~~~$R_P/R_*$\dotfill &Radius of planet in stellar radii \dotfill ']
    ars=[r'~~~~$a/R_*$\dotfill &Semi-major axis in stellar radii \dotfill ']
    deltas=[r'~~~~$\delta$\dotfill &$\left(R_P/R_*\right)^2$ \dotfill ']
    Tdepths=[r'~~~~$Depth_\rm{TESS}$\dotfill &Flux decrement at mid transit for \tess \dotfill ']
    taus=[r'~~~~$\tau$\dotfill &Ingress/egress transit duration (days)\dotfill ']
    t14s=[r'~~~~$T_{14}$\dotfill &Total transit duration (days)\dotfill ']
    tfwhms=[r'~~~~$T_{FWHM}$\dotfill &FWHM transit duration (days)\dotfill ']
    bs=[r'~~~~$b$\dotfill &Transit Impact parameter \dotfill ']
    cosis=[r'~~~~$\cos i$\dotfill &Cos of Inclination \dotfill']
    bss=[r'~~~~$b_S$\dotfill &Eclipse impact parameter \dotfill ']
    tauss=[r'~~~~$\tau_S$\dotfill &Ingress/egress eclipse duration (days)\dotfill ']
    ts14s=[r'~~~~$T_{S,14}$\dotfill &Total eclipse duration (days)\dotfill ']
    tfwhmss=[r'~~~~$T_{S,FWHM}$\dotfill &FWHM eclipse duration (days)\dotfill ']
    depth25s=[r'~~~~$\delta_{S,2.5\mu m}$\dotfill &Blackbody eclipse depth at 2.5$\mu$m (ppm)\dotfill ']
    depth50s=[r'~~~~$\delta_{S,5.0\mu m}$\dotfill &Blackbody eclipse depth at 5.0$\mu$m (ppm)\dotfill ']
    depth75s=[r'~~~~$\delta_{S,7.5\mu m}$\dotfill &Blackbody eclipse depth at 7.5$\mu$m (ppm)\dotfill ']
    rhops=[r'~~~~$\rho_P$\dotfill &Density (cgs)\dotfill ']
    #logps=[r'']
    loggps=[r'~~~~$logg_P$\dotfill &Surface gravity \dotfill ']
    safronovs=[r'~~~~$\Theta$\dotfill &Safronov Number \dotfill ']
    faves=[r'~~~~$\fave$\dotfill &Incident Flux (\fluxcgs)\dotfill ']
    tps=[r'~~~~$T_P$\dotfill &Time of Periastron (\bjdtdb)\dotfill ']
    tss=[r'~~~~$T_S$\dotfill &Time of eclipse (\bjdtdb)\dotfill ']
    tas=[r'~~~~$T_A$\dotfill &Time of Ascending Node (\bjdtdb)\dotfill ']
    tds=[r'~~~~$T_D$\dotfill &Time of Descending Node (\bjdtdb)\dotfill ']
    vcves=[r'~~~~$V_c/V_e$\dotfill ']
    ecosws=[r'~~~~$e\cos{\omega_*}$\dotfill & \dotfill']
    esinws=[r'~~~~$e\sin{\omega_*}$\dotfill & \dotfill']
    #secosws=[r'']
    #sesinws=[r'']
    msinis=[r'~~~~$M_P\sin i$\dotfill &Minimum mass (\mj)\dotfill ']
    qs=[r'~~~~$M_P/M_*$\dotfill &Mass ratio \dotfill ']
    drs=[r'~~~~$d/R_*$\dotfill &Separation at mid transit \dotfill ']
    pts=[r'~~~~$P_T$\dotfill &A priori non-grazing transit prob \dotfill ']
    ptgs=[r'~~~~$P_{T,G}$\dotfill &A priori transit prob \dotfill ']
    pss=[r'~~~~$P_S$\dotfill &A priori non-grazing eclipse prob ']
    psgs=[r'~~~~$P_{S,G}$\dotfill &A priori eclipse prob \dotfill ']
    
    
    for ii in range(len(target_list)):

        medians = grab_medians(target_list[ii], path, file_prefix=file_prefix)

        if Mstar ==True: 
            make_string(medians,'mstar',mstars)
        if Rstar == True:
            make_string(medians,'rstar',rstars)
        if RstarSED==True:
            make_string(medians,'rstarsed',rstarseds)
        if Lstar == True:
            make_string(medians,'lstar',lstars)
        if fbol == True:
            make_string(medians,'fbol',fbols)
        if rhostar == True:
            make_string(medians,'rhostar',rhostars)
        if logg == True:
            make_string(medians,'logg',loggs)
        if teff == True:
            make_string(medians,'teff',teffs)
        if teffsed == True:
            make_string(medians,'teffsed',teffseds)
        if feh == True:
            make_string(medians,'feh',fehs)
        if initfeh == True:
            make_string(medians,'initfeh',initfehs)
        if age == True:
            make_string(medians,'age',ages)
        if eep == True:
            make_string(medians,'eep',eeps)
        if logmstar == True:
            make_string(medians,'logmstar',logmstars)
        if Av == True:
            make_string(medians,'Av',avs)
        if errscale == True:
            make_string(medians,'errscale',errscales)
        if parallax==True:
            make_string(medians,'parallax',plaxes)
        if distance==True:
            make_string(medians,'distance',dists)
            
        if period == True:
            make_string(medians,'Period',periods)
        if rp == True:
            make_string(medians,'rp',rps)
        if mp==True:
            make_string(medians,'mp',mps)
        if mpsun == True:
            make_string(medians,'mpsun',mpsuns)
        if tc == True:
            make_string(medians,'tc',tcs)
        if tt == True:
            make_string(medians,'tt',tts)
        if t0 == True:
            make_string(medians,'t0',t0s)
        if semimajor == True:
            make_string(medians,'a',semimajors)
        if ideg == True:
            make_string(medians,'ideg',idegs)
        if ecc == True:
            make_string(medians,'e',eccs)
        if omegadeg == True:
            make_string(medians,'omegadeg',odegs)
        if teq==True:
            make_string(medians,'teq',teqs)
        if tcirc == True:
            make_string(medians,'tcirc',tcircs)
        if k == True:
            make_string(medians,'k',ks)
        if p == True:
            make_string(medians,'p',ps)
        if ar==True:
            make_string(medians,'ar',ars)          
        if delta ==True:
            make_string(medians,'delta',deltas)
        if depthTESS==True:
            make_string(medians,'depth_TESS',Tdepths)
        if tau==True:
            make_string(medians,'tau',taus)
        if t14==True:
            make_string(medians,'t14',t14s)
        if tfwhm==True:
            make_string(medians,'tfwhm',tfwhms)
        if b==True:
            make_string(medians,'b',bs)
        if cosi==True:
            make_string(medians,'cosi',cosis)
        if bs==True:
            make_string(medians,'bs',bss)
        if taus==True:
            make_string(medians,'taus',tauss)
        if ts14==True:
            make_string(medians,'t14s',ts14s)
        if tfwhms==True:
            make_string(medians,'tfwhms',tfwhmss)
        if depth25==True:
            make_string(medians,'eclipsedepth25',depth25s)
        if depth50==True:
            make_string(medians,'eclipsedepth50',depth50s)
        if depth75==True:
            make_string(medians,'eclipsedepth75',depth75s)   
        if rhop==True:
            make_string(medians,'rhop',rhops)
        #if logp==True:
        #    make_string(medians,'logp',logps)
        if loggp==True:
            make_string(medians,'loggp',loggps)
        if safronov==True:
            make_string(medians,'safronov',safronovs)
        if fave==True:
            make_string(medians,'fave',faves)        
        if tp==True:
            make_string(medians,'tp',tps)
        if ts==True:
            make_string(medians,'ts',tss)
        if ta==True:
            make_string(medians,'ta',tas)
        if td==True:
            make_string(medians,'td',tds)
        if vcve==True:
            make_string(medians,'vcve',vcves)
        if ecosw==True:
            make_string(medians,'ecosw',ecosws)
        if esinw==True:
            make_string(medians,'esinw',esinws)
        #if secosw==True:
        #    make_string(medians,'secosw',secosws)
        #if sesinw==True:
        #    make_string(medians,'sesinw',sesinws)
        if msini==True:
            make_string(medians,'msini',msinis)
        if q==True:
            make_string(medians,'q',qs)
        if dr==True:
            make_string(medians,'dr',drs)
        if pt==True:
            make_string(medians,'pt',pts)
        if ptg==True:
            make_string(medians,'ptg',ptgs)
        if ps==True:
            make_string(medians,'ps',pss)
        if psg==True:
            make_string(medians,'psg',psgs)

    colstring = 'lc'
    namestring = ''
    
    for ii in range(len(target_list)):
        colstring+='c'
        namestring += (' & ' + target_list[ii])
    

    with open(newfile, 'w') as fout: 
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
    r'\startlongtable'+ '\n'+
    r'\begin{deluxetable*}{' + colstring +'}\n'+
    r'\clearpage'+'\n'+
    r'\newpage'+'\n'+    
    r'\centering' + '\n' +  
    r'\tablecaption{Median values and 68\% confidence intervals.}'+'\n'+
    r'\tablehead{\omit}' + '\n' +
    r'\tabletypesize{\tiny}' + '\n' +               
    #r'Parameter & Units & Values' +
    #r'\tablehead{\colhead{~~~Parameter} & \colhead{Units} & \colhead{Values}' +
    #\multicolumn{' +   str(len(k2list))+'}{c}{Values}}'+
    #'\n'+
    r'\startdata'+'\n'+
    #r'\hline \\' + '\n' + 
    #r'\hline \\' + '\n' + 
    r'\multicolumn{2}{l}{Priors:} ' + namestring + r'\\' + '\n' + ###########put k2 ids here
    r'\hline \\' + '\n' +
    r'Gaussian & $\pi$ Gaia Parallax (mas) \dotfill & \\'  + '\n' +
    r'Gaussian & $[{\rm Fe/H}]$ Metallicity (dex) \dotfill & \\' + '\n' +
    r'Upper Limit & $A_V$ V-band extinction (mag) \dotfill & \\' + '\n' + 
    r'Gaussian$^{\prime}$& $D_T$ Dilution in \tess\  \dotfill & \\' + '\n' +
    r'\hline\hline' + '\n' +               
    r'Parameter & Units & Values \\' + '\n' +
    r'\textbf{Stellar Parameters}: \\' + '\n' )
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
        fout.write(r'\textbf{Planetary Parameters}:\\' + '\n')
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
        #    make_string(medians,'logp',logps)
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
        #    make_string(medians,'secosw',secosws)
        #if sesinw==True:
        #    make_string(medians,'sesinw',sesinws)
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
        
            
        fout.write(r'\enddata' + '\n'
                       + r'\label{}' + '\n'
                       + r'\end{deluxetable*}')