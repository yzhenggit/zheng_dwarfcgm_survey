def eagle_model_color():
    import seaborn as sns
    bin_cmap = sns.color_palette('colorblind', 6)
    bin_colors = [bin_cmap[0], bin_cmap[3], bin_cmap[2]]
    bin_als = [0.3, 0.4, 0.3]
    bins = ['bin1', 'bin2', 'bin3']
    bin_ranges = [[10., 10.5], [10.5, 11], [11, 11.5]]
    # bin1 is M_halo=10^10.0-10.5,
    # bin2 is M_halo=10^10.5-11.0,
    # bin3 is M_halo=10^11.0-11.5

    bin_labels = [r'$M_{\rm 200m}=10^{10.0-10.5}~M_\odot$',
                   r'$M_{\rm 200m}=10^{10.5-11}~M_\odot$',
                   r'$M_{\rm 200m}=10^{11-11.5}~M_\odot$']
    return bins, bin_colors, bin_als, bin_labels, bin_ranges

def plot_obsdata(ax, ion, use1sig=False, ecolor='k', color='k', fs=16,
                 al=1, add_legend_byhand=False):
    import numpy as np
    ## for symbol sizes
    ss = [6, 5, 8, 7, 6, 6, 7] 

    ## read in ioin measurement data 
    from load_cgm_dwarf_mod import read_ion_measurements_logN
    ion_label, lib_det, lib_sat, lib_nondet = read_ion_measurements_logN(ion, use1sig=use1sig)

    ## get corresponding det/sat/nondet values
    ref_det = lib_det['ref']
    rn_det = lib_det['r/rvir']
    logN_det = lib_det['logN']
    elogN_det = lib_det['elogN']
    logmhalo_det = lib_det['logmhalo']

    if len(lib_sat)!=0: # lib_sat
        rn_sat = lib_sat['r/rvir']
        logN_sat = lib_sat['logN']
        ref_sat = lib_sat['ref']
        logmhalo_sat = lib_sat['logmhalo']
    if len(lib_nondet)!=0: # lib_nondet
        rn_nondet = lib_nondet['r/rvir']
        logN_nondet = lib_nondet['logN']
        ref_nondet = lib_nondet['ref']
        logmhalo_nondet = lib_nondet['logmhalo']
    file_label = lib_nondet['file_label']

    # model colors
    from load_cgm_dwarf_mod import eagle_model_color
    bins, bin_colors, bin_als, bin_labels, bin_ranges = eagle_model_color()

    # plot
    from load_cgm_dwarf_mod import plt_symbols_standard
    refs, symbols = plt_symbols_standard()

    for i, iref in enumerate(refs):
        # detection ########
        if len(lib_det) != 0: 
            ind_ref = ref_det == iref
            for ibin, ibin_range in enumerate(bin_ranges):
                ind_mass = np.all([logmhalo_det>=ibin_range[0], logmhalo_det<=ibin_range[1]], axis=0)
                ind_comb = np.all([ind_ref, ind_mass], axis=0)
                if len(ref_det[ind_comb]) > 0: 
                    ax.errorbar(rn_det[ind_comb], logN_det[ind_comb],
                                yerr=elogN_det[ind_comb], color=bin_colors[ibin], 
                                markeredgecolor='k', markeredgewidth=0.5,
                                fmt=symbols[i], markersize=ss[i], label=None)
        # non detection ########
        if len(lib_nondet)!=0:
            ind_ref = ref_nondet == iref
            for ibin, ibin_range in enumerate(bin_ranges):
                ind_mass = np.all([logmhalo_nondet>=ibin_range[0], logmhalo_nondet<=ibin_range[1]], axis=0)
                ind_comb = np.all([ind_ref, ind_mass], axis=0)
                if len(ref_nondet[ind_comb]) > 0:
                    if ion == 'HI':
                        uplims = np.asarray([0.25]*len(ref_nondet[ind_comb]))
                    else:
                        uplims = np.asarray([0.15]*len(ref_nondet[ind_comb]))
                    ax.errorbar(rn_nondet[ind_comb], logN_nondet[ind_comb], fmt=symbols[i],
                                 yerr=uplims, markerfacecolor='none', markeredgecolor=bin_colors[ibin],
                                 markeredgewidth=1.5, capthick=0.5,  barsabove=True,
                                 markersize=ss[i], label=None, lw=0.5, uplims=uplims, ecolor=bin_colors[ibin])
        # saturation ########
        if len(lib_sat)!=0:
            ind_ref = ref_sat == iref
            for ibin, ibin_range in enumerate(bin_ranges):
                ind_mass = np.all([logmhalo_sat>=ibin_range[0], logmhalo_sat<=ibin_range[1]], axis=0)
                ind_comb = np.all([ind_ref, ind_mass], axis=0)
                if len(ref_sat[ind_comb]) > 0: 
                    if ion == 'HI':
                        lolims = np.asarray([0.35]*len(ref_sat[ind_comb]))
                    else:
                        lolims = np.asarray([0.25]*len(ref_sat[ind_comb]))
                    ax.errorbar(rn_sat[ind_comb], logN_sat[ind_comb], fmt=symbols[i],
                         yerr=lolims, markerfacecolor=bin_colors[ibin], markeredgecolor='k',
                         markersize=ss[i], label=None, lw=0.5, markeredgewidth=0.5,
                         lolims=lolims, ecolor=bin_colors[ibin])

    ### add a legend box for HI ###
    if ion=='HI' and add_legend_byhand==True:
        dx, dy = 0.06, 0.2
        xx = 0.03
        yy = 12.9
        gap = 0.4
        ax.fill_between([xx, xx+dx], yy, yy+dy, color=bin_colors[0], alpha=0.8)
        ax.text(xx+dx*1.1, yy-dy/5, bin_labels[0], fontsize=fs-7)
        ax.fill_between([xx, xx+dx], yy-gap, yy-gap+dy, color=bin_colors[1], alpha=0.8)
        ax.text(xx+dx*1.1, yy-gap-dy/5, bin_labels[1], fontsize=fs-7)
        ax.fill_between([xx, xx+dx], yy-2*gap, yy-2*gap+dy, color=bin_colors[2], alpha=0.8)
        ax.text(xx+dx*1.1, yy-2*gap-dy/5, bin_labels[2], fontsize=fs-7)
    return ax, file_label

def read_ion_measurements_logN(ion, use1sig=False):
    """
    Read in logN and b data from the table, change the 1sig limit to 3sigma if use1sig=False (default)
    """
    from astropy.table import Table
    import numpy as np
    # gal_qso_final = Table.read('tables/DwarfsCGM_combined_literature_refined_{}.csv'.format(date), format='ascii')
    gal_qso_final = Table.read('data/zheng_dwarfcgm-survey_cut-mstar9.5-snr8-br200m.csv', format='ascii')
    gal_logmhalo = gal_qso_final['logM200m']
    # sightline with detection
    ind_eq = gal_qso_final['Wr_flg_'+ion] == '='
    # sightline with saturation, only in HI
    ind_sat = gal_qso_final['Wr_flg_'+ion] == '>'
    # sightline with upper limit
    ind_nondet = (gal_qso_final['Wr_flg_'+ion] == '<=(1sig)') & np.isfinite(gal_qso_final['Wr_'+ion]) # censored data

    from load_cgm_dwarf_mod import get_values_det_sat
    # detections
    lib_det = get_values_det_sat(gal_qso_final, ind_eq, ion)
    lib_det['logmhalo'] = gal_logmhalo[ind_eq]

    # saturation if available
    if len(gal_qso_final[ind_sat]) != 0:
        has_sat = True
        lib_sat = get_values_det_sat(gal_qso_final, ind_sat, ion)
        lib_sat['logmhalo'] = gal_logmhalo[ind_sat]
    else:
        has_sat = False
        lib_sat = {}

    # non detections
    from load_cgm_dwarf_mod import get_values_nondet
    if len(gal_qso_final[ind_nondet]) != 0:
        has_nondet = True
        lib_nondet = get_values_nondet(gal_qso_final, ind_nondet, ion, use1sig=use1sig)# use 3sigma is default
        lib_nondet['logmhalo'] = gal_logmhalo[ind_nondet]
    else:
        has_nondet = False
        lib_nondet = {}

    ion_label_lib = {'CII': r'CII 1334$\rm \AA$',
                     'SiII': r'SiII 1260$\rm \AA$',
                     'SiIII': r'SiIII 1206$\rm \AA$',
                     'SiIV': r'SiIV 1393$\rm \AA$',
                     'CIV': r'CIV 1548$\rm \AA$',
                     'HI': r'HI 1215$\rm \AA$'}
    ion_label = ion_label_lib[ion]

    return ion_label, lib_det, lib_sat, lib_nondet

def get_values_nondet(gal_qso_final, ind, ion, use1sig=False):
    import numpy as np
    ### mainly for step4 pymc3
    r_nondet = np.array(gal_qso_final['impact_para_kpc'][ind])
    rvir_nondet = np.array(gal_qso_final['R200m_kpc'][ind])
    rn_nondet = r_nondet/rvir_nondet
    logrn_nondet = np.log10(rn_nondet)
    ref_nondet = np.array(gal_qso_final['References'][ind])

    # log column density
    logN_1sig = np.array(gal_qso_final['logN_{}'.format(ion)][ind])
    logN_3sig = np.log10((10**logN_1sig) * 3.)

    # equivalent widths and log values
    wr_1sig = np.array(gal_qso_final['Wr_{}'.format(ion)][ind])
    logwr_1sig = np.log10(wr_1sig)

    wr_3sig = wr_1sig*3.
    logwr_3sig = np.log10(wr_3sig)

    if use1sig == True:
        logN_nondet = logN_1sig
        wr_nondet = wr_1sig
        logwr_nondet = logwr_1sig
        use_3sig = False
        file_label = 'use1sig'
    else:
        logN_nondet = logN_3sig
        wr_nondet = wr_3sig
        logwr_nondet = logwr_3sig
        file_label = 'use3sig'

    return {'r': r_nondet,
            'rvir': rvir_nondet,
            'r/rvir': rn_nondet,
            'log(r/rvir)': logrn_nondet,
            'ref': ref_nondet,
            'logN': logN_nondet,
            'EW': wr_nondet,
            'logEW': logwr_nondet,
            'file_label': file_label}

def get_values_det_sat(gal_qso_final, ind, ion):
    import numpy as np
    # mainly for step 4 pymc3
    # radius, normalized by impact parameters, and log values
    r_det = np.array(gal_qso_final['impact_para_kpc'][ind])
    rvir_det = np.array(gal_qso_final['R200m_kpc'][ind])
    rn_det = r_det/rvir_det
    logrn_det = np.log10(rn_det)
    ref_det = np.asarray(gal_qso_final['References'][ind])

    # log column density
    logN_det = np.array(gal_qso_final['logN_{}'.format(ion)][ind])
    elogN_det = np.array(gal_qso_final['elogN_{}'.format(ion)][ind])

    # equivalent widths and log values
    wr_det = np.array(gal_qso_final['Wr_{}'.format(ion)][ind])
    ewr_det = np.array(gal_qso_final['eWr_{}'.format(ion)][ind])

    logwr_det = np.log10(wr_det)
    elogwr_det = ewr_det/wr_det/np.log(10)

    return {'r': r_det,
            'rvir': rvir_det,
            'r/rvir': rn_det,
            'log(r/rvir)': logrn_det,
            'ref': ref_det,
            'logN': logN_det,
            'elogN': elogN_det,
            'EW': wr_det,
            'eEW': ewr_det,
            'logEW': logwr_det,
            'elogEW': elogwr_det}

def linear_model(x, slope, intercept):
    """
    The power law model in log log space
    """
    y = x*slope + intercept
    return y

def plot_posterior(ax, chainfile, p50, model_form, nlines=500, c='steelblue', plt_linear=False, al=0.1, zorder=1):
    """
    randomly draw a sample from the posterior distribution (n=nlines)
    and plot the result on the logN vs b plot
    """
    import pandas as pd
    import numpy as np
    flat_chain = pd.read_feather(chainfile)
    para1_chain = np.asarray(flat_chain.slope)
    para2_chain = np.asarray(flat_chain.intercept)

    lw_dd = 4

    # plot posterior
    # x = np.logspace(logrn_data.min(), logrn_data.max(), 50)
    x = np.logspace(-1.5, 0, 20)
    logx = np.log10(x)
    for ni, i in enumerate(np.random.randint(0, len(para1_chain), nlines)):
        logy = model_form(logx, para1_chain[i], para2_chain[i])
        y = 10**logy
        if plt_linear == True:
            if ni == 0:
                ax.plot(x, y, color=c, lw=lw_dd, alpha=0.1, label='Random posterior draws', zorder=zorder)
            else:
                ax.plot(x, y, color=c, lw=lw_dd, alpha=0.1, label=None, zorder=zorder)
        else:
            if ni == 0:
                ax.plot(x, logy, color=c, lw=lw_dd, alpha=0.1, label='Random posterior draws', zorder=zorder)
            else:
                ax.plot(x, logy, color=c, lw=lw_dd, alpha=0.1, label=None, zorder=zorder)

    # also plot the 50th percentile value
    p50_logy = model_form(logx, p50[0, 0], p50[0, 1])
    p50_y = 10**p50_logy
    if plt_linear == True:
        ax.plot(x, p50_y, color='k', lw=2, label='50th percentile solution')
    else:
        ax.plot(x, p50_logy, color='k', lw=2, label='50th percentile solution')

def get_percentile(chainfile, labels=['slope', 'intercept'], do_print=False):
    """
    get the 16th, 50th, and 84th percentile

    return:
    p50[0] = 50th
    p50[1] = -(50th-16th), neg error
    p50[2] = +(84th-50th), pos error
    """
    import numpy as np
    import pandas as pd
    from IPython.display import display, Math
    chain_data = pd.read_feather(chainfile)
    para1_chain = np.asarray(chain_data.slope)
    para2_chain = np.asarray(chain_data.intercept)
    flat_chain = np.asarray([para1_chain, para2_chain])

    ndim = 2
    p50 = np.zeros((3, ndim))
    for i in range(ndim):
        mcmc = np.percentile(flat_chain[i], [16, 50, 84])
        p50[0, i] = mcmc[1]
        p50[1, i] = -(mcmc[1]-mcmc[0])
        p50[2, i]= mcmc[2]-mcmc[1]

        if do_print==True:
            q = np.diff(mcmc)
            txt = "\mathrm{{{3}}} = {0:.2f}_{{-{1:.2f}}}^{{+{2:.2f}}}"
            txt = txt.format(mcmc[1], q[0], q[1], labels[i])
            display(Math(txt))
    return p50

def plt_symbols_standard():
    refs = ['Bordoloi+2014', 'Johnson+2017',
            'Liang&Chen2014', 'Qu&Bregman2022', 'Zheng+2019b',
            'Zheng+2020b', 'NewObs/Archived']
    symbols = ['s', 'D', 'p', 'X', '*', '<', 'o']
    return refs, symbols

def ylims_set(ions):
    # decide the ylims for each ions, 3sigma and 1sigma version
    ylims_set = {'HI': {'use3sig':[11.9, 17.2],
                        'use1sig':[11.5, 17.2]},
                 #'CII': {'use3sig':[12.2, 15],
                 #        'use1sig':[11, 15.]},
                 'CII': {'use3sig':[11.5, 15],
                         'use1sig':[11, 15.]},
                 #'CIV': {'use3sig': [12., 14.7],
                 #        'use1sig': [11., 14.7]},
 		 'CIV': {'use3sig': [11.5, 14.7],
                         'use1sig': [11., 14.7]},
                 'SiII': {'use3sig': [11.2, 13.9],
                          'use1sig': [10, 14.1]},
                 'SiIII': {'use3sig': [11, 14.1],
                           'use1sig': [10, 14.1]},
                 'SiIV': {'use3sig': [11, 14.],
                          'use1sig': [10., 14]}
                }
    ylims = {}
    for ion in ions:
        ylims[ion] = ylims_set[ion]
    return ylims

def init_logN_rho_plot_2panel(fs=14, ions=['HI', 'CIV']):
    """
    initialize the logN b plot to put observational data points on
    only for 2 panels, default to HI and CIV
    """
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    mpl.rcParams['font.family'] = 'stixgeneral'

    from load_cgm_dwarf_mod import plt_symbols_standard
    refs, symbols = plt_symbols_standard()

    # ax1 for HI and ax2 for CIV
    fig = plt.figure(figsize=(7, 3.5))
    ax1 = fig.add_axes([0.1, 0.24, 0.4, 0.74])
    ax2 = fig.add_axes([0.58, 0.24, 0.4, 0.74])
    axes = [ax1, ax2]
    for ax in axes:
        ax.minorticks_on()
        ax.set_xlim(0, 1.02)
        ax.tick_params(labelsize=fs-4)
        ax.set_xlabel(r'$b/R_{\rm 200m}$', fontsize=fs)
    ax1.set_ylabel('log[N (cm$^{-2}$)]', fontsize=fs)

    ## add another axes for legend

    ax3 = fig.add_axes([0.03, 0.01, 0.9, 0.08])
    ax3.set_xlim(0, 1)
    ax3.set_ylim(0, 1)
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    ax3.spines['bottom'].set_visible(False)
    ax3.spines['left'].set_visible(False)
    ax3.get_xaxis().set_ticks([])
    ax3.get_yaxis().set_ticks([])
    ax3.set_facecolor('none')
    line1 = 0.7
    line2 = 0.23
    s_pp = 0.017 # point ga[]
    t_pp = 0.027 # text gap
    ll_gap = 0.2 # larger lagp between two columns
    start = 0.02

    # Bordoloi
    txt_fs = fs-8
    ss = 15
    i = 0
    line = line1
    ax3.scatter([start], [line], marker=symbols[i], color='k', s=ss)
    ax3.scatter([start+s_pp], [line], marker=symbols[i], edgecolor='k', facecolor='none', s=ss)
    ax3.text(start+t_pp, line-0.09, 'Bordoloi+2014,2018', fontsize=txt_fs)

    # johnson
    i = 1
    line = line2
    ax3.scatter([start], [line], marker=symbols[i], color='k', s=ss)
    ax3.scatter([start+s_pp], [line], marker=symbols[i], edgecolor='k', facecolor='none', s=ss)
    ax3.text(start+t_pp, line-0.09, 'Johnson+2017', fontsize=txt_fs)

    # LC2014
    i = 2
    line = line1
    ax3.scatter([start+ll_gap], [line], marker=symbols[i], color='k', s=ss)
    ax3.scatter([start+s_pp+ll_gap], [line], marker=symbols[i], edgecolor='k', facecolor='none', s=ss)
    ax3.text(start+t_pp+ll_gap, line-0.09, 'Liang & Chen (2014)', fontsize=txt_fs)

    # Qu & Bregman
    i = 3
    line = line2
    ax3.scatter([start+ll_gap], [line], marker=symbols[i], color='k', s=ss)
    ax3.scatter([start+s_pp+ll_gap], [line], marker=symbols[i], edgecolor='k', facecolor='none', s=ss)
    ax3.text(start+t_pp+ll_gap, line-0.09, 'Qu & Bregman (2022)', fontsize=txt_fs)

    # Zheng
    i = 5
    line = line1
    ax3.scatter([start+ll_gap*2+0.05], [line], marker=symbols[i], color='k', s=ss)
    ax3.scatter([start+s_pp+ll_gap*2+0.05], [line], marker=symbols[i],
                edgecolor='k', facecolor='none', s=ss)
    ax3.text(start+t_pp+ll_gap*2+0.05, line-0.09, 'Zheng+2020 (IC1613)', fontsize=txt_fs)

    # new wpairs
    i = 6
    line = line2
    ax3.scatter([start+ll_gap*2+0.05], [line], marker=symbols[i], color='k', s=ss)
    ax3.scatter([start+s_pp+ll_gap*2+0.05], [line], marker=symbols[i],
                edgecolor='k', facecolor='none', s=ss)
    ax3.text(start+t_pp+ll_gap*2+0.05, line-0.09, 'New Obs./Archival Search', fontsize=txt_fs)

    # det/nondet
    ax3.text(start+s_pp+ll_gap*3.7-.05, line1-0.09,
             'Open Symbols=Non-Detections', fontsize=txt_fs)
    ax3.text(start+s_pp+ll_gap*3.7-.05, line2-0.09,
             'Filled Symbols=Detections or Saturations', fontsize=txt_fs)

    ylims = ylims_set(ions)
    return fig, axes, ylims

def init_logN_rho_plot_6panel(fs=16, ions=['HI', 'CII', 'CIV', 'SiII', 'SiIII', 'SiIV']):
    """
    initialize the logN b plot to put observational data points on
    """
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    mpl.rcParams['font.family'] = 'stixgeneral'

    fig = plt.figure(figsize=(10, 6.5))
    ww = 0.08
    hh = 0.62
    dw = 0.26
    dh = 0.33
    dw_gap = 0.06
    dh_gap = 0.2

    ax1 = fig.add_axes([ww, hh, dw, dh])
    ax2 = fig.add_axes([ww+dw+dw_gap, hh, dw, dh])
    ax3 = fig.add_axes([ww+(dw+dw_gap)*2, hh, dw, dh])
    ax4 = fig.add_axes([ww, hh-dw-dh_gap, dw, dh])
    ax5 = fig.add_axes([ww+dw+dw_gap, hh-dw-dh_gap, dw, dh])
    ax6 = fig.add_axes([ww+(dw+dw_gap)*2, hh-dw-dh_gap, dw, dh])
    axes = [ax1, ax2, ax3, ax4, ax5, ax6]

    for ax in axes:
        ax.minorticks_on()
        ax.set_xlim(0, 1.02)
        ax.tick_params(labelsize=fs-4)
    ax1.set_ylabel('log[N (cm$^{-2}$)]', fontsize=fs)
    ax4.set_ylabel('log[N (cm$^{-2}$)]', fontsize=fs)

    ## add another axes for legend
    from load_cgm_dwarf_mod import plt_symbols_standard
    refs, symbols = plt_symbols_standard()
    ss = [6, 5, 8, 7, 6, 6, 7] ## for symbol sizes

    ax7 = fig.add_axes([0.03, 0.01, 0.9, 0.08])
    ax7.set_xlim(0, 1)
    ax7.set_ylim(0, 1)
    ax7.spines['top'].set_visible(False)
    ax7.spines['right'].set_visible(False)
    ax7.spines['bottom'].set_visible(False)
    ax7.spines['left'].set_visible(False)
    ax7.get_xaxis().set_ticks([])
    ax7.get_yaxis().set_ticks([])
    ax7.set_facecolor('none')
    line1 = 0.55
    line2 = 0.12
    s_pp = 0.017 # point ga[]
    t_pp = 0.027 # text gap
    ll_gap = 0.2 # larger lagp between two columns
    start = 0.05

    # Bordoloi
    i = 0
    line = line1
    ax7.scatter([start], [line], marker=symbols[i], color='k')
    ax7.scatter([start+s_pp], [line], marker=symbols[i], edgecolor='k', facecolor='none')
    ax7.text(start+t_pp, line-0.09, 'Bordoloi+2014,2018', fontsize=fs-5)

    # johnson
    i = 1
    line = line2
    ax7.scatter([start], [line], marker=symbols[i], color='k')
    ax7.scatter([start+s_pp], [line], marker=symbols[i], edgecolor='k', facecolor='none')
    ax7.text(start+t_pp, line-0.09, 'Johnson+2017', fontsize=fs-5)

    # LC2014
    i = 2
    line = line1
    ax7.scatter([start+ll_gap], [line], marker=symbols[i], color='k')
    ax7.scatter([start+s_pp+ll_gap], [line], marker=symbols[i], edgecolor='k', facecolor='none')
    ax7.text(start+t_pp+ll_gap, line-0.09, 'Liang & Chen (2014)', fontsize=fs-5)

    # Qu & Bregman
    i = 3
    line = line2
    ax7.scatter([start+ll_gap], [line], marker=symbols[i], color='k')
    ax7.scatter([start+s_pp+ll_gap], [line], marker=symbols[i], edgecolor='k', facecolor='none')
    ax7.text(start+t_pp+ll_gap, line-0.09, 'Qu & Bregman (2022)', fontsize=fs-5)

    # Zheng
    i = 5
    line = line1
    ax7.scatter([start+ll_gap*2], [line], marker=symbols[i], color='k')
    ax7.scatter([start+s_pp+ll_gap*2], [line], marker=symbols[i], edgecolor='k', facecolor='none')
    ax7.text(start+t_pp+ll_gap*2, line-0.09, 'Zheng+2020 (IC1613)', fontsize=fs-5)

    # ne wpairs
    i = 6
    line = line2
    ax7.scatter([start+ll_gap*2], [line], marker=symbols[i], color='k')
    ax7.scatter([start+s_pp+ll_gap*2], [line], marker=symbols[i], edgecolor='k', facecolor='none')
    ax7.text(start+t_pp+ll_gap*2, line-0.09, 'New Obs./Archival Search', fontsize=fs-5)

    # det/nondet
    ax7.text(start+s_pp+ll_gap*3.7-.1, line1-0.09, 'Open Symbols=Non-Detections', fontsize=fs-5)
    ax7.text(start+s_pp+ll_gap*3.7-.1, line2-0.09, 'Filled Symbols=Detections or Saturations', fontsize=fs-5)

    ylims = ylims_set(ions)

    return fig, axes, ylims

def init_logN_rho_plot_3panel(fs=16, ions=['HI', 'CII', 'CIV']):
    """
    initialize the logN b plot to put observational data points on
    """
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    mpl.rcParams['font.family'] = 'stixgeneral'

    fig = plt.figure(figsize=(10, 3.5))
    ww = 0.08
    hh = 0.28
    dw = 0.26
    dh = 0.62
    dw_gap = 0.06

    ax1 = fig.add_axes([ww, hh, dw, dh])
    ax2 = fig.add_axes([ww+dw+dw_gap, hh, dw, dh])
    ax3 = fig.add_axes([ww+(dw+dw_gap)*2, hh, dw, dh])
    axes = [ax1, ax2, ax3]

    for ax in axes:
        ax.minorticks_on()
        ax.set_xlim(0, 1.02)
        ax.tick_params(labelsize=fs-4)
    ax1.set_ylabel('log[N (cm$^{-2}$)]', fontsize=fs)

    ## add another axes for legend
    from load_cgm_dwarf_mod import plt_symbols_standard
    refs, symbols = plt_symbols_standard()
    ss = [6, 5, 8, 7, 6, 6, 7] ## for symbol sizes

    ax_ll = fig.add_axes([0.05, 0.01, 0.9, 0.12])
    ax_ll.set_xlim(0, 1)
    ax_ll.set_ylim(0, 1)
    ax_ll.spines['top'].set_visible(False)
    ax_ll.spines['right'].set_visible(False)
    ax_ll.spines['bottom'].set_visible(False)
    ax_ll.spines['left'].set_visible(False)
    ax_ll.get_xaxis().set_ticks([])
    ax_ll.get_yaxis().set_ticks([])
    ax_ll.set_facecolor('none')
    line1 = 0.58
    line2 = 0.15
    s_pp = 0.017 # point ga[]
    t_pp = 0.027 # text gap
    ll_gap = 0.2 # larger lagp between two columns
    start = 0.05

    # Bordoloi
    i = 0
    line = line1
    ax_ll.scatter([start], [line], marker=symbols[i], color='k')
    ax_ll.scatter([start+s_pp], [line], marker=symbols[i], edgecolor='k', facecolor='none')
    ax_ll.text(start+t_pp, line-0.09, 'Bordoloi+2014,2018', fontsize=fs-5)

    # johnson
    i = 1
    line = line2
    ax_ll.scatter([start], [line], marker=symbols[i], color='k')
    ax_ll.scatter([start+s_pp], [line], marker=symbols[i], edgecolor='k', facecolor='none')
    ax_ll.text(start+t_pp, line-0.09, 'Johnson+2017', fontsize=fs-5)

    # LC2014
    i = 2
    line = line1
    ax_ll.scatter([start+ll_gap], [line], marker=symbols[i], color='k')
    ax_ll.scatter([start+s_pp+ll_gap], [line], marker=symbols[i], edgecolor='k', facecolor='none')
    ax_ll.text(start+t_pp+ll_gap, line-0.09, 'Liang & Chen (2014)', fontsize=fs-5)

    # Qu & Bregman
    i = 3
    line = line2
    ax_ll.scatter([start+ll_gap], [line], marker=symbols[i], color='k')
    ax_ll.scatter([start+s_pp+ll_gap], [line], marker=symbols[i], edgecolor='k', facecolor='none')
    ax_ll.text(start+t_pp+ll_gap, line-0.09, 'Qu & Bregman (2022)', fontsize=fs-5)

    # Zheng
    i = 5
    line = line1
    ax_ll.scatter([start+ll_gap*2], [line], marker=symbols[i], color='k')
    ax_ll.scatter([start+s_pp+ll_gap*2], [line], marker=symbols[i], edgecolor='k', facecolor='none')
    ax_ll.text(start+t_pp+ll_gap*2, line-0.09, 'Zheng+2020 (IC1613)', fontsize=fs-5)

    # ne wpairs
    i = 6
    line = line2
    ax_ll.scatter([start+ll_gap*2], [line], marker=symbols[i], color='k')
    ax_ll.scatter([start+s_pp+ll_gap*2], [line], marker=symbols[i], edgecolor='k', facecolor='none')
    ax_ll.text(start+t_pp+ll_gap*2, line-0.09, 'New Obs./Archival Search', fontsize=fs-5)

    # det/nondet
    ax_ll.text(start+s_pp+ll_gap*3.7-.1, line1-0.09, 'Open Symbols=Non-Detections', fontsize=fs-5)
    ax_ll.text(start+s_pp+ll_gap*3.7-.1, line2-0.09, 'Filled Symbols=Detections or Saturations', fontsize=fs-5)

    # decide the ylims for each ions, 3sigma and 1sigma version
    ylims = ylims_set(ions)
    return fig, axes, ylims
