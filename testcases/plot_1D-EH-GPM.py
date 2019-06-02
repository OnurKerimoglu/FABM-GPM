# terminal call:
# python plot_1D_GPM.py 1D-GPM.nc plot_air(0/1) plot_sediment(0/1) num_years_to_plot(counting_backwards_from_the_last)

from pylab import *
import netCDF4 as nc4
import netcdftime
import sys
import warnings

def plot_main():
    
    plottype='wc_mean' #wc_int, wc_mean,middlerow
    colmap='viridis'
    #import pdb
    if len(sys.argv) < 2: #this means no arguments were passed
      fname='/home/onur/setups/test-BGCmodels/gpm-eh/1D-40m/test_GPM-EH/2019_05_21_odemet2_n20_2013/varsto/1D-40m_GPM-EH_varsto_dm.nc'
      disp('plotting default file:'+fname)
    else:
      disp('plotting file specified:'+sys.argv[1])
      fname=sys.argv[1]
    
    if len(sys.argv)<3: #no second argument was passed
      plotair=1
    else:
      plotair=int(sys.argv[2])
    disp('plotair:'+str(plotair))

    if len(sys.argv)<4: #no second argument was passed
      plotsed=1
    else:
      plotsed=int(sys.argv[3])
    disp('plotsed:'+str(plotsed))

    if len(sys.argv)<5: #no third argument was passed
      numyears=1 #plot only the last year
    else: 
      numyears=int(sys.argv[4])
    disp('plotting last '+str(numyears)+' year of the simulation')

    numcol = 3.0
    basewidth = 8+numyears*8
    baseheight = 20

    if len(sys.argv) < 6: #this means no arguments were passed
      varnames= [ 'temp','nuh','attenuation_coefficient_of_photosynthetic_radiative_flux_calculator_result',
                'EH_abioP_O2_percSat','total_chlorophyll_calculator_result','total_NPPR_calculator_result',
                'EH_abioP_DINO3', 'EH_abioP_DIP','EH_abioP_DISi',
                'EH_abioP_DIC','EH_abioP_DOC','EH_abioP_bacC',
                'EH_abioP_DINH4','EH_abioP_det1C', 'EH_abioP_det2C',
                'GPM_diat_C','GPM_nf_C','GPM_pha_C',
                'GPM_mixo_C','GPM_miczoo_C','GPM_meszoo_C']
    else:
      varnames=sys.argv[5].split(',')
      numcol=len(varnames)

    if plotair:
      # from plot_sediment import readsed
      # sediment variables
      pickled = 0
      varnames_air = ['airt', 'u10', 'v10']
      numairvars = len(varnames_air)
      h_air=3
    else:
      numairvars = 0
      h_air=0
        
    if plotsed: 
      #from plot_sediment import readsed
      #sediment variables
      pickled=0
      varnames_sed=['EH_abioS_PON','EH_abioS_POP','EH_abioS_POSi',
                    'EH_abioP_air_o2o', 'EH_abioS_o2o_brm', 'EH_abioS_sed_nn2']
      numsedvars=len(varnames_sed)
      h_sed = 3
    else:
      numsedvars=0
      h_sed = 0
    
    #figuresize=(20,10) #(25,15)
    figuresize = (basewidth, baseheight+h_air+h_sed)  # (25,15)
    
    #pelagic variables
    nc=nc4.Dataset(fname)
    ncv=nc.variables
    #print('available maecs variables:')        
    #disp(ncv)
    
    z=np.squeeze(ncv['z'][:-1]) #depths at layer centers (fabm variables, temp, salt, etc)
    zi=np.squeeze(ncv['zi'][:-1]) #depths at layer interfaces (diffusivities, fluxes, etc)

    tv = nc.variables['time']
    utime=netcdftime.utime(tv.units)
    # the very last element might be the first step of the next year, so crop that out
    tvec=utime.num2date(tv[:-1])

    #crop the data for the time period requested
    years=np.array([tvec[ti].year for ti in range(0,len(tvec))])
    years2plot=range(years[-2]+1-numyears, years[-2]+1)
    yeari=np.where((years>=years2plot[0]) * (years<=years2plot[-1]))
    tvecC=tvec[yeari[0]]

    f=figure(figsize=figuresize)
    f.subplots_adjust(top=0.95,bottom=0.05,hspace=0.5, wspace=0.5)

    numvar=numairvars+len(varnames)+numsedvars

    if plotair:
        for i in xrange(0, numairvars):
            ax = subplot(ceil(numvar / numcol), numcol, i + 1)

            if not (varnames_air[i] in ncv):
                ax.text(0.5, 0.5, varnames_air[i] + '\n\n was not found',
                        horizontalalignment='center',
                        verticalalignment='center',
                        transform=ax.transAxes)
                continue

            # print(ncv[varnames_air[i]].long_name)
            title(ncv[varnames_air[i]].long_name + ' [$%s$]' % ncv[varnames_air[i]].units, size=8.0)

            if ncv[varnames_air[i]].shape[1] > 1:
                if plottype == 'middlerow':
                    middlerow = int(round(len(ncv[varnames_air[i]][1, :]) / 2))
                    dat = squeeze(ncv[varnames_air[i]][:, middlerow])
                elif plottype == 'wc_int':
                    dat = sum(squeeze(ncv[varnames_air[i]][:, :]), 1)
                elif plottype == 'wc_mean':
                    dat = mean(squeeze(ncv[varnames_air[i]][:, :]), 1)
            else:
                dat = squeeze(ncv[varnames_air[i]][:])

            # crop the data for the time period requested
            datC = dat[yeari[0]]

            if varnames_air[i]=='airt':
                datC=datC-273.15 #K to C

            ax.plot(tvecC, datC, 'r-')

            # x-axis
            format_date_axis(ax, [tvecC[0], tvecC[-1]])
            ax.xaxis.grid(color='k', linestyle=':', linewidth=0.5)
            xlabel('')

            # y-axis
            yt = ax.get_yticks()
            ytl = ax.get_yticklabels()
            ax.set_yticks([yt[0], yt[-1]])
            ax.set_yticklabels([str(yt[0]), str(yt[-1])])

    for i,varn in enumerate(varnames):       
        print varn
        
        if varn in ['nuh','nus']:
            depth=zi #depth at interfaces
        else:
            depth=z

        #if depth vector is 2-D (vary with time)
        if len(depth.shape)==2:
            # repeat the tvecC to obtain a matrix
            t=np.transpose(array([tvecC,]*depth.shape[1]))
        else:
            t=tvecC

        ax=subplot(ceil(numvar/numcol),numcol,numairvars+i+1)

        if (varnames[i]=='skip'):
            continue
        if (not (varnames[i] in ncv)):
            ax.text(0.5,0.5,varnames[i]+'\n\n was not found',
                    horizontalalignment='center',
                    verticalalignment='center',
                    transform=ax.transAxes)
            continue


        #pdb.set_trace()
        dat=squeeze(ncv[varnames[i]][:,:])

        #crop the data for the time period requested
        datC=dat[yeari[0],:]
        longname=ncv[varnames[i]].long_name
        shortname=longname.split('hzg_maecs ')[-1]
        #datC[datC<-1e10] = np.nan

        if (np.max(datC)-np.min(datC)<1e-10):
            ax.text(0.5,0.5,varnames[i]+'\n\n all: %3.2f'%np.max(datC),
                    horizontalalignment='center',
                    verticalalignment='center',
                    transform=ax.transAxes)
            continue
        else:
            units=ncv[varnames[i]].units
            if units in ['%']:
                title(shortname + ' [%s]'%units, size=10.0)
            elif units== 'Celsius':
                title('GOTM Temperature [$^oC$]', size=10.0)
            else:
                title(shortname + ' [$%s$]'%units, size=10.0)

            if len(z.shape) == 2:
                pcf = ax.contourf(t, depth, datC, cmap=plt.get_cmap(colmap))
            else:
                pcf=ax.contourf(tvecC,depth,transpose(datC),cmap=plt.get_cmap(colmap))

        #x-axis
        format_date_axis(ax,[tvecC[0], tvecC[-1]])
        ax.xaxis.grid(color='k',linestyle=':',linewidth=0.5)
        xlabel('')

        #y-axis
        #yt  = gca().get_yticks()
        #ytl = gca().get_yticklabels()
        #gca().set_yticks([yt[0],yt[-1]])
        #gca().set_yticklabels([str(yt[0]),str(yt[-1])])
        ylabel('depth [m]')

        cbar = colorbar(pcf, shrink=0.8)
        #cbar.solids.set_edgecolor("face")
        #draw()
     
    
    if plotsed:
        for i in xrange(0, numsedvars):
            ax=subplot(ceil(numvar/numcol),numcol,i+numairvars+len(varnames)+1)

            if not (varnames_sed[i] in ncv):
                ax.text(0.5,0.5,varnames_sed[i]+'\n\n was not found',
                        horizontalalignment='center',
                        verticalalignment='center',
                        transform=ax.transAxes)
                continue

            #print(ncv[varnames_sed[i]].long_name)
            title(ncv[varnames_sed[i]].long_name+' [$%s$]'%ncv[varnames_sed[i]].units,size=8.0)

            if ncv[varnames_sed[i]].shape[1] > 1:
                if plottype=='middlerow':
                    middlerow=int(round(len(ncv[varnames_sed[i]][1,:])/2))
                    dat=squeeze(ncv[varnames_sed[i]][:,middlerow])
                elif plottype=='wc_int':
                    dat=sum(squeeze(ncv[varnames_sed[i]][:,:]),1)
                elif plottype=='wc_mean':
                    dat=mean(squeeze(ncv[varnames_sed[i]][:,:]),1)
            else:
                dat=squeeze(ncv[varnames_sed[i]][:])

            #crop the data for the time period requested
            datC=dat[yeari[0]]

            ax.plot(tvecC,datC,'r-')

            #x-axis
            format_date_axis(ax,[tvecC[0], tvecC[-1]])
            ax.xaxis.grid(color='k',linestyle=':',linewidth=0.5)
            xlabel('')

            #y-axis
            yt  = ax.get_yticks()
            ytl = ax.get_yticklabels()
            ax.set_yticks([yt[0],yt[-1]])
            ax.set_yticklabels([str(yt[0]),str(yt[-1])])
        
        
    nc.close()
    savefig(fname.split('.nc')[0]+'_cont_'+plottype+'.png')
    disp('python contour plot saved in: '+fname.split('.nc')[0]+'_cont_'+plottype+'.png')
        
    #if plotsed: return time,z,dz,data,datanames


def format_date_axis(ax,tspan):
    ax.set_xlim(tspan[0], tspan[1])
    if diff(tspan)[0].days<63:
        ax.xaxis.set_major_locator(matplotlib.dates.WeekdayLocator(byweekday=matplotlib.dates.MO) )
        ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%b\n%d'))
    elif diff(tspan)[0].days<367:
        ax.xaxis.set_major_locator(matplotlib.dates.MonthLocator(bymonthday=1, interval=12) )
        ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%Y'))
        ax.xaxis.set_minor_locator(matplotlib.dates.MonthLocator(bymonthday=1, interval=2) )
        ax.xaxis.set_minor_formatter(matplotlib.dates.DateFormatter('%b'))
        ax.xaxis.set_tick_params(which='major', pad=15)
    elif diff(tspan)[0].days<732:
        ax.xaxis.set_major_locator(matplotlib.dates.MonthLocator(bymonthday=1, interval=12) )
        ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%Y'))
        ax.xaxis.set_minor_locator(matplotlib.dates.MonthLocator(bymonthday=1, interval=3) )
        ax.xaxis.set_minor_formatter(matplotlib.dates.DateFormatter('%b'))
        ax.xaxis.set_tick_params(which='major', pad=16)
    elif diff(tspan)[0].days<1466:
        ax.xaxis.set_major_locator(matplotlib.dates.MonthLocator(bymonthday=1, interval=12) )
        ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%Y'))
        ax.xaxis.set_minor_locator(matplotlib.dates.MonthLocator(bymonthday=1, interval=6) )
        ax.xaxis.set_minor_formatter(matplotlib.dates.DateFormatter('%b'))
        ax.xaxis.set_tick_params(which='major', pad=16)
    elif diff(tspan)[0].days<3655:
        ax.xaxis.set_major_locator(matplotlib.dates.MonthLocator(bymonthday=1, interval=12) )
        ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%Y'))
    elif diff(tspan)[0].days<9130: #25*365=9125
        ax.xaxis.set_major_locator(matplotlib.dates.MonthLocator(bymonthday=1, interval=60) )
        ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%Y'))
    else:
        ax.xaxis.set_major_locator(matplotlib.dates.MonthLocator(bymonthday=1, interval=120) )
        ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%Y'))

if __name__ == "__main__":
    # if you call this script from the command line (the shell) it will
    # run the 'main' function
    plot_main()

