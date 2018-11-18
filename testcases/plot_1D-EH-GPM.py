"""from plot_maecs_omexdia import plot_maecs_omexdia 
    (time,z,dz,data,datanames)=plot_maecs_omexdia()"""
# python /home/onur/WORK/projects/GB/maecs/1d/BO//core/plot_maecs_medmac.py /home/onur/WORK/projects/GB/maecs/1d/BO//gotm-NS_scenarios/GB25_NS_m1024_p0930_medmac_base 1 3
    
from pylab import *
import netCDF4 as nc4
import netcdftime
import sys

def plot_maecs_Bpoolx2_phy():
    
    """from plot_maecs_omexdia import plot_maecs_omexdia 
    (time,z,dz,data,datanames)=plot_maecs_omexdia()"""    
    
    plottype='wc_mean' #wc_int, wc_mean,middlerow
    colmap='viridis'
    #import pdb
    if len(sys.argv) < 2: #this means no arguments were passed      
      fname='/home/onur/setups/test-BGCmodels/gpm-eh/1d-hlg/test/1D-HLG_GPM-EH'
      disp('plotting default file:'+fname)
    else:
      disp('plotting file specified:'+sys.argv[1])
      fname=sys.argv[1]
    
    if len(sys.argv)<3: #no second argument was passed
      plotsed=1
    else:
      plotsed=int(sys.argv[2])
    disp('plotsed:'+str(plotsed))
      
    if len(sys.argv)<4: #no third argument was passed
      numyears=1
    else: 
      numyears=int(sys.argv[3])
    disp('plotting last '+str(numyears)+' year of the simulation')
    
    if len(sys.argv) < 5: #this means no arguments were passed
	varnames= [ 'temp','EH_abioP_O2_percSat','EH_abioP_DINH4',
                'EH_abioP_DINO3', 'EH_abioP_DIP','EH_abioP_DISi',
                'EH_abioP_DOC','EH_abioP_det1C', 'EH_abioP_det2C',
                'GPM_diat_C','GPM_nf_C','GPM_pha_C',
                'GPM_mixo_C','GPM_miczoo_C','GPM_meszoo_C'
                #,'total_NPPR_calculator_result','total_chlorophyll_calculator_result'
		      ]
	numcol=3.0
    else: 
	varnames=sys.argv[4].split(',')
	numcol=length(varnames)
	
    if plotsed: 
        #from plot_sediment import readsed
        #sediment variables
        pickled=0
        varnames2=['EH_abioP_air_o2o', 'EH_abioS_o2o_brm', 'EH_abioS_sed_nn2']
        numsedvars=len(varnames2)
        figuresize=(16,20)
    else:
        numsedvars=0
        figuresize=(20,10) #(25,15)
    
    #pelagic variables
    nc=nc4.Dataset(fname+'.nc')
    ncv=nc.variables
    #print('available maecs variables:')        
    #disp(ncv)
    
    depth=np.squeeze(ncv['z'][:])
    if len(depth.shape)==2:
        if all(depth[0]==depth[-1]):
            depth=depth[0,:]
        else:
            depth = depth[0, :]
            #raise(Warning('plotting for varying-depths is not implemented')) #eg., it needs interpolation on a fixed grid

    #time=ncv['time'][:]/86400.

    tv = nc.variables['time']
    utime=netcdftime.utime(tv.units)
    tvec=utime.num2date(tv[:])

    #crop the data for the time period requested
    years=np.array([tvec[ti].year for ti in range(0,len(tvec))])
    years2plot=range(years[-1]+1-numyears, years[-1]+1)
    yeari=np.where((years>=years2plot[0]) * (years<=years2plot[-1]))
    tvecC=tvec[yeari[0]]

    f=figure(figsize=figuresize)
    f.subplots_adjust(top=0.95,bottom=0.05,hspace=0.5, wspace=0.5)

    numvar=len(varnames)+numsedvars

    for i,varn in enumerate(varnames):       
        print varn
        ax=subplot(ceil(numvar/numcol),numcol,i+1)

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
            ax=subplot(ceil(numvar/numcol),numcol,i+len(varnames)+1)

            if not (varnames2[i] in ncv):
                ax.text(0.5,0.5,varnames2[i]+'\n\n was not found',
                        horizontalalignment='center',
                        verticalalignment='center',
                        transform=ax.transAxes)
                continue

            #print(ncv[varnames2[i]].long_name)
            title(ncv[varnames2[i]].long_name+' [$%s$]'%ncv[varnames2[i]].units,size=8.0)

            if ncv[varnames2[i]].shape[1] > 1:
                if plottype=='middlerow':
                    middlerow=int(round(len(ncv[varnames2[i]][1,:])/2))
                    dat=squeeze(ncv[varnames2[i]][:,middlerow])
                elif plottype=='wc_int':
                    dat=sum(squeeze(ncv[varnames2[i]][:,:]),1)
                elif plottype=='wc_mean':
                    dat=mean(squeeze(ncv[varnames2[i]][:,:]),1)
            else:
                dat=squeeze(ncv[varnames2[i]][:])

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
    savefig(fname+'_cont_'+plottype+'_phy.png')
    disp('python contour plot saved in: '+fname+'_cont'+plottype+'_phy.png')
    #show()
        
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
    plot_maecs_Bpoolx2_phy()
