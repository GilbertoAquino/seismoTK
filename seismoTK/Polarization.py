class Polarization:
    def __init__(self):
        from pathlib import Path
        BASE_DIR = Path(__file__).resolve().parent
        self.__rootfort=(str(BASE_DIR)+'/Pol/prog/')
        self.libsac="/usr/local/sac/lib/libsacio.a"
        self.fort = "gfortran"
        self.Z=None
        self.E=None
        self.N=None
        self.pow=3
        self.wlenf=20
        self.dopm=0.6
        self.cycle=3
        self.f1=None
        self.f2=None
        self.nflen=3
        self.nfr=None
        self.nsp=None
        self.ave="ave"
        self.remove_npts=0
        self.max_npts=8192
    
    def __version__(self):
        import os
        os.system("cat "+self.__rootfort+"version.txt")
        print()

    def __compile__(self):
        import os
        os.system(self.fort+" -fno-pie -no-pie "+self.__rootfort+"polfre_s1.66el.f -o "+self.__rootfort+"/polfre_s1.66el "+self.libsac)
    
    def get_pol(self):
        import os
        import pandas as pd
        self.Z.write("V.sac")
        self.E.write("E.sac")
        self.N.write("N.sac")
        par1=" wlenf="+str(self.wlenf)+" pow="+str(self.pow)+" "+self.ave+" dopm="+str(self.dopm)+" nsp="+str(self.nsp)+" "
        par2=" f1="+str(self.f1)+" f2="+str(self.f2)+" nflen="+str(self.nflen)+" cycle="+str(self.cycle)+" nfr="+str(self.nfr)+""
        os.system(self.__rootfort+"polfre_s1.66el V.sac N.sac E.sac "+par1+par2+" flog")
        os.system("awk 'NR>1{print $1,$2,$3,$4,$5}' azi_dopm.asc > temp.asc")
        self.Pol=pd.read_table("temp.asc",sep=" ",header=None)
        self.Pol=self.Pol.rename(columns={0:"BAZ",1:"FREQ",2:"DOP",3:"TIME",4:"LIN"})
        os.system("rm V.sac N.sac E.sac azi_dopm.asc temp.asc")

    def read_station(self,data_root,station):
        from obspy import read
        try:
            self.Z=read(str(data_root)+str(station)+".V.sac")
            self.E=read(str(data_root)+str(station)+".T.sac")
            self.N=read(str(data_root)+str(station)+".R.sac")
        except:
            print("No rotated files. Using normal files")
            self.Z=read("../"+str(data_root)+str(station)+".V.sac")
            self.E=read("../"+str(data_root)+str(station)+".E.sac")
            self.N=read("../"+str(data_root)+str(station)+".N.sac")

    def preprocessing(self,dec=1,remove_fp=0):
        if self.Z != None or self.N != None or self.E != None:
            self.Z.decimate(int(dec))
            self.E.decimate(int(dec))
            self.N.decimate(int(dec))
            self.Z[0].data=self.Z[0].data[int(remove_fp):len(self.Z[0].data)]
            self.N[0].data=self.N[0].data[int(remove_fp):len(self.N[0].data)]
            self.E[0].data=self.E[0].data[int(remove_fp):len(self.E[0].data)]
        else:
            raise ValueError('No data in components')
    
    def plot_tz(self,x,y,z,ranges=[None,None],dz=20,s=50,log=True,set_limits=[[None,None],[0,None]],xlabel='Time [s]',ylabel='Freq [Hz]',zlabel='BAZ [DEG]',pshow = True,colormap='hsv'):
        import matplotlib.pyplot as plt
        import numpy as np
        if None in ranges:
            ranges[0] = z.min()
            ranges[1] = z.max()
        if None in set_limits[0] or None in set_limits[1]:
            set_limits[0][0] = y.min()
            set_limits[0][1] = y.max()
            set_limits[1][1] = x.max()
        fig, ax = plt.subplots(2,1,sharex=True,figsize=(9,16),gridspec_kw={'height_ratios': [1, 4]})
        cm = plt.cm.get_cmap(colormap)
        sc = ax[1].scatter(x,y, c=z,s=s, vmin=ranges[0], vmax=ranges[1], cmap=cm,zorder=100)
        cbarticks=np.arange(ranges[0],ranges[1]+1,dz)
        cbar=plt.colorbar(sc,ticks=cbarticks,orientation='horizontal',pad=0.1)
        cbar.set_label(zlabel)
        cbar.ax.tick_params(labelsize=8)
        plt.subplots_adjust(hspace=0)
        if log == True:
            ax[1].set_yscale("log")
        ax[1].set_ylim(set_limits[0])
        ax[1].set_xlim(set_limits[1])
        plt.xlabel(xlabel)
        ax[1].set_ylabel(ylabel)
        time=np.arange(0,len(self.Z[0].data)/self.Z[0].stats.sampling_rate,1/self.Z[0].stats.sampling_rate)
        ax[0].plot(time,self.Z[0].data)
        if pshow:
            plt.show()
    
    def plot_tBazDop(self,x,y,z,dop,ranges=[[None,None],[None,None]],s=2,log=True,set_limits=[[None,None],[0,None]],xlabel='Time [s]',ylabel='Freq [Hz]',zlabel='BAZ [DEG]',pshow = True,colormap1='hsv',colormap2='gist_stern_r',seismogram=None):
        """
            x-> pandas.Series: x.axis; Time
            y-> pandas.Series: y.axis; Frequency
            z-> pandas.Series: first color scatter plot; Baz
            dop-> pandas.Series: second color scatter plot; DOP

        """
        if None in ranges[0] or None in ranges[1]:
            ranges[0][0] = z.min()
            ranges[0][1] = z.max()
            ranges[1][0] = dop.min()
            ranges[1][1] = dop.max()
        if None in set_limits[0] or None in set_limits[1]:
            set_limits[0][0] = y.min()
            set_limits[0][1] = y.max()
            set_limits[1][1] = x.max()
        import matplotlib.pyplot as plt
        import numpy as np
        fig, ax = plt.subplots(3,1,sharex=True,figsize=(9,16),gridspec_kw={'height_ratios': [1, 2, 2]})
        cm = plt.cm.get_cmap(colormap1)
        cm2 = plt.cm.get_cmap(colormap2)
        sc = ax[1].scatter(x,y, c=z,s=s, vmin=ranges[0][0], vmax=ranges[0][1], cmap=cm,zorder=100)
        sc2 = ax[2].scatter(x,y, c=dop,s=s, vmin=ranges[1][0], vmax=ranges[1][1], cmap=cm2,zorder=100)
        cbarticks=np.arange(ranges[0][0],ranges[0][1]+1,(ranges[0][1]-ranges[0][0])/3)
        cbarticks2=np.arange(ranges[1][0],ranges[1][1]+1,(ranges[1][1]-ranges[1][0])/3)
        cbaxes = fig.add_axes([0.12, 0.05, 0.25, 0.02])
        cbaxes2 = fig.add_axes([0.62, 0.05, 0.25, 0.02])  
        cbar=plt.colorbar(sc,cax = cbaxes,ticks=cbarticks,orientation='horizontal',pad=-0.01,shrink=0.3)
        cbar2=plt.colorbar(sc2,cax = cbaxes2,ticks=cbarticks2,orientation='horizontal',shrink=0.3)
        cbar.set_label(zlabel)
        cbar.ax.tick_params(labelsize=8)
        cbar2.set_label('DOP')
        cbar2.ax.tick_params(labelsize=8)
        plt.subplots_adjust(hspace=0,top=0.995,bottom=0.1)
        if log == True:
            ax[1].set_yscale("log")
            ax[2].set_yscale("log")
        ax[1].set_ylim(set_limits[0])
        ax[1].set_xlim(set_limits[1])
        ax[2].set_ylim(set_limits[0])
        ax[2].set_xlim(set_limits[1])
        ax[2].set_xlabel(xlabel)
        ax[1].set_ylabel(ylabel)
        ax[2].set_ylabel(ylabel)
        if seismogram == None:
            seismogram = self.Z[0]
        time=np.arange(0,len(seismogram.data)/seismogram.stats.sampling_rate,1/seismogram.stats.sampling_rate)
        ax[0].plot(time,seismogram.data,lw=0.4,color='black')
        ax[0].legend([seismogram.stats.station])
        if pshow:
            plt.show()
    
    def filter_data(self,column,min=None,max=None):
        if min != None:
            self.Pol = self.Pol[self.Pol[column] > min]
        if max != None:
            self.Pol = self.Pol[self.Pol[column] < max]