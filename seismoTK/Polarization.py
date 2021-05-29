class Polarization:
    def __init__(self,f1=0.1,f2=1,nfr = 100,nsp = 8192,cycle = 3,wlenf=25,dopm=0.6,max_npts=8192):
        from pathlib import Path
        BASE_DIR = Path(__file__).resolve().parent
        self._rootfort=(str(BASE_DIR)+'/Pol/prog/')
        self.libsac="/usr/local/sac/lib/libsacio.a"
        self.fort = "gfortran"
        self.flags = "-fno-pie -no-pie"
        self.Z=None
        self.E=None
        self.N=None
        self.pow=3
        self.wlenf=wlenf
        self.dopm=dopm
        self.cycle=cycle
        self.f1=f1
        self.f2=f2
        self.nflen=3
        self.nfr=nfr
        self.nsp=nsp
        self.ave="ave"
        self.remove_npts=0
        self.max_npts=max_npts
    
    def version(self):
        import os
        os.system("cat "+self._rootfort+"version.txt")
        print()

    def compile(self):
        import os
        os.system(self.fort+" "+self.flags+" "+self._rootfort+"polfre_s1.66el.f -o "+self._rootfort+"/polfre_s1.66el "+self.libsac)
    
    def get_pol(self):
        import os
        import pandas as pd
        try:
            self.Z.write("V.sac")
            self.E.write("E.sac")
            self.N.write("N.sac")
            par1=" wlenf="+str(self.wlenf)+" pow="+str(self.pow)+" "+self.ave+" dopm="+str(self.dopm)+" nsp="+str(self.nsp)+" "
            par2=" f1="+str(self.f1)+" f2="+str(self.f2)+" nflen="+str(self.nflen)+" cycle="+str(self.cycle)+" nfr="+str(self.nfr)+""
            os.system(self._rootfort+"polfre_s1.66el V.sac N.sac E.sac "+par1+par2+" flog")
            os.system("awk 'NR>1{print $1,$2,$3,$4,$5}' azi_dopm.asc > temp.asc")
            self.Pol=pd.read_table("temp.asc",sep=" ",header=None)
            self.Pol=self.Pol.rename(columns={0:"BAZ",1:"FREQ",2:"DOP",3:"TIME",4:"LIN"})
            os.system("rm V.sac N.sac E.sac azi_dopm.asc temp.asc")
        except:
            os.system("rm V.sac N.sac E.sac azi_dopm.asc temp.asc")
            raise ValueError('get_pol failed')
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
        if self.Z != None and self.N != None and self.E != None:
            self.Z.decimate(int(dec))
            self.E.decimate(int(dec))
            self.N.decimate(int(dec))
            self.Z[0].data=self.Z[0].data[int(remove_fp):len(self.Z[0].data)]
            self.N[0].data=self.N[0].data[int(remove_fp):len(self.N[0].data)]
            self.E[0].data=self.E[0].data[int(remove_fp):len(self.E[0].data)]
        else:
            raise ValueError('No data in components')
    
    def plot_tz(self,x,y,z,ranges=[None,None],dz=20,s=50,log=True,set_limits=[[None,None],[0,None]],xlabel='Time [s]',ylabel='Freq [Hz]',zlabel='BAZ [DEG]',pshow = True,colormap='hsv',seismogram=None):
        import matplotlib.pyplot as plt
        import numpy as np
        if None in ranges:
            ranges[0] = z.min()
            ranges[1] = z.max()
        if None in set_limits[0] or None in set_limits[1]:
            set_limits[0][0] = y.min()
            set_limits[0][1] = y.max()
            set_limits[1][1] = x.max()
        if seismogram == None:
            seismogram = self.Z[0]
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
        ax[0].plot(time,seismogram.data,lw=0.4,color='black')
        try:
            ax[0].legend([seismogram.stats.station])
        except:
            pass
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
    
    def freq_baz(self,xmin=-180,xmax=180,ymin=None,ymax=None,xn=50,yn=20,xlabel='BAZ[deg]',ylabel='FREQ[Hz]'):
        import numpy as np
        import matplotlib.pyplot as plt
        x = self.Pol['BAZ'].to_numpy()
        y = self.Pol['FREQ'].to_numpy()
        h = np.zeros([yn,xn])
        n=0
        xy=np.zeros([yn*xn,3])
        x_axis = np.zeros(xn)
        y_axis = np. zeros(yn)
        if ymin == None:
            ymin = y.min()
        if ymax == None:
            ymax = y.max()
        if len(x) == len(y):
            for i in range(0,len(x)):
                if x[i] >= xmin and x[i] <= xmax and y[i] >= ymin and y[i] <= ymax:
                    xi = int(xn * (x[i] - xmin) / (xmax - xmin))
                    if xi >= xn: 
                        xi = xn -1
                    yi = int(yn * (y[i] - ymin) / (ymax - ymin))
                    if yi >= yn: 
                        yi = yn -1
                    h[yi,xi]=h[yi,xi] + 1
                    n=n+1
            xscale = (xmax -xmin) / xn
            yscale = (ymax -ymin) / yn
            if n > 0:
                percent = 100.0 / n
            else:
                percent = 0.0
            k=0
            for j in range(0,yn):
                for i in range(0,xn):
                    xy[k][0] = np.round(xmin+xscale * (i + 0.5),8)
                    xy[k][1] = np.round(ymin+yscale * (j + 0.5),8)
                    h[j,i] = percent * h[j,i]
                    xy[k][2] = percent * h[j,i]
                    k=k+1
            for i in range(0,xn):
                x_axis[i] = np.round(xmin+xscale * (i + 0.5),8)
            for j in range(0,yn):
                y_axis[j] = np.round(ymin+yscale * (j + 0.5),8)
            fig = plt.figure()
            ax = fig.add_subplot(111)
            plt.contourf(x_axis,y_axis,h)
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            plt.show()
        else:
            raise ValueError('REALLY? HOW THIS IS POSSIBLE? x and y length missmatch')

    def xyz2grd(self,x,y,z,xmin=None,xmax=None,ymin=None,ymax=None,zmin=None,zmax=None,xn=100,yn=20):
        import numpy as np
        if xmin == None:
            xmin = x.min()
        if ymin == None:
            ymin = y.min()
        if zmin == None:
            zmin = z.min()
        if xmax == None:
            xmax = x.max()
        if ymax == None:
            ymax = y.max()
        if zmax == None:
            zmax = z.max()
        h = np.zeros([yn,xn])
        x_axis = np.zeros(xn)
        y_axis = np. zeros(yn)
        if len(x) == len(y):
            for i in range(0,len(x)):
                if x[i] >= xmin and x[i] <= xmax and y[i] >= ymin and y[i] <= ymax:
                    xi = int(xn * (x[i] - xmin) / (xmax - xmin))
                    if xi >= xn: 
                        xi = xn -1
                    yi = int(yn * (y[i] - ymin) / (ymax - ymin))
                    if yi >= yn: 
                        yi = yn -1
                    h[yi,xi]=h[yi,xi] + z[i]
                    if h[yi,xi] > zmax:
                        h[yi,xi] = zmax
            xscale = (xmax -xmin) / xn
            yscale = (ymax -ymin) / yn
            for i in range(0,xn):
                x_axis[i] = np.round(xmin+xscale * (i + 0.5),8)
            for j in range(0,yn):
                y_axis[j] = np.round(ymin+yscale * (j + 0.5),8)
        return x_axis,y_axis,h

    def save_pol(self,name):
        self.Pol.to_csv(name)

def read_pol(name,S_Filter=False):
    import pandas as pd
    if S_Filter:
        from . import S_Filter
        a = S_Filter()
        a.Pol = pd.read_csv(name)
        return a
    print('HOLA')
    from . import Polarization
    a = Polarization()
    a.Pol = pd.read_csv(name)
    return a