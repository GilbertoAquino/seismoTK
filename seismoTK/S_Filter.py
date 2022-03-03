from matplotlib.colors import Colormap
from . import Polarization
class S_Filter(Polarization):
    def S(self):
        from scipy.interpolate import interp2d
        self.SF = self.Pol.drop(columns=["LIN","BAZ"])
        self.SF['TIME'].round(decimals=2)
        x,y,z=self.xyz2grd(self.SF["TIME"],self.SF["FREQ"],self.SF["DOP"],xn=3600,yn=100)
        self.S_Tr()
        inter = interp2d(x, y, z)
        dop_inter = inter(self.tsig,self.fsig)
        self.filter(dop_inter)
        

    def plotspec(self,psx, cmap, lofreq=None, hifreq=None, t1=None, t2=None):
        import matplotlib.pyplot as plt
        extent = [0,0,0,0]
        if t1 != None and t2 != None:
            extent[0] = t1
            extent[1] = t2
        if lofreq != None:
            extent[2] = lofreq
        if hifreq != None:
            extent[3] = hifreq
        sc = plt.imshow(psx, cmap=cmap, extent=extent, aspect='auto', origin='lower')
        plt.yscale("log")
        cbar=plt.colorbar(sc,orientation='horizontal')
        plt.show()
    
    def S_Tr(self):
        import stockwell as smt
        import numpy as np
        print()
        dt = self.Z[0].stats.delta
        df = 1/(dt*self.Z[0].stats.npts)
        self.Nfmin = int(self.f1/df)
        self.Nfmax = int(self.f2/df)
        print(dt,df,self.Nfmin,self.Nfmax)
        self.S_Z = smt.st(self.Z[0].data,self.Nfmin,self.Nfmax)
        self.S_E = smt.st(self.E[0].data,self.Nfmin,self.Nfmax)
        self.S_N = smt.st(self.N[0].data,self.Nfmin,self.Nfmax)
        self.fsig = np.arange(self.Nfmin*df,self.Nfmax*df+df,df)
        self.tsig = np.arange(0,self.Z[0].stats.npts*dt,dt)
    
    def filter(self,dop_inter,):
        import numpy as np
        import stockwell as smt
        if self.S_Z.shape[0] == dop_inter.shape[0] and self.S_Z.shape[1] == dop_inter.shape[1]:    
            Sfilt_Z = np.zeros([self.S_Z.shape[0],self.S_Z.shape[1]],dtype=complex)
            Sfilt_E = np.zeros([self.S_Z.shape[0],self.S_Z.shape[1]],dtype=complex)
            Sfilt_N = np.zeros([self.S_Z.shape[0],self.S_Z.shape[1]],dtype=complex)
            for i in range(0,self.S_Z.shape[0]):
                for j in range(0,self.S_Z.shape[1]):
                    Sfilt_Z[i][j] = self.S_Z[i][j]*dop_inter[i][j]
                    Sfilt_E[i][j] = self.S_E[i][j]*dop_inter[i][j]
                    Sfilt_N[i][j] = self.S_N[i][j]*dop_inter[i][j]
        else:
            raise ValueError ("S_signal and dop_inter are not equal!!!!")
        Z_filtered = smt.ist(Sfilt_Z,self.Nfmin,self.Nfmax)
        E_filtered = smt.ist(Sfilt_E,self.Nfmin,self.Nfmax)
        N_filtered = smt.ist(Sfilt_N,self.Nfmin,self.Nfmax)
        self.Z_filt = self.Z.Copy()
        self.E_filt = self.E.Copy()
        self.N_filt = self.N.Copy()
        self.Z_filt[0].data = Z_filtered
        self.E_filt[0].data = E_filtered
        self.N_filt[0].data = N_filtered