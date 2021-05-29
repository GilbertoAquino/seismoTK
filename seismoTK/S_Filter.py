from matplotlib.colors import Colormap
from . import Polarization
class S_Filter(Polarization):
    def S(self):
        #import numpy as np
        import matplotlib.pyplot as plt
        from scipy.interpolate import griddata
        self.SF = self.Pol.drop(columns=["LIN","BAZ"])
        x,y,z=self.xyz2grd(self.SF["TIME"],self.SF["FREQ"],self.SF["DOP"],xn=3600,yn=100)
        

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
        plt.imshow(psx, cmap=cmap, extent=extent, aspect='auto', origin='lower')
        plt.yscale("log")
        cbar=plt.colorbar(sc,orientation='horizontal')
        plt.show()