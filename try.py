from seismoTK import Polarization
from seismoTK import S_Filter
"""pl = Polarization(f1 = 0.1,f2 = 5,nfr = 200,nsp = 3600,cycle = 3,wlenf=25)
pl.read_station('../Registros_Sismicos_RACM_2010-2020/Formato_SAC/19sep2017/','LV17')
pl.preprocessing(remove_fp=5000)
pl.get_pol()
pl.filter_data('LIN',0.5)
pl.plot_tz(pl.Pol["TIME"], pl.Pol["FREQ"], pl.Pol["LIN"],ranges=[0,1],s=5,dz=0.1,colormap='turbo')
pl.plot_tBazDop(pl.Pol["TIME"], pl.Pol["FREQ"], pl.Pol["BAZ"],pl.Pol["DOP"])
pl.freq_baz(yn=10,xn=10)
pl.save_pol('hola.csv')"""
pl2 = S_Filter(f1 = 0.1,f2 = 5,nfr = 100,nsp = 3600,cycle = 3,wlenf=25)
pl2.read_station('../Registros_Sismicos_RACM_2010-2020/Formato_SAC/19sep2017/','LV17')
pl2.preprocessing(remove_fp=5000)
pl2.get_pol()
pl2.plot_tBazDop(pl2.Pol["TIME"], pl2.Pol["FREQ"], pl2.Pol["BAZ"],pl2.Pol["DOP"])
pl2.S()
