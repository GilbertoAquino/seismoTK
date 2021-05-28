from seismoTK import Polarization
pl = Polarization()
pl.f1 = 0.1
pl.f2 = 5
pl.nfr = 200
pl.nsp = 3600
pl.cycle = 3
pl.wlenf=25
pl.read_station('../Registros_Sismicos_RACM_2010-2020/Formato_SAC/19sep2017/','LV17')
pl.preprocessing(remove_fp=5000)
pl.get_pol()
pl.filter_data('LIN',0.5)
pl.plot_tz(pl.Pol["TIME"], pl.Pol["FREQ"], pl.Pol["LIN"],ranges=[0,1],s=5,dz=0.1,colormap='turbo')
pl.plot_tBazDop(pl.Pol["TIME"], pl.Pol["FREQ"], pl.Pol["BAZ"],pl.Pol["DOP"])
pl.freq_baz()