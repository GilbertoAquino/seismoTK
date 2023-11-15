from .helpers.ASA2SAC_helpers import *
from seismoTK.Detection.Detection import *
from ..Detection.Detection import *
from tqdm import tqdm
import numpy as np
import os
from obspy import read

class RACM:
    def __init__(self,name=None,root=None,EarthquakeDate=None):
        self.name=name
        self.root=root
        self.date=EarthquakeDate
        self.__nlength = 45000
        try:
            self.read()
        except:
            self.V = None
            self.E = None
            self.N = None
            self.T = None
            self.R = None
    
    def read(self):
        rootsave='./'+str(self.name)+'/'
        os.chdir(rootsave)
        self.V=read('*.V.*')
        self.E=read('*.E.*')
        self.N=read('*.N.*')
        try:
            self.T = read('*.T.*')
            self.R = read('*.R.*')
        except:
            self.T = None
            self.R = None
        os.chdir("../")

    def write(self,CR=False):
        rootsave='./'+str(self.name)+'/'
        os.chdir(rootsave)
        for i in range(0,len(self.V)):
            self.V[i].write(self.V[i].stats.station+'.V.sac')
            self.E[i].write(self.V[i].stats.station+'.E.sac')
            self.N[i].write(self.V[i].stats.station+'.N.sac')
        if CR:
            self.R[i].write(self.V[i].stats.station+'.R.sac')
            self.T[i].write(self.V[i].stats.station+'.T.sac')
        os.chdir("../")

    def ASA2SAC(self):
        if (self.V != None and self.N != None and self.E != None):
            raise ValueError("Stream objects are not NULL (None)")
        name=self.name
        root=self.root
        rootsave='./'+str(name)+'/'
        FECHA=self.date
        os.mkdir(str(name))
        Estaciones = os.listdir(root)
        print("Converting ASA files to SAC:")
        for estacion in tqdm(Estaciones):
            aa=[]
            next_line_longitud = False
            next_line_longitud_sis = False
            start_data = False
            count = 0
            with open(os.path.join(root,estacion)) as est:
                for line in est.readlines():
                    if next_line_longitud:
                        Longitud = float(line.strip().split(":")[1].strip().split()[0])
                        next_line_longitud = False
                        continue
                    if next_line_longitud_sis:
                        lonsis = float(line.strip().split(":")[1].strip().split()[0])
                        next_line_longitud_sis = False
                        continue
                    if "HORA DE LA PRIMERA MUESTRA" in line:
                        HoraDeIn = line.strip().split(" : ")[1].strip()
                        horasis = float(HoraDeIn.split(":")[0])
                        minsis = float(HoraDeIn.split(":")[1])
                        segsis = float(HoraDeIn.split(":")[2])
                        continue
                    if "COORDENADAS DE LA ESTACION" in line:
                        Latitud = float(line.strip().split(" : ")[1].strip().split()[0])
                        next_line_longitud = True
                        continue
                    if "CLAVE DE LA ESTACION" in line:
                        Clave = line.strip().split(" : ")[1].strip()
                        continue
                    if "INSTITUCION RESPONSABLE" in line:
                        Institucion = line.strip().split(" : ")[1].strip()
                        continue
                    if "ORIENTACION C1-C6" in line:
                        Orientacion = [int for int in line.strip().split(" : ")[1].strip().split("/") if int != '']
                        continue
                    if "INTERVALO DE MUESTREO, C1" in line:
                        Intervalo = float([int for int in line.strip().split(" : ")[1].strip().split("/") if int != ''][0])
                        continue
                    if "COORDENADAS DEL EPICENTRO" in line:
                        latsis = float(line.strip().split(" : ")[1].strip().split()[0])
                        next_line_longitud_sis = True
                        continue
                    #if "HORA EPICENTRO" in line:
                    #    HoraSis = line.strip().split(" : ")[1].strip()
                    #    horasis = float(HoraSis.split(":")[0])
                    #    minsis = float(HoraSis.split(":")[1])
                    #    segsis = float(HoraSis.split(":")[2])
                    #    continue
                    if "---------+---------+---------+" in line:
                        start_data = True
                        count += 1
                        continue
                    if start_data:
                        if count == 2:
                            #print([data for data in line.strip().split(" ") if data != ""])
                            aa.append([data for data in line.strip().split(" ") if data != ""])
                if aa[len(aa)-1]==[]:
                    aa.pop(len(aa)-1)
                comp1=np.zeros(len(aa))
                comp2=np.zeros(len(aa))
                comp3=np.zeros(len(aa))
                lengthoffile=len(aa)
                for jj in range(0,lengthoffile):
                    comp=np.zeros(3)
                    example_list = aa[jj]
                    ii=0
                    for x1 in range(0,len(example_list)):
                        if example_list[x1]=='':
                            pass
                        else:
                            try:
                                comp[ii]=example_list[x1]
                                ii=ii+1
                            except:
                                if example_list[x1][0]=='-':
                                    do=example_list[x1][1:]
                                    example_list2 = [k for k in do.split('-')]
                                    comp[ii]=float(example_list2[0])*-1.0
                                    comp[ii+1]=float(example_list2[1])*-1.0
                                else:
                                    do=example_list[x1]
                                    example_list2 = [k for k in do.split('-')]
                                    comp[ii]=float(example_list2[0])
                                    comp[ii+1]=float(example_list2[1])*-1.0
                    comp1[jj]=comp[0]
                    comp2[jj]=comp[1]
                    comp3[jj]=comp[2]
                O1=Orientacion[0]
                O2=Orientacion[1]
                O3=Orientacion[2]
                O1,comp1=AsignacionDeOrientacion(O1,comp1)
                O2,comp2=AsignacionDeOrientacion(O2,comp2)
                O3,comp3=AsignacionDeOrientacion(O3,comp3)
                ASATOSAC(comp1,O1,Clave,Intervalo,Latitud,-Longitud,latsis,-lonsis,horasis,minsis,segsis,FECHA,rootsave)
                ASATOSAC(comp3,O3,Clave,Intervalo,Latitud,-Longitud,latsis,-lonsis,horasis,minsis,segsis,FECHA,rootsave)
                ASATOSAC(comp2,O2,Clave,Intervalo,Latitud,-Longitud,latsis,-lonsis,horasis,minsis,segsis,FECHA,rootsave)
        self.read()
        

    def CheckDelta(self):
        import os
        import statistics
        rootsave='./'+str(self.name)+'/'
        self.V.detrend(type="linear")
        self.E.detrend(type="linear")
        self.N.detrend(type="linear")
        delta=[]
        for i in range(0,len(self.V)):
            delta.append(self.V[i].stats.delta)
        maxdelta=max(delta)
        mode=statistics.mode(delta)
        print(maxdelta,mode)
        print("..Revisando Intervalos de muestreo..")
        for i in range(0,len(self.V)):
            deltaisok = False
            while deltaisok == False:
                if delta[i]<0.01:
                    print("Estación: "+str(self.V[i].stats.station)+" Delta: "+str(delta[i]))
                    self.V[i].decimate(2,strict_length=False, no_filter=True)
                    self.E[i].decimate(2,strict_length=False, no_filter=True)
                    self.N[i].decimate(2,strict_length=False, no_filter=True)
                    print("Nuevo delta: "+str(self.V[i].stats.delta)+" "+str(self.E[i].stats.delta)+" "+str(self.N[i].stats.delta))
                    delta[i] = self.V[i].stats.delta
                    self.write()
                elif delta[i]>0.01:
                    print("Estación: "+str(self.V[i].stats.station)+" Delta: "+str(delta[i]))
                    self.V[i].interpolate(sampling_rate=100)
                    self.E[i].interpolate(sampling_rate=100)
                    self.N[i].interpolate(sampling_rate=100)
                    print("Nuevo delta: "+str(self.V[i].stats.delta)+" "+str(self.E[i].stats.delta)+" "+str(self.N[i].stats.delta))
                    delta[i] = self.V[i].stats.delta
                    self.write()
                elif delta[i] == 0.01:
                    if self.V[i].stats.delta == self.E[i].stats.delta and self.V[i].stats.delta == self.N[i].stats.delta:
                        deltaisok = True
                    else:
                        print(self.V[i].stats.station,self.E[i].stats.station,self.N[i].stats.station)
                        raise ValueError('DELTA VALUES MISMATCH!')

    def plot_trigger1(self,trace, cft, thr_on, thr_off, df, lin=False, show=True):
        import matplotlib.pyplot as plt
        import numpy as np
        from obspy.signal.trigger import trigger_onset
        npts = len(trace)
        t = np.arange(npts, dtype=np.float32)
        fig = plt.figure()
        ax1 = fig.add_subplot(211)
        if lin==True:
            ax1.semilogy(t, trace.data, 'k')
        else:
            ax1.plot(t,trace.data,'k')    
        ax2 = fig.add_subplot(212, sharex=ax1)
        ax2.plot(t, cft, 'k')
        on_off = np.array(trigger_onset(cft, thr_on, thr_off))
        i, j = ax1.get_ylim()
        try:
            ax1.vlines(on_off[:, 0], i, j, color='r', lw=2,
                       label="Trigger On")
            ax1.vlines(on_off[:, 1], i, j, color='b', lw=2, label="Trigger Off")
            ax1.legend()
        except IndexError:
            pass
        ax2.axhline(thr_on, color='red', lw=1, ls='--')
        ax2.axhline(thr_off, color='blue', lw=1, ls='--')
        ax2.set_xlabel("Tiempo [muestras]")
        fig.suptitle("STA/LTA")
        fig.canvas.draw()
        return on_off

    def STALTA_for_10secods_pulse(self,V,show=True):
        from obspy.signal.trigger import recursive_sta_lta as stalta
        import matplotlib.pyplot as plt
        df = V.stats.sampling_rate
        cft = stalta(V.data,int(15 * df) ,int(45 * df))
        a=self.plot_trigger1(V, cft, 1.5, 0.9,df)
        if show:
            plt.show()
        else:
            plt.close()
        if a.shape[0] == 1:
            return a[0]
        elif a.shape[0] == 0:
            plt.show()
            return [0,0]
        else:
            print("Se selecciono el primer pulso")
            return a[0]

    def Alinear(self,stationname=[None]):
        name=self.name
        root=self.root
        rootsave='./'+str(name)+'/'
        FECHA=self.date
        os.chdir(rootsave)
        V=read('*.V.sac')
        E=read('*.E.sac')
        N=read('*.N.sac')
        try:
            R=read('*.R.sac')
            T=read('*.T.sac')
            Componentes_Rotadas = True
        except:
            Componentes_Rotadas = False
        V.detrend(type="linear")
        E.detrend(type="linear")
        N.detrend(type="linear")
        if Componentes_Rotadas:
            R.detrend(type="linear")
            T.detrend(type="linear")
        VF=V.copy()
        VF.filter("bandpass",corners=4, freqmin=0.0833, freqmax=0.125, zerophase=True)
        clave=[]
        dist=[]
        velfreq=[]
        j=None
        for i in range(0,len(V)):
            clave.append(V[i].stats.station)
            dist.append(V[i].stats.sac.dist)
            velfreq.append(V[i].stats.sampling_rate)
            #plt.plot(VF[i].data)
            if len(stationname) == 1:
                if V[i].stats.station == stationname[0]:
                    j=i
            else:
                if V[i].stats.station in stationname:
                    j=i
        if j==None:
            mindist=min(dist)
            mindistindex=dist.index(min(dist))
            print("Minima distancia: "+clave[mindistindex])

        else:
            mindist=dist[j]
            mindistindex=j
        time_min = V[mindistindex].stats.starttime
        time_max = V[mindistindex].stats.endtime
        #plt.show()
        pulse = [int(VF[mindistindex].stats.sac.a*100),int(VF[mindistindex].stats.sac.t1*100)]
        pulsemenos = pulse[0]
        pulsemas = pulse[1]
        data2correlatemenos=[]
        data2correlatemas=[]
        for i in range(0,len(V)):
            data2correlatemenos.append(int(V[i].stats.sac.a*100))
            data2correlatemas.append(int(V[i].stats.sac.t1*100))
            #pp = self.STALTA_for_10secods_pulse(VF[i],show=False)
            #if pp[1] == 0:
            #    print('STA/LTA no encontro un pulso')
                #data2correlatemenos.append(0)
            #    data2correlatemas.append(len(VF[i].data))
            #else:
                #data2correlatemenos.append(int(pp[0])-500)
            #    data2correlatemas.append(int(pp[1])+250)
        d2menos =  0 #int(min(data2correlatemenos))
        #data2correlatemas=28000
        #print(data2correlatemenos,data2correlatemas)
        cor=np.correlate(VF[mindistindex].data[d2menos:data2correlatemas[mindistindex]],VF[mindistindex].data[pulsemenos:pulsemas],mode='full')
        #cor=np.correlate(VF[mindistindex].data,VF[mindistindex].data[pulsemenos:pulsemas],mode='full')
        cor=cor/max(cor)
        maxargument=np.argmax(cor)

        print('Minima distancia: ',mindist,'Estación: ',clave[mindistindex])
        times = []
        if Componentes_Rotadas:
            for i in range(0,len(VF)):
                #a=np.correlate(VF[i].data[data2correlatemenos:data2correlatemas],VF[mindistindex].data[pulsemenos:pulsemas],mode='full')
                a=np.correlate(VF[i].data[d2menos:data2correlatemas[i]],VF[mindistindex].data[pulsemenos:pulsemas],mode='full')
                if np.abs(max(a)) < np.abs(min(a)):
                    #print("Estacion: ",clave[i]," Componente vertical invertida. Corrigiendo...")
                    VF[i].data= VF[i].data*-1
                    V[i].data = V[i].data*-1
                    a=a*-1
                a=a/max(a)
                mindex1=np.argmax(a)
                print(clave[i],mindex1,maxargument)
                while mindex1 < maxargument:
                    a=np.insert(a,0,0)
                    VF[i].data=np.insert(VF[i].data,0,0)
                    V[i].data=np.insert(V[i].data,0,0)
                    E[i].data=np.insert(E[i].data,0,0)
                    N[i].data=np.insert(N[i].data,0,0)
                    T[i].data=np.insert(T[i].data,0,0)
                    R[i].data=np.insert(R[i].data,0,0)
                    mindex1=np.argmax(a)

                while mindex1 > maxargument:
                    a=np.delete(a,0)
                    VF[i].data=np.delete(VF[i].data,0)
                    V[i].data=np.delete(V[i].data,0)
                    E[i].data=np.delete(E[i].data,0)
                    N[i].data=np.delete(N[i].data,0)
                    R[i].data=np.delete(E[i].data,0)
                    T[i].data=np.delete(N[i].data,0)
                    mindex1=np.argmax(a)
        else:
            for i in range(0,len(VF)):
                time = time_min
                #a=np.correlate(VF[i].data[data2correlatemenos:data2correlatemas],VF[mindistindex].data[pulsemenos:pulsemas],mode='full')
                a=np.correlate(VF[i].data[d2menos:data2correlatemas[i]],VF[mindistindex].data[pulsemenos:pulsemas],mode='full')
                if np.abs(max(a)) < np.abs(min(a)):
                    #print("Estacion: ",clave[i]," Componente vertical invertida. Corrigiendo...")
                    VF[i].data= VF[i].data*-1
                    V[i].data = V[i].data*-1
                    a=a*-1
                a=a/max(a)
                mindex1=np.argmax(a)
                print(clave[i],mindex1,maxargument)
                while mindex1 < maxargument:
                    a=np.insert(a,0,0)
                    VF[i].data=np.insert(VF[i].data,0,0)
                    V[i].data=np.insert(V[i].data,0,0)
                    E[i].data=np.insert(E[i].data,0,0)
                    N[i].data=np.insert(N[i].data,0,0)
                    mindex1=np.argmax(a)
                    #time += 0.01

                while mindex1 > maxargument:
                    a=np.delete(a,0)
                    VF[i].data=np.delete(VF[i].data,0)
                    V[i].data=np.delete(V[i].data,0)
                    E[i].data=np.delete(E[i].data,0)
                    N[i].data=np.delete(N[i].data,0)
                    mindex1=np.argmax(a)
                    #time -= 0.01
                times.append(time)

        rmenos = pulse[0] #int(input("-->Limite inferior de la fase a alinear: "))
        rmas = pulse[1] #input("-->Limite superior de la fase a alinear: "))
        maximo=np.argmax(VF[mindistindex].data[rmenos:rmas])
        rmenos = maximo + pulse[0] - 100
        rmas = maximo + pulse[0] + 100
        maximo=np.argmax(VF[mindistindex].data[rmenos:rmas])
        print('Alineacion de maximos')
        if Componentes_Rotadas:
            for i in range(0,len(V)):
                maximoi = np.argmax(VF[i].data[rmenos:rmas])
                print(clave[i],maximoi,maximo)
                if maximoi == 0:
                    maximoi = maximo
                while maximoi < maximo:
                    VF[i].data=np.insert(VF[i].data,0,0)
                    V[i].data=np.insert(V[i].data,0,0)
                    E[i].data=np.insert(E[i].data,0,0)
                    N[i].data=np.insert(N[i].data,0,0)
                    T[i].data=np.insert(T[i].data,0,0)
                    R[i].data=np.insert(R[i].data,0,0)
                    maximoi = np.argmax(VF[i].data[rmenos:rmas])
                while maximoi > maximo:
                    VF[i].data=np.delete(VF[i].data,0)
                    V[i].data=np.delete(V[i].data,0)
                    E[i].data=np.delete(E[i].data,0)
                    N[i].data=np.delete(N[i].data,0)
                    T[i].data=np.delete(T[i].data,0)
                    R[i].data=np.delete(R[i].data,0)
                    maximoi = np.argmax(VF[i].data[rmenos:rmas])
            for i in range(0,len(V)):
                V[i].write(str(clave[i])+'.'+'V.sac')
                E[i].write(str(clave[i])+'.'+'E.sac')
                N[i].write(str(clave[i])+'.'+'N.sac')
                T[i].write(str(clave[i])+'.'+'T.sac')
                R[i].write(str(clave[i])+'.'+'R.sac')
        else:
            for i in range(0,len(V)):
                time = times[i]
                maximoi = np.argmax(VF[i].data[rmenos:rmas])
                print(clave[i],maximoi,maximo)
                if maximoi == 0:
                    maximoi = maximo
                while maximoi < maximo:
                    VF[i].data=np.insert(VF[i].data,0,0)
                    V[i].data=np.insert(V[i].data,0,0)
                    E[i].data=np.insert(E[i].data,0,0)
                    N[i].data=np.insert(N[i].data,0,0)
                    maximoi = np.argmax(VF[i].data[rmenos:rmas])
                    #time += 0.01
                while maximoi > maximo:
                    VF[i].data=np.delete(VF[i].data,0)
                    V[i].data=np.delete(V[i].data,0)
                    E[i].data=np.delete(E[i].data,0)
                    N[i].data=np.delete(N[i].data,0)
                    maximoi = np.argmax(VF[i].data[rmenos:rmas])
                    #time -= 0.01
                times[i] = time
            for i in range(0,len(V)):
                os.remove(str(clave[i])+'.V.sac')
                os.remove(str(clave[i])+'.E.sac')
                os.remove(str(clave[i])+'.N.sac')
                V[i].stats.sac["nzhour"] = times[i].hour
                V[i].stats.sac["nzmin"] = times[i].minute
                V[i].stats.sac["nzsec"] = times[i].second
                V[i].stats.sac["nzmsec"] = times[i].microsecond
                N[i].stats.sac["nzhour"] = times[i].hour
                N[i].stats.sac["nzmin"] = times[i].minute
                N[i].stats.sac["nzsec"] = times[i].second
                N[i].stats.sac["nzmsec"] = times[i].microsecond
                E[i].stats.sac["nzhour"] = times[i].hour
                E[i].stats.sac["nzmin"] = times[i].minute
                E[i].stats.sac["nzsec"] = times[i].second
                E[i].stats.sac["nzmsec"] = times[i].microsecond
                V[i].stats.starttime = times[i]
                E[i].stats.starttime = times[i]
                N[i].stats.starttime = times[i]
                print("Guardadon datos!")
                print(V[i].stats)
                V[i].write(str(clave[i])+'.'+'V.2.sac', format='SAC')
                E[i].write(str(clave[i])+'.'+'E.2.sac', format='SAC')
                N[i].write(str(clave[i])+'.'+'N.2.sac', format='SAC')
        os.chdir('../')

    def CheckAlinear(self):
        import matplotlib.pyplot as plt
        import numpy as np
        import os
        from obspy import read
        import statistics
        name=self.name
        os.chdir(name)
        V=read('*.V.sac')
        V.filter("bandpass",corners=4, freqmin=0.0833, freqmax=0.125, zerophase=True)
        beginTrigger=[]
        endTrigger=[]
        s=False
        for i in range(0,len(V)):
        	beginTrigger.append(int(V[i].stats.sac.a*100))
        	endTrigger.append(int(V[i].stats.sac.t1*100))
        """    if V[i].stats.station in [None]:
                a = self.STALTA_for_10secods_pulse(V[i],show=True)
            else:
                a = self.STALTA_for_10secods_pulse(V[i],show=s)
            beginTrigger.append(a[0])
           endTrigger.append(a[1])
        """
        bmean=statistics.mean(beginTrigger)
        bmode=statistics.mode(beginTrigger)
        bmedian=statistics.median(beginTrigger)
        bstdev=statistics.stdev(beginTrigger)
        emean=statistics.mean(endTrigger)
        emode=statistics.mode(endTrigger)
        emedian=statistics.median(endTrigger)
        estdev=statistics.stdev(endTrigger)
        print('Imprimiendo resultados!!!!')
        print("Promedio del inicio: ",bmean,"Moda del inicio: ",bmode)
        print("Media del inicio: ",bmedian,"Desviacion Estandar del inicio: ",bstdev)
        print('#####################################################################')
        print("Promedio del final: ",emean,"Moda del final: ",emode)
        print("Media del final: ",emedian,"Desviacion Estandar del final: ",estdev)
        """
        badT=[]
        for i in range(0,len(V)):
            if endTrigger[i] > emedian + estdev or endTrigger[i] < emedian - estdev:
                badT.append(i)
            elif V[i].stats.station in ['CUP5']:
                badT.append(i)
        plottraces(V,badT,beginTrigger,endTrigger)
        """
        Pos = np.arange(0,len(V),1)
        plottraces(V,beginTrigger,endTrigger,Pos)
        for i in range(0,len(V)):
            plt.plot(V[i],color='black')
        plt.show()
        os.chdir('../')

    def Velocidad(self):
        import numpy as np
        import os
        from obspy import read
        import statistics
        name=self.name
        root=self.root
        rootsave='./'+str(name)+'/'
        FECHA=self.date
        os.chdir(rootsave)
        V=read('*.V.sac')
        E=read('*.E.sac')
        N=read('*.N.sac')
        try:
            R=read('*.R.sac')
            T=read('*.T.sac')
            Componentes_Rotadas = True
        except:
            Componentes_Rotadas = False
        VF=V.copy()
        EF=E.copy()
        NF=N.copy()
        if Componentes_Rotadas:
            TF=T.copy()
            RF=R.copy()
        clave=[]
        dist=[]
        velfreq=[]
        claveOrden=[]
        distOrden=[]
        velfreqOrden=[]
        l=len(V)
        if Componentes_Rotadas:
            for i in range(0,len(V)):
                clave.append(V[i].stats.station)
                dist.append(V[i].stats.sac.dist)
                velfreq.append(V[i].stats.sampling_rate)
                a1=VF.select(station=str(clave[i]))
                a2=EF.select(station=str(clave[i]))
                a3=NF.select(station=str(clave[i]))
                a4=RF.select(station=str(clave[i]))
                a5=TF.select(station=str(clave[i]))
                VF.remove(a1[0])
                EF.remove(a2[0])
                NF.remove(a3[0])
                RF.remove(a4[0])
                TF.remove(a5[0])
        else:
            for i in range(0,len(V)):
                clave.append(V[i].stats.station)
                dist.append(V[i].stats.sac.dist)
                velfreq.append(V[i].stats.sampling_rate)
                a1=VF.select(station=str(clave[i]))
                a2=EF.select(station=str(clave[i]))
                a3=NF.select(station=str(clave[i]))
                VF.remove(a1[0])
                EF.remove(a2[0])
                NF.remove(a3[0])
        distmean = statistics.mean(dist)
        if Componentes_Rotadas:
            for i in range(0,len(V)):
                mindist=min(dist)
                mindistindex=dist.index(min(dist))
                claveOrden.append(clave[mindistindex])
                distOrden.append(mindist)
                velfreqOrden.append(velfreq[mindistindex])
                clave.pop(mindistindex)
                dist.pop(mindistindex)
                velfreq.pop(mindistindex)
                VF.append(V[mindistindex])
                EF.append(E[mindistindex])
                NF.append(N[mindistindex])
                RF.append(R[mindistindex])
                TF.append(T[mindistindex])
                a1=V.select(station=str(claveOrden[i]))
                a2=E.select(station=str(claveOrden[i]))
                a3=N.select(station=str(claveOrden[i]))
                a4=R.select(station=str(claveOrden[i]))
                a5=T.select(station=str(claveOrden[i]))
                V.remove(a1[0])
                E.remove(a2[0])
                N.remove(a3[0])
                R.remove(a4[0])
                T.remove(a5[0])
            V=VF.copy()
            E=EF.copy()
            N=NF.copy()
            R=RF.copy()
            T=TF.copy()
        else:
            for i in range(0,len(V)):
                mindist=min(dist)
                mindistindex=dist.index(min(dist))
                claveOrden.append(clave[mindistindex])
                distOrden.append(mindist)
                velfreqOrden.append(velfreq[mindistindex])
                clave.pop(mindistindex)
                dist.pop(mindistindex)
                velfreq.pop(mindistindex)
                VF.append(V[mindistindex])
                EF.append(E[mindistindex])
                NF.append(N[mindistindex])
                a1=V.select(station=str(claveOrden[i]))
                a2=E.select(station=str(claveOrden[i]))
                a3=N.select(station=str(claveOrden[i]))
                V.remove(a1[0])
                E.remove(a2[0])
                N.remove(a3[0])
            V=VF.copy()
            E=EF.copy()
            N=NF.copy()

        vel=3.2
        timemean = distmean/vel
        pp=round(distOrden[0]/vel,2)
        if Componentes_Rotadas:
            for i in range(1,l):
                p=round(distOrden[i]/vel,2)
                Nmas=int((p-pp)*velfreqOrden[i])
                for j in range(0,Nmas):
                    V[i].data=np.insert(V[i].data,0,0)
                    E[i].data=np.insert(E[i].data,0,0)
                    N[i].data=np.insert(N[i].data,0,0)
                    R[i].data=np.insert(R[i].data,0,0)
                    T[i].data=np.insert(T[i].data,0,0)
        else:
            for i in range(1,l):
                p=round(distOrden[i]/vel,2)
                Nmas=int((p-pp)*velfreqOrden[i])
                for j in range(0,Nmas):
                    V[i].data=np.insert(V[i].data,0,0)
                    E[i].data=np.insert(E[i].data,0,0)
                    N[i].data=np.insert(N[i].data,0,0)
            #print(Nmas)
        Numpts=[]
        for i in range(0,l):
            Numpts.append(V[i].stats.npts)
        Nmax=int(statistics.mean(Numpts))
        print(Nmax)
        #print(Nmax)
        if Componentes_Rotadas:
            for i in range(0,l):
                while(len(V[i].data)<Nmax):
                    V[i].data=np.insert(V[i].data,len(V[i].data),0)
                    E[i].data=np.insert(E[i].data,len(E[i].data),0)
                    N[i].data=np.insert(N[i].data,len(N[i].data),0)
                    R[i].data=np.insert(R[i].data,len(R[i].data),0)
                    T[i].data=np.insert(T[i].data,len(T[i].data),0)
                while (len(V[i].data)>Nmax):
                    V[i].data=np.delete(V[i].data,len(V[i].data)-1)
                    E[i].data=np.delete(E[i].data,len(E[i].data)-1)
                    N[i].data=np.delete(N[i].data,len(N[i].data)-1)
                    R[i].data=np.delete(R[i].data,len(R[i].data)-1)
                    T[i].data=np.delete(T[i].data,len(T[i].data)-1)
                print(i+1,"Estación: "+str(claveOrden[i])+' Ajustada con: ',len(V[i].data))
        else:
            for i in range(0,l):
                while(len(V[i].data)<Nmax):
                    V[i].data=np.insert(V[i].data,len(V[i].data),0)
                    E[i].data=np.insert(E[i].data,len(E[i].data),0)
                    N[i].data=np.insert(N[i].data,len(N[i].data),0)
                while (len(V[i].data)>Nmax):
                    V[i].data=np.delete(V[i].data,len(V[i].data)-1)
                    E[i].data=np.delete(E[i].data,len(E[i].data)-1)
                    N[i].data=np.delete(N[i].data,len(N[i].data)-1)
                print(i+1,"Estación: "+str(claveOrden[i])+' Ajustada con: ',len(V[i].data))
                
            '''if len(V[i].data)==Nmax:
                V[i].stats.starttime = V[i].stats.starttime - timemean
                E[i].stats.starttime = E[i].stats.starttime - timemean
                N[i].stats.starttime = N[i].stats.starttime - timemean
                if V[i].stats.sac.nzsec < timemean:
                    V[i].stats.sac.nzsec = V[i].stats.sac.nzsec - timemean + 60
                    V[i].stats.sac.nzmin = V[i].stats.sac.nzmin - 1
                    E[i].stats.sac.nzsec = E[i].stats.sac.nzsec - timemean + 60
                    E[i].stats.sac.nzmin = E[i].stats.sac.nzmin - 1
                    N[i].stats.sac.nzsec = N[i].stats.sac.nzsec - timemean + 60
                    N[i].stats.sac.nzmin = N[i].stats.sac.nzmin - 1
                else:
                    V[i].stats.sac.nzsec = V[i].stats.sac.nzsec - timemean
                    E[i].stats.sac.nzsec = E[i].stats.sac.nzsec - timemean
                    N[i].stats.sac.nzsec = N[i].stats.sac.nzsec - timemean'''   
            

        for i in range(0,l):
            V[i].write(str(claveOrden[i])+'.'+'V.sac')
            E[i].write(str(claveOrden[i])+'.'+'E.sac')
            N[i].write(str(claveOrden[i])+'.'+'N.sac')
            try:
                R[i].write(str(claveOrden[i])+'.'+'R.sac')
                T[i].write(str(claveOrden[i])+'.'+'T.sac')
            except:
                pass
        os.chdir('../')

    def back(self):
        import os
        os.chdir('..')

    def rotar_sac(self):
        import os
        os.chdir(self.name)
        os.system("ls *.V.sac | awk -F[.] '{print $1}' > stations.txt")
        stations = self.readcsv('stations.txt')
        os.system("rm stations.txt")
        self.back()
        os.system("echo 'message station-$1$;r ./$2$/$1$.E.sac;ch CMPINC 90;ch CMPAZ 90;wh;r ./$2$/$1$.N.sac;ch CMPINC 90;ch CMPAZ 0;wh;r ./$2$/$1$.N.sac ./$2$/$1$.E.sac;rotate to GCP;w ./$2$/$1$.R.sac ./$2$/$1$.T.sac;message Rotate-Complete;q' > rotate.m")
        for i in stations:  
            os.system('sac rotate.m '+str(i[0])+' '+self.name)
        os.system('rm rotate.m')

    def readcsv(self, name):
        import csv
        with open(name) as cvsfile:
            file=csv.reader(cvsfile,delimiter=' ')
            array=[]
            for i in file:
                array.append(i)
        return array

    def fit_pulses(self,show=False):
        from obspy import read
        from tensorflow import keras
        import numpy as np
        from pathlib import Path
        BASE_DIR = Path(__file__).resolve().parent.parent
        if (self.R == None or self.T == None):
            Componentes_Rotadas = False
        else:
            Componentes_Rotadas = True
        clave = []
        for i in self.V:
            clave.append(i.stats.station)
        ReduceTime(self.V,N=self.__nlength)
        ReduceTime(self.N,N=self.__nlength)
        ReduceTime(self.E,N=self.__nlength)
        if Componentes_Rotadas:
            ReduceTime(self.R,N=self.__nlength)
            ReduceTime(self.T,N=self.__nlength)
        VF=self.V.copy()
        VF.filter("bandpass",corners=4, freqmin=0.0833, freqmax=0.125, zerophase=True)
        X = spectrograms_pulses(VF)
        X = np.array(X)
        model = keras.models.load_model(BASE_DIR / 'ANN_Model_pulses2.tf')
        y_pred = model.predict(X)
        yb = []
        ye = []
        j=0
        for i in y_pred:
            if i[0] < 0:
            	i[0] = 0
            yb.append(i[0]*self.__nlength)
            ye.append(i[1]*self.__nlength)
            self.V[j].stats.sac['a'] = i[0]*450
            self.V[j].stats.sac['t1'] = i[1]*450
            j=j+1
        Pos = np.arange(0,len(self.V),1)
        if show:
            plottraces(self.V,yb,ye,Pos)
        self.write(CR=Componentes_Rotadas)
    
    def get_picks(self,stream,normalized = True):
        import numpy as np
        yb = []
        ye = []
        y = []
        dum = np.zeros([len(stream),2])
        for i in stream:
            try:
                yb.append(i.stats.sac.a)
            except:
                print(i.stats.sac)
                raise ValueError('Couldnt find AMarker in sac file: '+i.stats.station)
            try:
                ye.append(i.stats.sac.t1)
            except:
                raise ValueError('Couldnt find T1Marker in sac file: '+i.stats.station)
        for i in range(0,len(yb)):
            if normalized:
                yb[i]=yb[i]/(self.__nlength/100)
                ye[i]=ye[i]/(self.__nlength/100)
            dum[i][0] = yb[i]
            dum[i][1] = ye[i]
            y.append(dum[i])
        return y

    def CheckTraces(self):
        rootsave='./'+str(self.name)+'/'
        print('Componentes Verticales leidas: ',len(self.V))
        print('Componentes Norte leidas: ',len(self.N))
        print('Componentes Este leidas: ',len(self.E))
        if len(self.V) == len(self.N) and len(self.V) == len(self.E):
            print('TODO EN ORDEN')
        else:
            print('REVISAR TRAZAS')

    def readcsv(self,name,delimit=' '):
        import csv
        with open(name) as cvsfile:
            file=csv.reader(cvsfile,delimiter=delimit)
            array=[]
            for i in file:
                array.append(i)
        return array

    def makeall(self):
        import os
        import pandas as pd
        os.system('mkdir SAC')
        os.system('mkdir ASA')
        os.system('ls > datos.txt')
        leidos = self.readcsv('datos.txt')
        os.system('rm datos.txt')
        sismos=[]
        j=0
        for i in range(0,len(leidos)):
            if leidos[i][0] in ['datos.txt','NoOrganizados','SAC','ASA']:
                pass
            else:
                sismos.append(leidos[i][0])
        
        for j in range(0,len(sismos)):
            if sismos[j][len(sismos[j])-3] + sismos[j][len(sismos[j])-2] + sismos[j][len(sismos[j])-1]  == 'sac':
                sismos.pop(j)
        print(sismos)    
        print('SISMOS LEIDOS!!! ----> ',len(sismos),' carpetas con archivos formato ASA!!!!')
        df = self.readcsv('../Control.csv',delimit=',')
        for sismo in sismos:
            self.name = sismo+'-sac'
            self.root = sismo
            self.date = sismo
            print('INICIO ASA2SAC: ',sismo)
            self.ASA2SAC()
            print('Finalizo con exito')
            print('INICIO CHECKDELTA: ',sismo)
            self.CheckDelta()
            print('Finalizo con exito')
            print('INICIO AJUSTE DE FASES: ',sismo)
            self.fit_pulses()
            print('Finalizo con exito')
            print('INICIO ALINEAR: ',sismo)
            self.Alinear(['CUP5','CUP4'])
            print('Finalizo con exito')
            print('INICIO ASIGNACION DE VELOCIDAD: ',sismo)
            self.Velocidad()
            print('Finalizo con exito')
            print('INICIO AJUSTE DE FASES: ',sismo)
            self.fit_pulses()
            print('Finalizo con exito')
            try:
                df = self.readcsv('../Control.csv',delimit=',')
                df.append([self.root,self.name,self.root,self.date])
                self.writecvs(df,'../Control.csv')
            except:
                pass
            os.system('mv '+self.name+' SAC')
            os.system('mv '+self.name+'Est.txt'+' SAC')
            os.system('mv '+self.root+' ASA')

    def writecvs(self,Aimprimir,name):
        import csv
        with open(str(name), mode='w') as p:
            p2 = csv.writer(p)
            p2.writerows(Aimprimir)