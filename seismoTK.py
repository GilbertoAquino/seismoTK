from tkinter import ANCHOR
from unittest import skip

from numpy import pad


class RACM:
    def __init__(self,name=None,root=None,EarthquakeDate=None):
        self.name=name
        self.root=root
        self.date=EarthquakeDate
        self.nlength = 45000

    def ASA2SAC(self):
        import csv
        import numpy as np
        import os
        name=self.name
        root=self.root
        rootsave='./'+str(name)+'/'
        FECHA=self.date
        os.system("mkdir ./"+str(name))
        dirf='dirf *'
        cd='cd '+root
        os.chdir(root)
        os.system('rm ../Estacionesleidas.txt ../hpm.txt ../SDI.txt ../N.txt ../No.txt ../La.txt ../Lat.txt ../Lo.txt ../Lon.txt ../In.txt ../Insti.txt ../Ori.txt ../Inter.txt')
        os.system('rm ../Orientacion.txt ../Muestreo.txt')
        os.system("rm ../"+str(name)+"Est.txt")
        os.system("ls >> ../Estacionesleidas.txt")
        #os.os.system("awk '{print $3}' filenr.lis > ../Estacionesleidas.txt")
        Estaciones=[]
        with open('../Estacionesleidas.txt') as cvsfile:
            file=csv.reader(cvsfile,delimiter=' ')
            for i in file:
                Estaciones.append(i)
        #---------------------------------------------------------------------------#
        #---------------------------------------------------------------------------#
        #------------Tiempo de inicio-------------------#
        #---------------------------------------------------------------------------#
        #---------------------------------------------------------------------------#
        os.system('touch ../hpm.txt')
        for i in range(0,len(Estaciones)-1):
            a=str(Estaciones[i][0])
            #print(a)
            grep='grep "HORA DE LA PRIMER" '+str(a)+' >> ../hpm.txt'
            os.system(grep)
        os.system("awk '{print $8}' ../hpm.txt > ../SDI.txt")
        HoraDeIn=[]
        with open('../SDI.txt') as cvsfile:
            file=csv.reader(cvsfile,delimiter=':')
            for i in file:
                HoraDeIn.append(i)
        #segdein=[]
        #mindein=[]
        #for i in range(0,len(HoraDeIn)):
        #   mindein.append(float(HoraDeIn[i][1]))
        #   segdein.append(float(HoraDeIn[i][2]))
        #---------------------------------------------------------------------------#
        #---------------------------------------------------------------------------#
        #-----Nombre-LATITUD-LONGITUD e Institución-----#
        #---------------------------------------------------------------------------#
        #---------------------------------------------------------------------------#
        os.system('touch ../N.txt')
        os.system('touch ../No.txt')
        os.system('touch ../La.txt')
        os.system('touch ../Lat.txt')
        os.system('touch ../Lo.txt')
        os.system('touch ../Lon.txt')
        os.system('touch ../In.txt')
        os.system('touch ../Insti.txt')
        for i in range(0,len(Estaciones)):
            a=str(Estaciones[i][0])
            grep1='grep "COORDENADAS DE LA ESTACION" '+str(a)+' >> ../La.txt'
            grep2="awk 'NR<30' "+str(a)+" | awk '/LONG/' >> ../Lo.txt"
            grep3='grep "CLAVE DE LA ESTACION" '+str(a)+' >> ../N.txt'
            grep4='grep "INSTITUCION RESPONSABLE" '+str(a)+' >> ../In.txt'
            os.system(grep1)
            os.system(grep2)
            os.system(grep3)
            os.system(grep4)
        os.system("awk '{print $6}' ../La.txt > ../Lat.txt")
        os.system("awk '{print $2}' ../Lo.txt > ../Lon.txt")
        os.system("awk '{print $6}' ../N.txt > ../No.txt")
        os.system("awk '{print $4}' ../In.txt > ../Insti.txt")
        Latitud=[]
        Longitud=[]
        Clave=[]
        Institucion=[]
        with open('../Lat.txt') as cvsfile:
            file=csv.reader(cvsfile,delimiter=' ')
            for i in file:
                Latitud.append(i)
        with open('../Lon.txt') as cvsfile:
            file=csv.reader(cvsfile,delimiter=' ')
            for i in file:
                Longitud.append(i)
        with open('../No.txt') as cvsfile:
            file=csv.reader(cvsfile,delimiter=' ')
            for i in file:
                Clave.append(i)
        with open('../Insti.txt') as cvsfile:
            file=csv.reader(cvsfile,delimiter=' ')
            for i in file:
                Institucion.append(i)
        for i in range(0,len(Latitud)):
            Latitud[i]=float(Latitud[i][0])
            Longitud[i]=float(Longitud[i][0])
            Clave[i]=str(Clave[i][0])
            Institucion[i]=str(Institucion[i][0])
            #print(Clave[i],Institucion[i])
        with open(str(name)+'Est.txt', mode='w') as p:
            p2 = csv.writer(p,delimiter=' ')
            p2.writerow(['# Clave Institucion Latitud Longitud'])
            for i in range(0,len(Latitud)):
                p2.writerow([i+1,Clave[i],Institucion[i],Latitud[i], -Longitud[i]])
        os.system("mv "+str(name)+"Est.txt ../")
        #---------------------------------------------------------------------------#
        #---------------------------------------------------------------------------#
        #-------------------ORIENTACIÓN, DURACION Y MUESTREO------------------------#
        #---------------------------------------------------------------------------#
        #---------------------------------------------------------------------------#
        os.system('touch ../Ori.txt')
        os.system('touch ../Inter.txt')
        for i in range(0,len(Estaciones)):
            a=str(Estaciones[i][0])
            grep1='grep "INTERVALO DE MUESTREO, C1" '+str(a)+' >> ../Inter.txt'
            grep2='grep "ORIENTACION C1-C6" '+str(a)+' >> ../Ori.txt'
            os.system(grep1)
            os.system(grep2)
        os.system("awk -F/ '{print $3}' ../Inter.txt >> ../Muestreo.txt")
        os.system("awk -F/ '{print $2,$3,$4}' ../Ori.txt >> ../Orientacion.txt")
        Intervalo=[]
        Orientacion=[]
        VelInt=[]
        with open('../Muestreo.txt') as cvsfile:
            file=csv.reader(cvsfile,delimiter=' ')
            for i in file:
                Intervalo.append(i)
        with open('../Orientacion.txt') as cvsfile:
            file=csv.reader(cvsfile,delimiter=' ')
            for i in file:
                Orientacion.append(i)
        for i in range(0,len(Intervalo)):
            Intervalo[i]=float(Intervalo[i][0])
            VelInt.append(1.0/Intervalo[i])
        #-----------------------------------------------------------------------------#
        #-----------------------Caracteristicas del sismos-----------------------------#
        #-----------------------------------------------------------------------------#
        a=str(Estaciones[(len(Estaciones))//2][0])
        grep1='grep "COORDENADAS DEL EPICENTRO" '+str(a)+' >> ../SISMOLA.txt'
        grep2="awk 'NR>50' "+str(a)+" | awk '/LON/' >> ../SISMOLo.txt"
        grep3='grep "HORA EPICENTRO" '+str(a)+' >> ../SISMOHORA.txt'
        os.system(grep1)
        os.system(grep2)
        os.system(grep3)
        os.system("awk -F: '{print $2,$3,$4}' ../SISMOHORA.txt >> ../SISMO.txt")
        os.system("awk '{print $5}' ../SISMOLA.txt >> ../SISMO.txt")
        os.system("awk '{print $2}' ../SISMOLo.txt >> ../SISMO.txt")
        sismo=[]
        with open('../SISMO.txt') as cvsfile:
            file=csv.reader(cvsfile,delimiter=' ')
            for i in file:
                sismo.append(i)
        horasis=float(sismo[0][1])
        minsis=float(sismo[0][2])
        segsis=float(sismo[0][3][0:1])
        latsis=float(sismo[1][0])
        lonsis=float(sismo[2][0])
        os.system("rm ../SISMOLA.txt ../SISMOLo.txt ../SISMOHORA.txt ../SISMO.txt")
        #-----------------------------------------------------------------------------#
        #-----------------------------------------------------------------------------#
        #-----------------------------------------------------------------------------#
        for i in range(0,len(Estaciones)):
            print("---------------------------------------------------------")
            print('Estación: ',Clave[i])
            a=str(Estaciones[i][0])
            datos=[]
            os.system("awk 'NR>109{print $0}' "+str(a)+" > ../"+str(a)+".dat")
            aa=[]
            j=0
            lec=0
            with open('../'+str(a)+'.dat') as cvsfile:
                file=csv.reader(cvsfile,delimiter='	')
                while lec==0:
                    try:
                        for x1 in file:
                            aa.append(x1)
                        lec=1
                    except:
                        next(cvsfile, None)
                        j=j+1
            print(len(aa))
            if aa[len(aa)-1]==[]:
                aa.pop(len(aa)-1)
            comp1=np.zeros(len(aa))
            comp2=np.zeros(len(aa))
            comp3=np.zeros(len(aa))
            #if Clave[i] == 'CYK2':
            #    aa.pop(0)
            #    aa = aa[0:46800]
            #if Clave[i] == 'UKK2':
            #    aa.pop(0)
            #    aa = aa[0:37200]
            lengthoffile=len(aa)
            #if Clave[i] == 'CUP5':
            #       lengthoffile = 30099-109
            for jj in range(0,lengthoffile):
                b=aa[jj][0]
                comp=np.zeros(3)
                example_list = [k for k in b.split(' ')]
                ii=0
                for x1 in range(0,len(example_list)):
                    if example_list[x1]=='':
                        uselesse=1
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
            O1=Orientacion[i][0]
            O2=Orientacion[i][1]
            O3=Orientacion[i][2]
            O1,comp1=self.AsignacionDeOrientacion(O1,comp1)
            O2,comp2=self.AsignacionDeOrientacion(O2,comp2)
            O3,comp3=self.AsignacionDeOrientacion(O3,comp3)
            print(i+1,"Orientación de las componentes: ",O1,O2,O3)
            self.ASATOSAC(comp1,O1,Clave[i],Intervalo[i],Latitud[i],-Longitud[i],latsis,-lonsis,horasis,minsis,segsis,FECHA,rootsave)
            self.ASATOSAC(comp2,O2,Clave[i],Intervalo[i],Latitud[i],-Longitud[i],latsis,-lonsis,horasis,minsis,segsis,FECHA,rootsave)
            self.ASATOSAC(comp3,O3,Clave[i],Intervalo[i],Latitud[i],-Longitud[i],latsis,-lonsis,horasis,minsis,segsis,FECHA,rootsave)
            os.system("rm ../"+str(a)+".dat")
        os.system('rm ../Estacionesleidas.txt ../hpm.txt ../SDI.txt ../N.txt ../No.txt ../La.txt ../Lat.txt ../Lo.txt ../Lon.txt ../In.txt ../Insti.txt ../Ori.txt ../Inter.txt')
        os.system('rm ../Orientacion.txt ../Muestreo.txt')
        os.chdir('../')

    def ASATOSAC(self,Datos,Orientacion,CLAVE,Delta,stla,stlo,evla,evlo,hora,min,seg,FECHA,rootsave):
        from obspy import read
        from obspy.core import UTCDateTime
        st = read()
        st.remove(st[2])
        st.remove(st[1])
        st[0].data=Datos
        st[0].stats.station=str(CLAVE)
        st[0].stats.delta=Delta
        st[0].stats.channel=str(Orientacion)
        st[0].stats.network = 'RACM'
        U=UTCDateTime(str(FECHA)+'T'+str(int(hora))+':'+str(int(min))+':'+str(int(seg)))
        st[0].stats.starttime=U
        saq={
            'stla':stla,
            'stlo':stlo,
            'evla':evla,
            'evlo':evlo,
            }
        st[0].stats.sac=saq
        if Orientacion in ['Vert','East','North']:
            st.write('../'+str(rootsave)+str(CLAVE)+'.'+str(Orientacion[0])+'.sac')
        else:
            orint=Orientacion[0]+Orientacion[3]
            print(orint)
            if orint.lower() in ['sw','ne']:
                print('Orientacion: '+Orientacion+' Guardado como N')
                st.write('../'+str(rootsave)+str(CLAVE)+'.N.sac')
            elif orint.lower() in ['se','nw']:
            	print('Orientacion: '+Orientacion+' Guardado como E')
            	st.write('../'+str(rootsave)+str(CLAVE)+'.E.sac')

    def AsignacionDeOrientacion(self,O,DatosO):
        dum=1
        if O[0:4]=='N00E' or O[0:4]=='N00W':
            O='North'
            return O,DatosO
        elif O[0:4]=='S00E' or O[0:4]=='S00W':
            O='North'
            DatosO = DatosO*-1.0
            print(O,"mult*-1")
            return O,DatosO
        elif O[0:4]=='N90E' or O[0:4]=='S90E':
            O='East'
            return O,DatosO
        elif O[0:4]=='N90W' or O[0:4]=='S90W' or O[0:4]=='S90O' or O[0:4]=='N90O':
            O='East'
            DatosO = DatosO*-1.0
            print(O,"mult*-1")
            return O,DatosO
        elif O=='V' or O=='+V' or O=='+V;V' or O == 'V;+V' or O == '+V;+V' or O == 'V;V':
            O='Vert'
            return O,DatosO
        else:
            print('No se identifico correctamente esta componente. Se asignó el valor: '+str(O))
            #O='UNKNOWN'+str(dum)
            dum=dum+1
            return O,DatosO

    def CheckDelta(self):
        import os
        import statistics
        from obspy import read
        name=self.name
        root=self.root
        rootsave='./'+str(name)+'/'
        FECHA=self.date
        os.chdir(rootsave)
        V=read('*.V.*')
        E=read('*.E.*')
        N=read('*.N.*')
        V.detrend(type="linear")
        E.detrend(type="linear")
        N.detrend(type="linear")
        delta=[]
        for i in range(0,len(V)):
            delta.append(V[i].stats.delta)
        maxdelta=max(delta)
        mindelta=min(delta)
        mode=statistics.mode(delta)
        print(maxdelta,mode)
        print("..Revisando Intervalos de muestreo..")
        for i in range(0,len(V)):
            deltaisok = False
            while deltaisok == False:
                if delta[i]<0.01:
                    print("Estación: "+str(V[i].stats.station)+" Delta: "+str(delta[i]))
                    V[i].decimate(2,strict_length=False, no_filter=True)
                    E[i].decimate(2,strict_length=False, no_filter=True)
                    N[i].decimate(2,strict_length=False, no_filter=True)
                    print("Nuevo delta: "+str(V[i].stats.delta)+" "+str(E[i].stats.delta)+" "+str(N[i].stats.delta))
                    delta[i] = V[i].stats.delta
                    V[i].write(str(V[i].stats.station)+'.'+'V.sac')
                    E[i].write(str(E[i].stats.station)+'.'+'E.sac')
                    N[i].write(str(N[i].stats.station)+'.'+'N.sac')
                elif delta[i]>0.01:
                    print("Estación: "+str(V[i].stats.station)+" Delta: "+str(delta[i]))
                    V[i].interpolate(sampling_rate=100)
                    E[i].interpolate(sampling_rate=100)
                    N[i].interpolate(sampling_rate=100)
                    print("Nuevo delta: "+str(V[i].stats.delta)+" "+str(E[i].stats.delta)+" "+str(N[i].stats.delta))
                    delta[i] = V[i].stats.delta
                    V[i].write(str(V[i].stats.station)+'.'+'V.sac')
                    E[i].write(str(E[i].stats.station)+'.'+'E.sac')
                    N[i].write(str(N[i].stats.station)+'.'+'N.sac')
                elif delta[i] == 0.01:
                    if V[i].stats.delta == E[i].stats.delta and V[i].stats.delta == N[i].stats.delta:
                        deltaisok = True
                    else:
                        os.chdir('../')
                        print(V[i].stats.station,E[i].stats.station,N[i].stats.station)
                        raise ValueError('DELTA VALUES MISMATCH!')
        os.chdir('../')

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
        import matplotlib.pyplot as plt
        import numpy as np
        import os
        import sys
        from obspy import read
        import statistics as stcs
        name=self.name
        root=self.root
        rootsave='./'+str(name)+'/'
        FECHA=self.date
        os.chdir(rootsave)
        V=read('*.V.sac')
        E=read('*.E.sac')
        N=read('*.N.sac')
        V.detrend(type="linear")
        E.detrend(type="linear")
        N.detrend(type="linear")
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
                    j=i;
            else:
                if V[i].stats.station in stationname:
                    j=i;
        if j==None:
            mindist=min(dist)
            mindistindex=dist.index(min(dist))
            print("Minima distancia: "+clave[mindistindex])

        else:
            mindist=dist[j]
            mindistindex=j
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
        d2menos = int(min(data2correlatemenos))
        #data2correlatemas=28000
        #print(data2correlatemenos,data2correlatemas)
        cor=np.correlate(VF[mindistindex].data[d2menos:data2correlatemas[mindistindex]],VF[mindistindex].data[pulsemenos:pulsemas],mode='full')
        #cor=np.correlate(VF[mindistindex].data,VF[mindistindex].data[pulsemenos:pulsemas],mode='full')
        cor=cor/max(cor)
        maxargument=np.argmax(cor)

        print('Minima distancia: ',mindist,'Estación: ',clave[mindistindex])
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
                mindex1=np.argmax(a)

            while mindex1 > maxargument:
                a=np.delete(a,0)
                VF[i].data=np.delete(VF[i].data,0)
                V[i].data=np.delete(V[i].data,0)
                E[i].data=np.delete(E[i].data,0)
                N[i].data=np.delete(N[i].data,0)
                mindex1=np.argmax(a)

        rmenos = pulse[0] #int(input("-->Limite inferior de la fase a alinear: "))
        rmas = pulse[1] #input("-->Limite superior de la fase a alinear: "))
        maximo=np.argmax(VF[mindistindex].data[rmenos:rmas])
        rmenos = maximo + pulse[0] - 100
        rmas = maximo + pulse[0] + 100
        maximo=np.argmax(VF[mindistindex].data[rmenos:rmas])
        print('Alineacion de maximos')
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
                maximoi = np.argmax(VF[i].data[rmenos:rmas])
            while maximoi > maximo:
                VF[i].data=np.delete(VF[i].data,0)
                V[i].data=np.delete(V[i].data,0)
                E[i].data=np.delete(E[i].data,0)
                N[i].data=np.delete(N[i].data,0)
                maximoi = np.argmax(VF[i].data[rmenos:rmas])
        for i in range(0,len(V)):
            V[i].write(str(clave[i])+'.'+'V.sac')
            E[i].write(str(clave[i])+'.'+'E.sac')
            N[i].write(str(clave[i])+'.'+'N.sac')
        os.chdir('../')

    def plottraces(self,traces,bT,eT,Postraces=[0]):
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(len(Postraces),sharex=True)
        j=0
        if (len(Postraces) > 1):
            for i in Postraces:
                ax[j].plot(traces[i].data)
                ax[j].axvline(bT[i],color='red')
                ax[j].axvline(eT[i],color='blue')
                ax[j].set_ylabel(traces[i].stats.station)
                j=j+1
        else:
            ax.plot(traces.data)
            ax.axvline(bT,color='red')
            ax.axvline(eT,color='blue')
            ax.set_ylabel(traces.stats.station)
        plt.subplots_adjust(top=0.95, bottom=0.05,hspace=0)
        plt.show()

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
        self.plottraces(V,beginTrigger,endTrigger,Pos)
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
        VF=V.copy()
        EF=E.copy()
        NF=N.copy()
        clave=[]
        dist=[]
        velfreq=[]
        claveOrden=[]
        distOrden=[]
        velfreqOrden=[]
        l=len(V)
        #print(l)
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

        for i in range(0,l):
            while(len(V[i].data)<Nmax):
                V[i].data=np.insert(V[i].data,len(V[i].data),0)
                E[i].data=np.insert(E[i].data,len(E[i].data),0)
                N[i].data=np.insert(N[i].data,len(N[i].data),0)
            while (len(V[i].data)>Nmax):
                V[i].data=np.delete(V[i].data,len(V[i].data)-1)
                E[i].data=np.delete(E[i].data,len(E[i].data)-1)
                N[i].data=np.delete(N[i].data,len(N[i].data)-1)
                
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
            print(i+1,"Estación: "+str(claveOrden[i])+' Ajustada con: ',len(V[i].data))

        for i in range(0,l):
            V[i].write(str(claveOrden[i])+'.'+'V.sac')
            E[i].write(str(claveOrden[i])+'.'+'E.sac')
            N[i].write(str(claveOrden[i])+'.'+'N.sac')
        os.chdir('../')

    def back(self):
        import os
        os.chdir('..')

    def Rotar_sac(self):
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

    def ReduceTime(self,Traces,N="mean"):
        import statistics as stat
        import numpy as np
        try:
            from tqdm import tqdm
            tqdmisinstalled=True
        except:
            print('Tqdm module not installed')
            tqdmisinstalled=False
        if N == 'mean':
            duration=[]
            for i in Traces:
                duration.append(len(i.data))
            N = int(stat.mean(duration))
            duration=[]
            print(N)
        elif N == 'max':
            duration=[]
            for i in Traces:
                duration.append(len(i.data))
            N = int(max(duration))
            duration=[]
            print(N)
        if tqdmisinstalled:
            for i in tqdm(range(0,len(Traces))):
                while(len(Traces[i].data)<N):
                    Traces[i].data=np.insert(Traces[i].data,len(Traces[i].data),0)
                while (len(Traces[i].data)>N):
                    Traces[i].data=np.delete(Traces[i].data,len(Traces[i].data)-1)
        else:
            for i in range(0,len(Traces)):
                while(len(Traces[i].data)<N):
                    Traces[i].data=np.insert(Traces[i].data,len(Traces[i].data),0)
                while (len(Traces[i].data)>N):
                    Traces[i].data=np.delete(Traces[i].data,len(Traces[i].data)-1)
            #print("Estación: ",Traces[i].stats.station,'Ajustada con: ',Traces[i].stats.npts)
    
    def spectrograms_pulses(self,Stream,show=False,reshape=False):
        import scipy as sc
        from tqdm import tqdm
        import matplotlib.pyplot as plt
        X=[]
        for i in tqdm(range(0,len(Stream))):
            fs = Stream[i].stats.sampling_rate
            nps = 5000
            try:
                f,t,Sxx = sc.signal.spectrogram(Stream[i].data,fs,nperseg=nps,noverlap=nps-99,nfft=8192*3)
            except:
                from scipy.signal import spectrogram
                f,t,Sxx = spectrogram(Stream[i].data,fs,nperseg=nps,noverlap=nps-99,nfft=8192*3)
            f=f[0:50]
            Sxx = Sxx[0:50]
            dum = []
            for j in range(0,50):
                dum.append(max(Sxx[j]))
            Sxx = Sxx/max(dum)
            if reshape:
                X.append(Sxx.reshape(50*405))
            else:
                X.append(Sxx)
            if show:
                print(i,len(Sxx),len(Sxx[0]),len(f),len(t))
                c = plt.pcolormesh(t,f, Sxx,shading='auto',cmap='seismic')
                plt.colorbar(c)
                plt.title(Stream[i].stats.station)
                #plt.axvline(St[i].stats.sac.a,color='red')
                #plt.axvline(St[i].stats.sac.t1,color='red')
                plt.show()
        return X

    def fit_pulses(self,show=False):
        from obspy import read
        from tensorflow import keras
        import os
        import numpy as np
        os.chdir(self.name)
        V = read('*V.sac')
        E = read('*E.sac')
        N = read('*N.sac')
        clave = []
        for i in V:
            clave.append(i.stats.station)
        self.ReduceTime(V,N=self.nlength)
        self.ReduceTime(N,N=self.nlength)
        self.ReduceTime(E,N=self.nlength)
        VF=V.copy()
        VF.filter("bandpass",corners=4, freqmin=0.0833, freqmax=0.125, zerophase=True)
        X = self.spectrograms_pulses(VF)
        X = np.array(X)
        model = keras.models.load_model('/home/gilberto/ANN_Model_pulses2.tf')
        y_pred = model.predict(X)
        yb = []
        ye = []
        j=0
        for i in y_pred:
            if i[0] < 0:
            	i[0] = 0
            yb.append(i[0]*self.nlength)
            ye.append(i[1]*self.nlength)
            V[j].stats.sac['a'] = i[0]*450
            V[j].stats.sac['t1'] = i[1]*450
            j=j+1
        Pos = np.arange(0,len(V),1)
        if show:
            self.plottraces(V,yb,ye,Pos)
        for i in range(0,len(V)):
            V[i].write(clave[i]+'.V.sac')
            E[i].write(clave[i]+'.E.sac')
            N[i].write(clave[i]+'.N.sac')
        os.chdir('../')

    def get_fecha(self):
    	fecha = ['enero','febrero','marzo','abril','mayo','junio','julio','agosto','septiembre','octubre','noviembre','diciembre']
    	pass
    
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
                yb[i]=yb[i]/(self.nlength/100)
                ye[i]=ye[i]/(self.nlength/100)
            dum[i][0] = yb[i]
            dum[i][1] = ye[i]
            y.append(dum[i])
        return y

    def CheckTraces(self):
        import os
        from obspy import read
        rootsave='./'+str(self.name)+'/'
        os.chdir(rootsave)
        V=read('*.V.sac')
        E=read('*.E.sac')
        N=read('*.N.sac')
        print('Componentes Verticales leidas: ',len(V))
        print('Componentes Norte leidas: ',len(N))
        print('Componentes Este leidas: ',len(E))
        if len(V) == len(N) and len(V) == len(E):
            print('TODO EN ORDEN')
        else:
            print('REVISAR TRAZAS')
        os.chdir('../')

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
    
    def plot_tz(self,x,y,z,ranges=[-180,180],dt=20,s=50,log=True,set_limits=[[0.070,0.6],[0,300]],xlabel='Time [s]',ylabel='Freq [Hz]',zlabel='BAZ [DEG]',pshow = True,colormap='hsv'):
        import matplotlib.pyplot as plt
        import numpy as np
        fig, ax = plt.subplots(2,1,sharex=True,figsize=(9,16),gridspec_kw={'height_ratios': [1, 4]})
        cm = plt.cm.get_cmap(colormap)
        sc = ax[1].scatter(x,y, c=z,s=s, vmin=ranges[0], vmax=ranges[1], cmap=cm,zorder=100)
        cbarticks=np.arange(ranges[0],ranges[1]+1,dt)
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