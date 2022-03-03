def ASATOSAC(Datos,Orientacion,CLAVE,Delta,stla,stlo,evla,evlo,hora,min,seg,FECHA,rootsave):
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
        try:
            orint=Orientacion[0]+Orientacion[3]
        except:
            orint=Orientacion
        print(orint)
        if orint.lower() in ['sw','ne']:
            print('Orientacion: '+Orientacion+' Guardado como N')
            st.write('../'+str(rootsave)+str(CLAVE)+'.N.sac')
        elif orint.lower() in ['se','nw']:
            print('Orientacion: '+Orientacion+' Guardado como E')
            st.write('../'+str(rootsave)+str(CLAVE)+'.E.sac')
        elif orint == Orientacion:
            st.write('../'+str(rootsave)+str(CLAVE)+'.'+Orientacion+'.sac')

def AsignacionDeOrientacion(O,DatosO):
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
        dum=dum+1
        return O,DatosO