import sys
sys.path.insert(1,"../")
from seismoTK import RACM
#Instanciamos a la clase RACM
#Nombre del directorio a guardar, nombre del directorio con los datos, Fecha del sismo.
RC=RACM('2021-09-08-sac-2','2021-09-08','2021-09-08')
#Convierte los datos de formato ASA a SAC
RC.ASA2SAC()
#Revisa que todos los datos tengan el mismo intervalo de muestreo.
RC.CheckDelta()
#print(RC.V[0].stats)
RC.CheckTraces()
#Rota las componentes con un script de sac. (Solo ubuntu)
#RC.Rotar_sac()
#Delimita los pulsos superficiales a partir de un modelo de ML. Si no quieres utilizarlo, deberas limitar los pulsos picando en sac con a para el inicio y t1 para el final.
RC.fit_pulses()
#Grafica los pulsos con su respectivo picado. Revisa que la identificacion del pulso se realizo de manera correcta.
RC.CheckAlinear()
#Alinea los pulsos con los picados anteriores.
RC.Alinear(stationname=["CUP5"])
#Revisa que la alineacion se llevo a cabo de manera correcta. Este metodo NO mueve los picados, por lo tanto, ya no corresponden al pulso. Si quieres que correspondan, ejecuta el metodo fit_pulses nuevamente.
RC.CheckAlinear()
#Asigna la velocidad de 3.2 segun la distancia estacion-epicentro.
RC.fit_pulses()
#RC.Velocidad()
#Los metodos se encuentran en el script llamad RACM.py dentro  de la carpeta seismoTK.