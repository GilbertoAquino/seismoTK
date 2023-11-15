import sys
import matplotlib.pyplot as plt
import numpy as np
sys.path.insert(1,"../")
from seismoTK import RACM
RC=RACM('2021-09-08-sac','2021-09-08','2021-09-08')
for i in range(0,len(RC.V)):
    time = np.arange(0,RC.V[i].stats["delta"]*RC.V[i].stats["npts"],RC.V[i].stats["delta"])
    fig, ax = plt.subplots(2,sharex=True)
    ax[0].plot(time,RC.V[i].data,color="black",lw=0.6)
    Vf = RC.V[i].copy()
    Vf.filter("bandpass",corners=4, freqmin=0.0833, freqmax=0.125, zerophase=True)
    ax[1].plot(time,Vf,color="black",lw=0.8)
    ax[0].axvline(RC.V[i].stats["sac"]["a"],color='red')
    ax[0].axvline(RC.V[i].stats["sac"]["t1"],color='blue')
    ax[1].axvline(RC.V[i].stats["sac"]["a"],color='red')
    ax[1].axvline(RC.V[i].stats["sac"]["t1"],color='blue')
    ax[0].set_ylabel("Amplitud [m/s/s]")
    ax[1].set_xlabel("Tiempo [s]")
    plt.subplots_adjust(top=0.8, bottom=0.2,hspace=0)
    plt.title(RC.V[i].stats.station)
    plt.show()