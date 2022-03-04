def fit_pulses(V,show=False,N=4500):
    from obspy import read
    from tensorflow import keras
    import os
    import numpy as np
    from pathlib import Path
    BASE_DIR = Path(__file__).resolve().parent.parent
    clave = []
    for i in V:
        clave.append(i.stats.station)
    ReduceTime(V,N=N)
    VF=V.copy()
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
        yb.append(i[0]*N)
        ye.append(i[1]*N)
        V[j].stats.sac['a'] = i[0]*450
        V[j].stats.sac['t1'] = i[1]*450
        j=j+1
    Pos = np.arange(0,len(V),1)
    if show:
        plottraces(V,yb,ye,Pos)
    return V

def ReduceTime(Traces,N="mean"):
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
        print("Reduciendo el tiempo de componente: ",)
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

def plottraces(traces,bT,eT,Postraces=[0]):
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

def spectrograms_pulses(Stream,show=False,reshape=False):
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
            plt.show()
    return X