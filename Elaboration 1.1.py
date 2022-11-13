import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import hilbert
from math import *
import time


a = 'b02'

print(a)

f = open('intervalli risultanti '+a+'\\'+a+'_rr_interval.txt',"r+") #a01 - a02 - a03 - a04 - a06 - a08 - a09 / b01 - b02 - b03 - b04 - b05/ c01

l = f.readlines()

f.close()
#------------------------------------------------------------------------------------
#Normalizzaione degli intervalli: vengono rimossi gli intervalli <400ms e >2000ms



out=0

i=0
while i in range(0,len(l)):
    if (l[i]=='\n' or l[i]==''):                #Si sono creati dall'elbaorazione di più file!
        l.remove(l[i])
    if (float(l[i])<400 or float(l[i])>2000): #Problema: nel calcolo della soglia prende i punti di massimo come gli intervalli massimi
        l.remove(l[i])
        out+=1
    else:
        i+=1

print("Dati fuori dal range:"+str(out) + '\n' + "Dati presenti:" + str(len(l)))

i=0

for i in range(0,len(l)):
    l[i]=(float(l[i]))

print("Spete un moment...")


#----------------------------------------------------------------------------
#Creazione vettore temporale
        


t = [0]*(len(l))
ts = 0

for i in range(0,len(l)):
    ts += l[i]
    t[i] = ts


xpoints = np.array(t)
ypoints = np.array(l)

plt.figure()

plt.plot(xpoints, ypoints)

#--------------------------------------------------------------------------
#Da millisecondi a secndi

for i in range(0,len(t)):
    t[i]=t[i]/1000
    l[i]=l[i]/1000

#------------------------------------------------------------------------------
#Filtro intervalli irregolari



sumi = 0
filt = 0.2
hwin = 20
win = hwin*2

t_filt = []
l_filt = []

j = hwin

for i in range(0,win+1):
    sumi += l[i]

    
sumi -= l[hwin]
av = sumi/win
sumi += l[hwin] - l[0]

filtmax = filtmin = filt * av;

if (l[hwin] < av+filtmax and l[hwin] > av-filtmin):
    t_filt.append(t[hwin])
    l_filt.append(l[hwin])

e = win+1
while (e < len(l)):

    sumi += l[e] - l[e-hwin]   #Aggiungo il prossimo dato e tolgo il nuovo centrale
    av = sumi/win   
        
    sumi += l[e-hwin] - l[e-win]

    filtmax = filtmin = filt * av

    if (l[e-hwin] < av+filtmax and l[e-hwin] > av-filtmin):
        t_filt.append(t[e-hwin])
        l_filt.append(l[e-hwin])
    e+=1

t = list(t_filt)
l = list(l_filt)

plt.figure()

xpoints = np.array(t)
ypoints = np.array(l)

plt.plot(xpoints, ypoints)


#--------------------------------------------------------------
#Interpolazione



l_intrep = []
time_intrep = []

x0 = t[0] + 1;

i = 0

while (i<(len(l)-2)):

    while (x0 > t[i] and i<(len(l)-2)):

        i+=1

    b = (l[i+1] - l[i]) / (t[i+1] - t[i])
    a = l[i] - b * t[i]

    y0 = (b * x0 + a)
    

    l_intrep.append(y0)
    time_intrep.append(x0)

    x0+=1


xpoints = np.array(time_intrep)
ypoints = np.array(l_intrep)

#plt.figure()

#plt.plot(xpoints, ypoints)



#--------------------------------------------------------------
#Rimozione della tendenza



l_detrend = []
time_detrend = []

hwin = 40;
win = 2*hwin +1;


sumx = sumy = sumxy = sumx2 = 0;
        

for i in range(0,win): 
    sumx += time_intrep[i]
    sumy += l_intrep[i]
    sumxy += time_intrep[i]*l_intrep[i]
    sumx2 += time_intrep[i]*time_intrep[i]
        

b = (sumxy - sumx*sumy/win) / (sumx2 - sumx*sumx/win)                                #Calcolo la tendenza
a = sumy/win - b*sumx/win

for i in range (0,hwin+1):                                                           #Sottraggo la tendenza
    l_detrend.append(l_intrep[i] - (a + b * time_intrep[i]))
    time_detrend.append(time_intrep[i])
        

for i in range (win,len(l_intrep)):
    sumx += time_intrep[i]-time_intrep[i-win]                                        #Rimuvo i vecchi acampioni e aggiungo i nuovi
    sumy += l_intrep[i]-l_intrep[i-win]
    sumxy += time_intrep[i]*l_intrep[i]-time_intrep[i-win]*l_intrep[i-win]
    sumx2 += time_intrep[i]*time_intrep[i]-time_intrep[i-win]*time_intrep[i-win]

    b = (sumxy - sumx*sumy/win) / (sumx2 - sumx*sumx/win)                            #Calcolo la tendenza
    a = sumy/win - b*sumx/win

    l_detrend.append(l_intrep[i-hwin] - (a + b*time_intrep[i-hwin]))              #Sottraggo la tendenza
    time_detrend.append(time_intrep[i-hwin])
       

for i in range(i-hwin,len(l_intrep)):                                                 #Sottraggo la tendenza
    l_detrend.append(l_intrep[i] - (a + b * time_intrep[i]))
    time_detrend.append(time_intrep[i])


xpoints = np.array(time_detrend)
ypoints = np.array(l_detrend)

#plt.figure()

#plt.plot(xpoints, ypoints)

#--------------------------------------------------------------
#Smooth



i = 0

l_smooth = []
time_smooth = []


win = 5
sumx = sumy = 0

while (i<win):
    sumx += time_detrend[i]
    sumy += l_detrend[i]
    i+=1
        

l_smooth.append(sumy/win)
time_smooth.append(sumx/win)

sumx -= time_detrend[0]
sumy -= l_detrend[0]

while (i<len(time_detrend)):

    sumx += time_detrend[i]
    sumy += l_detrend[i]

    l_smooth.append(sumy/win)
    time_smooth.append(sumx/win)

    sumx -= time_detrend[i-win+1]
    sumy -= l_detrend[i-win+1]

    i+=1

xpoints = np.array(time_smooth)
ypoints = np.array(l_smooth)

plt.figure()

plt.plot(xpoints, ypoints)
#--------------------------------------------------------------
#Trasformata di Hilbert



lfilt = 128;

xh =[0]*len(time_smooth)
ampl =[0]*len(time_smooth)
phase =[0]*len(time_smooth)
omega =[0]*len(time_smooth)
phase = [0]*len(time_smooth)
hilb = [0]*(lfilt+1)

pi2 = 2*pi
npt=len(time_smooth)-1



for i in range (1,lfilt+1):
    hilb[i] = (float) (1 / ((i - lfilt / 2) - 0.5) / pi)

#Convoluzione
for e in range (1,npt-lfilt+2):
    yt = 0
    for i in range (1,lfilt+1):
        yt = (yt+l_smooth[e+i-1]*hilb[lfilt+1-i])
    xh[e] = yt


for i in range (1,npt-lfilt+1):
    xh[i] = (float) (0.5*(xh[i]+xh[i+1]))
        

# shifting lfilt/1+1/2 points #
for i in range (npt-lfilt,0,-1):
    xh[i+lfilt//2]=xh[i]


# Ampl and phase 
for i in range(lfilt//2+1,(npt-(lfilt//2))+1):       #Da sistemare nell'app java
    xt = l_smooth[i]
    xht = xh[i]
    ampl[i] = (float)(sqrt(xt*xt+xht*xht))
    phase[i] = (float)(atan2(xht ,xt))
    if (phase[i] < phase[i-1]):
        omega[i]=(phase[i]-phase[i-1]+pi2)
    else:
        omega[i]=(phase[i]-phase[i-1])
        

for i in range(0,len(omega)):
    omega[i] = omega[i] / pi2

time_ht = time_smooth[lfilt//2+2:npt-lfilt//2+1]
omega = omega[lfilt//2+2:npt-lfilt//2+1]
ampl = ampl[lfilt//2+2:npt-lfilt//2+1]

analytic_signal = hilbert(l_smooth)

amplitude_envelope = np.abs(analytic_signal)
instantaneous_phase = np.unwrap(np.angle(analytic_signal))
instantaneous_frequency = (np.diff(instantaneous_phase) /(2.0*np.pi))

#plt.figure()

#plt.plot(time_ht,omega)

plt.figure()

plt.plot(time_ht,ampl)



#-----------------------------------------------------------------
#Filtro a finestra mobile trasformata di Hilber



win=60
v_sx =[0]*win
v_sy =[0]*win
sx =[]
sy = []

omega_filt = []
amp_filt = []
time_ht_filt = []


u = 0
j = hwin = (win//2)-1

for k in range(0,win):
    sx.append(ampl[k])
    sy.append(omega[k])
        
v_sx = list(sx)
v_sy = list(sy)


sx.sort()
sy.sort()

time_ht_filt.append(time_ht[j])
amp_filt.append(sx[hwin])
omega_filt.append(sy[hwin])


for e in range(win,len(time_ht)):

    for k in range(1,win):           #Shift
        v_sx[k-1] = v_sx[k]          #
        v_sy[k-1] = v_sy[k]          #


    v_sx[win-1] = ampl[e]            #Aggiungo valore
    v_sy[win-1] = omega[e]           #
            
    j+=1

    sx = list(v_sx)
    sy = list(v_sy)
            
    sx.sort()
    sy.sort()
    

    time_ht_filt.append(time_ht[j])  
    amp_filt.append(sx[hwin])
    omega_filt.append(sy[hwin])

    

#plt.figure()

#plt.plot(time_ht_filt,omega_filt)

plt.figure()

plt.plot(time_ht_filt,amp_filt)
        

#--------------------------------------------------------------
#Normalizzazione delle ampiezze


amp_norm = []

av_amp = 0

for i in range (0,len(time_ht_filt)):
    av_amp += amp_filt[i]   

av = av_amp/len(time_ht_filt)

for i in range(0,len(time_ht_filt)):
    amp_norm.append(amp_filt[i]/av)

plt.figure()

plt.plot(time_ht_filt,amp_norm)


#--------------------------------------------------------------
#Calcola il miniomo e la soglia



ord_ampl= list(amp_norm)

ord_ampl.sort()

miin = ord_ampl[0]
maax = ord_ampl[len(ord_ampl)-100]

mid = (maax + miin)/2

thres = (-0.555 + 1.3*(mid+1)/2)

print(thres)

#--------------------------------------------------------------
#Deviazione stadard ecc...


start = []
avy = []
sdy = []
ydetw = []
avz  = []
sdz = []
zdetw = []

e = 1

incr = 60     #Incremento

win = 300

v_ydet = v_zdet = 0


v_start = time_ht_filt[0]
sumy = amp_norm[0]
sumz = omega_filt[0]
sumyy = amp_norm[0]*amp_norm[0]
sumzz = sumz*sumz

if (amp_norm[0] >= thres):
    v_ydet+=1
            
if (omega_filt[0] <= 0.06):
    v_zdet+=1


for i in range(1,win):          #Prende la prima finestra e calcola il primo dato

    sumy += amp_norm[i]
    sumz += omega_filt[i]
    sumyy += amp_norm[i]*amp_norm[i]
    sumzz += omega_filt[i]*omega_filt[i]

    if (amp_norm[i] >= thres):
        v_ydet+=1
    if (omega_filt[i] <= 0.0555):
        v_zdet+=1
                

v_avy = sumy/win
v_avz = sumz/win
v_sdy = sqrt((sumyy - sumy*sumy/win)/(win-1))
v_sdz = sqrt((sumzz - sumz*sumz/win)/(win-1))

start.append(v_start)
avy.append(v_avy)
sdy.append(v_sdy)
ydetw.append(v_ydet/win)
avz.append(v_avz)
sdz.append(v_sdz)
zdetw.append(v_zdet/win)


for j in range(0,incr):      #Sottrae il contributo dell'incremento
    sumy -= amp_norm[j]
    sumz -= omega_filt[j]
    sumyy -= amp_norm[j]*amp_norm[j]
    sumzz -= omega_filt[j]*omega_filt[j]

    if (amp_norm[j] >= thres):
        v_ydet-=1
    if (omega_filt[j] <= 0.06):
        v_zdet-=1
        

v_start += incr;                       #Aggionra in tempo di start


i=win
e = 1
while (i<len(time_ht_filt)):

    while (time_ht_filt[i]< v_start and i<len(time_ht_filt)):  #Serve davvero ? non è sempre vera ?
        i+=1

    if (len(time_ht_filt)-i<incr):
        break

    sumy += amp_norm[i]                   #Aggiungo i contributi di un nuovo elemento
    sumz += omega_filt[i]
    sumyy += amp_norm[i]*amp_norm[i]
    sumzz += omega_filt[i]*omega_filt[i]

    if (amp_norm[i] >= thres):
        v_ydet+=1;
    if (omega_filt[i] <= 0.06):
        v_zdet+=1;

    for j in range(i+1,i+incr):                #Aggiunge il contributo di un incremento meno l'elemento preso sopra

        sumy += amp_norm[j]
        sumz += omega_filt[j]
        sumyy += amp_norm[j]*amp_norm[j]
        sumzz += omega_filt[j]*omega_filt[j]

        if (amp_norm[j] >= thres):
            v_ydet+=1
        if (omega_filt[j] <= 0.06):
            v_zdet+=1

    i+=incr

    v_avy = sumy/win                                       #Calcola i parametri
    v_avz = sumz/win
    v_sdy = sqrt((sumyy - sumy*sumy/win)/(win-1))
    v_sdz = sqrt((sumzz - sumz*sumz/win)/(win-1))

    start.append(v_start)
    avy.append(v_avy)
    sdy.append(v_sdy)
    ydetw.append(v_ydet/(win))
    avz.append(v_avz)
    sdz.append(v_sdz)
    zdetw.append(v_zdet/(win))

    e +=1 
    
    for j in range(i-win,(i-win)+incr):                                 #Rimuove l'incremetno 
        sumy -= amp_norm[j]
        sumz -= omega_filt[j]
        sumyy -= amp_norm[j]*amp_norm[j]
        sumzz -= omega_filt[j]*omega_filt[j]
        if (amp_norm[j] >= thres):
            v_ydet-=1
        if (omega_filt[j] <= 0.06):
            v_zdet-=1
            
    v_start += incr;

#--------------------------------------------------------
#Applicazione dei limiti


AVAMP0=0.65
AVAMP1=2.5
SDAMP0=0
SDAMP1=0.6
AMPTIME0=0.006
AMPTIME1=1
AVFREQ0=0.01
AVFREQ1=0.055
SDFREQ0=0
SDFREQ1=0.01
FREQTIME0=0.7
FREQTIME1=1

p_times = []


for i in range (0,len(start)):
    if (avy[i] >= AVAMP0 and avy[i] <= AVAMP1) and (sdy[i] >= SDAMP0 and sdy[i]<= SDAMP1) and (ydetw[i] >= AMPTIME0 and ydetw[i] <= AMPTIME1) and (avz[i] >= AVFREQ0 and avz[i] <= AVFREQ1) and (sdz[i] >= SDFREQ0 and sdz[i] <= SDFREQ1) and (zdetw[i] >= FREQTIME0 and zdetw[i] <= FREQTIME1):
        p_times.append(start[i])
    #print(str(AVAMP0) + ' < '  + str(avy[i]) + ' > ' + str(AVAMP1) + ' -- ' + str(SDAMP0) + ' < '  + str(sdy[i]) + ' > ' +str(SDAMP1)+ ' -- ' +str(AMPTIME0) + ' < '  + str(ydetw[i]) + ' > ' +str(AMPTIME1)+ ' -- ' +str(AVFREQ0) + ' < '  + str(avz[i]) + ' > ' +str(AVFREQ1)+ ' -- '+ str(SDFREQ0) + ' < '  + str(sdz[i]) + ' > ' + str(SDFREQ1)+ ' -- ' +str(FREQTIME0) + ' < '  + str(zdetw[i]) + ' > ' +str(FREQTIME1)+'\n') 

#--------------------------------------------------------------
#Determinazione intervalli



incr = 60
win = 300
mini = 900


start_time = []
end_time = []

runflag = runflag0 = False

runstart = p_times[0]
lasttime = p_times[0]

runend0 = runstart0 = sumi = 0

i = 1


while (i<len(p_times)): 

    timep = p_times[i]

    if (timep - lasttime != incr):                           #Mi porta alla fine della sequanza di intervalli positivi consecutivi
        if (lasttime - runstart + win >= mini):
            if runflag0 is False:
                runflag0 = True
                runstart0 = runstart
                runend0 = lasttime
            elif (runstart <= runend0 + win):
                runend0 = lasttime
            else:
                start_time.append(time.strftime('%H:%M:%S', time.gmtime(runstart0)))
                end_time.append(time.strftime('%H:%M:%S', time.gmtime(runend0 + win)))
                sumi += runend0-runstart0+win
                runstart0 = runstart
                runend0 = lasttime
        runstart = timep
    lasttime = timep

    i+=1


if (lasttime - runstart + win >= mini):
    runflag = True

if (runflag0):
    if (runflag and mini > win and runstart <= runend0 + win):
        runend0 = lasttime
        runflag = False
    start_time.append(time.strftime('%H:%M:%S', time.gmtime(runstart0)))
    end_time.append(time.strftime('%H:%M:%S', time.gmtime(runend0 + win)))
    sumi += runend0-runstart0+win

if (runflag):
    start_time.append(time.strftime('%H:%M:%S', time.gmtime(runstart)))
    end_time.append(time.strftime('%H:%M:%S', time.gmtime(lasttime + win)))
    sumi += lasttime-runstart+win
    
somma = time.strftime('%H:%M:%S', time.gmtime(sumi))


for i in range(0,len(start_time)):
    print("Start: "+ start_time[i]+ "End: " + end_time[i])
print("Somma:" + somma, sumi)

#plt.show()
