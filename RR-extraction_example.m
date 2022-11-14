clear

Fs  = 500;

pkg load signal

#Lettura

c = dlmread(["C:\\Users\\Denis\\Documents\\Uni\\Apnea sleep detector\\Filtraggio\\Dati\\spalla.txt"])

c = c.-mean(c)

figure

plot(c)

#Filtro passa-alto

f1 = 15 

f2 = 25 

d_f = f2-f1;

dB = 40;

N = dB*Fs/(22*d_f);

f = f1/(Fs/2);

low_filt = fir1(round(N)-1, f,'low');

figure

#plot((-0.5:1/4096:0.5-1/4096)*Fs,20*log10(abs(fftshift(fft(low_filt,4096)))))

y2 = filter(low_filt,1,c);

plot(y2)

#Filtro passa-alto

f1 = 1;

f2 = 5;

d_f = f2- f1;

dB = 20;

N = dB*Fs/(22*d_f);

f = f2/(Fs/2);

high_filt = fir1(round(N)-1,f,'high');

c = filter(high_filt,1,y2);

figure;

plot(c);

#Filtro derivatore

for i=5:(length(c)-1)  
    Y(i-4) = (c(i-4) + 2*c(i-3) - 2*c(i+1) - c(i))/8;
end


figure;

plot(Y)


#Squadratura

x = Y .^ 2;

figure;

plot(x)

#Finestra mobile

xn= x;

N = (150/(1/Fs))/1000;
xf= zeros(size(xn));
xfgnu = zeros(size(xn));
fk = 1/N*ones(N,1);

xfgnu = filter(fk,1,xn); 

figure
plot(xfgnu)


#Picchi

d = xfgnu

n_to_test = 2000#2/(1/500)

e = d(1:n_to_test)


n_i = mean(e)
s_i = max(e)*0.875



t = 2*pi*linspace(0,1,1024)';
dt = t(2)-t(1)

[pks loc] = findpeaks(d,"DoubleSided","MinPeakDistance",round(0.5/dt))  #La distanza minima camabia in base alla frequenza di campionamento - 0.5 500hz - 0.1 100hz


pk = []
lc = []
e = 1
i = 1

for i=1:(length(pks))

  ts = n_i + 0.25*(s_i-n_i)

  if (pks(i)>ts)
    s_i = 0.125*pks(i) + 0.875*s_i
    pk(e) = pks(i)
    lc(e) = loc(i)
    e = e + 1
  endif

  if (pks(i)<=ts)
    n_i = 0.125*pks(i) + 0.875*n_i
  endif
  
endfor

hold on

scatter(lc,pk) #Plotta i picchi R

#-------------------------------------------------------------------------------

for i=1:length(pk)-1 
  rr_interval(i) = (lc(i+1)-lc(i))*10          #[ms]
endfor

rr_interval = []
#-------------------------------------------------------------------------------
