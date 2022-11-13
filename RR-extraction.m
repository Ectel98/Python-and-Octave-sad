clear

Fs  = 100;

pkg load signal

for u=0:400

c = dlmread(["C:\\Users\\Denis\\Documents\\Uni\\Apnea sleep detector\\Filtraggio\\Dati\\dati physionet\\divisi a19\\a19-"  num2str(u)  ".txt"])

c = c.-mean(c)

#figure

#plot(c)

f1 = 15 #100-125    15

f2 = 26 #125-130    16

d_f = f2-f1;

dB = 40;

N = dB*Fs/(22*d_f);

f = f1/(Fs/2);

low_filt = fir1(round(N)-1, f,'low');

#figure

#plot((-0.5:1/4096:0.5-1/4096)*Fs,20*log10(abs(fftshift(fft(low_filt,4096)))))

y2 = filter(low_filt,1,c);

#plot(y2)

#---------------------

f1 = 1;

f2 = 5;

d_f = f2- f1;

dB = 20;

N = dB*Fs/(22*d_f);

f = f2/(Fs/2);

high_filt = fir1(round(N)-1,f,'high');

y3 = filter(high_filt,1,y2);

#figure;

#plot(y3);

#-------------------------------------------------------------

#f = 50/250 = 0.20 

#f1 = f-0.02;
 
#f2 = f+0.02;

#N = 100

#notch = fir1(N,[0.23 0.26],'stop'); 

#y = filter(notch,1,y3);

#figure

#plot(y)

c = y3

#c = dlmread('C:\Users\Denis\Desktop\spallafilt.txt')

#-----------------------------------

for i=5:(length(c)-1) #6 -> 
    Y(i-4) = (c(i-4) + 2*c(i-3) - 2*c(i+1) - c(i))/8;
end


#figure;

#plot(Y)


#-------------------------------------

#n=5
#for i=5+1:length(c) #6 -> 
#    Y(i-n) = (c(i) - c(i-n))/n; #6-1
    #Y(i-4) = (2*c(i) + c(i-1) - c(i-3) - 2*c(i-4))
#end

#figure;

#plot(Y)

#------------------------------------

x = Y .^ 2;

#figure;

#plot(x)

#-------------------------------------

xn= x;

N = (150/(1/Fs))/1000;
xf= zeros(size(xn));
xfgnu = zeros(size(xn));
fk = 1/N*ones(N,1);


#for idx = N:length(xn)
#  xf(idx) = sum(xn(idx-N+1:idx))/N; % Average calculation!
#end

xfgnu = filter(fk,1,xn); % Z-Transform approach

#figure
#plot(xfgnu)
#figure

#-----------------------------------------

d = xfgnu

n_to_test = 2000#2/(1/500)

e = d(1:n_to_test)


n_i = mean(e)
s_i = max(e)*0.875



t = 2*pi*linspace(0,1,1024)';
dt = t(2)-t(1)

[pks loc] = findpeaks(d,"DoubleSided","MinPeakDistance",round(0.1/dt))  #La distanza minima camabia in base alla frequenza di campionamento! 0.5 500hz - 0.1 100hz


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

#hold on

#scatter(lc,pk) #Plotta i picchi R

#-------------------------------------------------------------------------------

uno = []

for i=1:length(pk)-1 
  rr_interval(i) = (lc(i+1)-lc(i))*10          #[ms]
endfor

save "-append" "-ascii" 'C:\Users\Denis\Documents\Uni\Apnea sleep detector\Filtraggio\Dati\dati physionet\intervalli risultanti a19\a19_rr_interval.txt' rr_interval

rr_interval = []

endfor
#-------------------------------------------------------------------------------
