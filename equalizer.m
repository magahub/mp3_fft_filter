clear;
clc;
%% MAIN CODE
filename = 'bg.mp3';
[X,Fs] = audioread(filename);
info = audioinfo(filename);
T = 1/Fs;             % Sampling period       
L = info.TotalSamples;             % Length of signal
t = (0:L-1)*T;        % Time vector

Y = fft(X);
s = 1:20000;

para = round(L/numel(s));

transfer1=zeros(numel(s),2);
transfer1(:,1)=10.*(s./500.5 +1)./(s./50.5+1)./(s./2121.5 +1);
transfer1(:,2)=10.*(s./500.5 +1)./(s./50.5+1)./(s./2121.5 +1);

transfer2=zeros(numel(s),2);
transfer2(:,1)=1./10./(s./500.5 +1).*(s./50.5+1).*(s./2121.5 +1);
transfer2(:,2)=1./10./(s./500.5 +1).*(s./50.5+1).*(s./2121.5 +1);

filter1 = zeros(numel(Y)/2,2);
w = zeros(numel(Y)/2,1);
j=1;
for i = 1:numel(Y)/2
    
    filter1(i,1)=transfer1(j,1);
    filter1(i,2)=transfer1(j,1);
    w(i) = j;
    if mod(i,para) == 0
        if j < numel(s)
            j=j+1;
        end
    end
end

filter2 = zeros(numel(Y)/2,2);
j=1;
for i = 1:numel(Y)/2
    
    filter2(i,1)=transfer2(j,1);
    filter2(i,2)=transfer2(j,1);
    if mod(i,para) == 0
        if j < numel(s)
            j=j+1;
        end
    end
end

filtered1 = filter1.*Y;
filtered2 = filter2.*Y;
equalized = filter2.*filtered1;


%?????1
K=Y;
P2 = abs(K/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
subplot(2,2,1)
semilogx(f,P1) 
xlim([20 20000])
title('Original')
xlabel('f (Hz)')
ylabel('amplitude')

%?????2
K=filtered1;
P2 = abs(K/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
subplot(2,2,2)
semilogx(f,P1) 
xlim([20 20000])
title('only RIAA')
xlabel('f (Hz)')
ylabel('amplitude')

%?????3
K=filtered2;
P2 = abs(K/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
subplot(2,2,3)
semilogx(f,P1) 
xlim([20 20000])
title('only Equalize')
xlabel('f (Hz)')
ylabel('amplitude')

%?????4
K=equalized;
P2 = abs(K/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
subplot(2,2,4)
semilogx(f,P1) 
xlim([20 20000])
title('RIAA + Equalize')
xlabel('f (Hz)')
ylabel('amplitude')

comeback1 = ifft(filtered1);
comeback2 = ifft(filtered2);
comeback3 = ifft(equalized);
audiowrite(strcat(filename,'_RIAA.wav'),comeback1,Fs);
audiowrite(strcat(filename,'_Only_Equalize.wav'),comeback1,Fs);
audiowrite(strcat(filename,'_RIAA+Equalize.wav'),comeback3,Fs);


%% RIAA Filter plot code

P2 = filter1;
P1 = P2(1:L/2+1);
P1(2:end-1) = P1(2:end-1);
f = Fs*(0:(L/2))/L;
P1 = 20.*log10(P1);
semilogx(f,P1) 
axis([20,22000,-20,20])
title('RIAA Filter')
xlabel('f (Hz)')
ylabel('amplitude')

%% EQUALIZER Filter plot code

P2 = filter2;
P1 = P2(1:L/2+1);
P1(2:end-1) = P1(2:end-1);
f = Fs*(0:(L/2))/L;
P1 = 20.*log10(P1);
semilogx(f,P1) 
axis([20,22000,-20,20])
title('Equalizer Filter')
xlabel('f (Hz)')
ylabel('amplitude')
