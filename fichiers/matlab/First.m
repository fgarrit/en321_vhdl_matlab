%% Projet EN321 - Simulation d'une TX/RX sans modulation OFDM - Version double codeur de canal
% La chaine de communication est compos�e de 2 codeurs de canal s�par�s par
% un entrelaceur convolutif. Le premier codeur est un codeur en bloc, le
% second est un codeur en ligne.

clear all;
close all;
clc;

%Initialisation

nb=2;
NbSymbOFDM=1000;
N=256;
L=4;
h=[randn(1,L)+1i*randn(1,L)];
H=fft(h,N);
sigmab2=0;
SNR = 30;

%TX
bits=rand(1,nb*N*NbSymbOFDM)>0.5;
for i=1:N*NbSymbOFDM
    if bits((i-1)*nb+1:i*nb)==[0 0]
        symb(1,i)=exp(1i*pi/4);
    elseif bits((i-1)*nb+1:i*nb)==[0 1]
        symb(1,i)=exp(1i*3*pi/4);
    elseif bits((i-1)*nb+1:i*nb)==[1 1]
        symb(1,i)=exp(1i*5*pi/4);
    else
        symb(1,i)=exp(1i*7*pi/4);
    end
end
        
%scatterplot(symb)

%% Modulation OFDM

%matrix = reshape(symb,[N,NbSymbOFDM]);
trame_OFDM=[];
IG = zeros(1,L-1);

for i = 1:NbSymbOFDM
    symb_OFDM = ifft(symb(1,(i-1)*N+1:i*N),N);
    trame_OFDM=[trame_OFDM symb_OFDM(1,end-length(IG)+1:end) symb_OFDM];
end



%% Canal
sortie = filter(h,1,trame_OFDM);
Pin = 0;
Pin = mean(abs(sortie).^2);
sigmab2 = Pin/10^(SNR/10);

%% Rx
bruit = sqrt(sigmab2/2)*(randn(size(sortie))+1i*randn(size(sortie)));
y = sortie + bruit;


symb_estime=[];
for i=1:NbSymbOFDM
    temp=fft(y(1,i*length(IG)+(i-1)*N+1:i*N+i*length(IG)),N)./H;
    symb_estime=[symb_estime temp];
end


for i=1:length(symb_estime)
    if real(symb_estime(1,i))>0
        if imag(symb_estime(1,i))>0
            bits_estime(1,(i-1)*nb+1:i*nb)=[0 0];
        else
            bits_estime(1,(i-1)*nb+1:i*nb)=[1 0];
        end
    else
        if imag(symb_estime(1,i))>0
            bits_estime(1,(i-1)*nb+1:i*nb)=[0 1];
        else
            bits_estime(1,(i-1)*nb+1:i*nb)=[1 1];
        end
    end
end

figure(2)
plot(real(trame_OFDM));
figure(3)
grid on
plot(real(symb_estime),imag(symb_estime),'ok')
hold on
plot(real(symb),imag(symb),'*r')

BER = mean(abs(bits-bits_estime))



