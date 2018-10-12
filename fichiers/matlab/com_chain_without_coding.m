%% Projet EN321 - Simulation d'une TX/RX sans modulation OFDM - Version double codeur de canal
% La chaine de communication est compos�e de 2 codeurs de canal s�par�s par
% un entrelaceur convolutif. Le premier codeur est un codeur en bloc, le
% second est un codeur en ligne.

clear all;
close all;

%%%%%%%%%%%
%% INITIALIZATION
%%%%%%%%%%%

Fe=20e6; % sampling frequency
Te=1/Fe;
TC=40; % temperature (Celcius)
TK=274+TC; % temperature (Kelvin)
f0=2e9; % carrier frequency
kboltzman=1.3806400e-23; % Boltzmann constant
N0=kboltzman*Fe*TK; % Noise power
N0dB=10*log10(N0);
sprintf('Noise power : %d dB',N0dB)
Ptx=0.01; % transmitted signal power (Watt)
Ptx_dB=10*log10(Ptx);
d=100*rand+10; % distance between Tx and Rx (meters)
c=3e8; % speed of light
Prx=Ptx*(c/(f0*4*pi*d))^2; % recevied signal power (Watt)
Prx_dB=10*log10(Prx);

sprintf('Distance Tx/Rx (m): %d',d)
sprintf('Transmitted signal power  : %d dB',Ptx_dB)
sprintf('Received signal power  : %d dB',Prx_dB)
sprintf('SNR at the receiver side : %d dB',Prx_dB-N0dB)

NFFT=64;
% bch_k=52; % sub-carrier number

%W=bch_k/NFFT*Fe; % transmitted signal bandwidth

nb=2; % number of bits per symbol
disp('Code MAC binaire correspondant � la modulation num�rique : ')
b0_b2=de2bi(nb,3,'left-msb')

M=2^nb;
type_mod='psk';

L=30; % size of the channel
freq_axis = [-1/(2*Te):1/(NFFT*Te):1/(2*Te)-1/(NFFT*Te)];
noise_variance=N0; % noise variance

%% Data reading/generation

data_mode = 'rand_binary_image'; % generation of a random binary image
%data_mode = 'color_image';

if(data_mode == 'rand_binary_image')
    
    % Generation de donnees aleatoire et disposition dans une image
    Nb_ligne_IMG=100;
    Nb_colonne_IMG=100;
    U_soft_size=Nb_ligne_IMG*Nb_colonne_IMG;   % nombre de bits utiles codés
    % génération aléatoire de donn?es binaires
    rng(654354)
    tmp=(randi(2,U_soft_size)-1);
    U_soft = tmp(1,:)
    % on place les données dans une matrice qui sera affichée comme une image
    img2send=reshape(U_soft,Nb_ligne_IMG,Nb_colonne_IMG)
    
elseif(data_mode == 'color_image')
    
    % Lecture d'une image
    img2send=imread('./bdd_image/logo.jpg'); % l'image est retourn?e sous la forme d'une matrice 3D RGB
    U_soft_R=reshape(de2bi(reshape(img2send(:,:,1),[],1),8,'left-msb').',[],1); % flux binaire du rouge
    U_soft_G=reshape(de2bi(reshape(img2send(:,:,2),[],1),8,'left-msb').',[],1); % flux binaire du vert
    U_soft_B=reshape(de2bi(reshape(img2send(:,:,3),[],1),8,'left-msb').',[],1); % flux binaire du bleu    
    U_soft=[U_soft_R;U_soft_G;U_soft_B].';
    U_soft_size=length(U_soft);
    Nb_ligne_IMG=size(img2send,1);
    Nb_colonne_IMG=size(img2send,2);
    U_soft=[U_soft_R;U_soft_G;U_soft_B].';
end

U_soft_size=length(U_soft);

%%%--------------------------------------------------------------------%%%%
%%- DIGITAL MODULATION
%%%---------------------------------------------------------------------%%%
X=bi2de(reshape(U_soft.',length(U_soft)/nb,nb),'left-msb').'; % bit de poids fort � gauche
init_phase=0;
if type_mod=='psk'
    if nb==2
        init_phase=pi/4;
    end
       symb_utiles = pskmod(X,M,init_phase,'gray');
elseif type_mod=='qam'
       symb_utiles = qammod(X,M,0,'gray');
else
    sprintf('Erreur modulation inconnue')
    s=[];
end


%%
%%%--------------------------------------------------------------------%%%%
%%- CHANNEL (normalized channel : average power)
%%%---------------------------------------------------------------------%%%
h = 1; % discrete channel without multi-path
% h=sqrt(1/(2*L))*(randn(1,L)+1i*randn(1,L)); % discrete channel with multi-path
y = filter(h,1,symb_utiles);
       
%%
%%%--------------------------------------------------------------------%%%%
%% RECEIVER
%%%---------------------------------------------------------------------%%%
%noise_variance = 0.2; noisy channel
noise_variance = 0; % noise-less channel

noise = sqrt(noise_variance/(2))*(randn(size(y))+1i*randn(size(y)));
z = y + noise;

%% OFDM Demodulator 
% No OFDM here

%% Channel equalizer
% No channel equalization

%% Demodulation

symb_U_Rx = z;

init_phase = 0;
if type_mod=='psk'
    if nb==2
        init_phase=pi/4;
    end
       s = pskdemod(symb_U_Rx,M,init_phase,'gray');
       X=de2bi(s,log2(M),'left-msb').'; % bit de poids fort � gauche   
       
else
       s = qamdemod(symb_U_Rx,M,0,'gray');
       X=de2bi(s,log2(M),'left-msb').'; % bit de poids fort � gauche
       
end

U_r_soft=reshape(X.',1,[]);



BER_U = mean(abs(U_r_soft-U_soft));

%% Image reconstruction

if(data_mode == 'rand_binary_image')
    imgRx=reshape(U_r_soft,Nb_ligne_IMG,Nb_colonne_IMG);
elseif(data_mode == 'color_image')
    bitsRx=reshape(U_r_soft,[],3);
    intRx_R=uint8(bi2de(reshape(bitsRx(:,1),8,[]).','left-msb'));
    intRx_G=uint8(bi2de(reshape(bitsRx(:,2),8,[]).','left-msb'));
    intRx_B=uint8(bi2de(reshape(bitsRx(:,3),8,[]).','left-msb'));

    imgRx(:,:,1)=reshape(intRx_R,Nb_ligne_IMG,Nb_colonne_IMG);
    imgRx(:,:,2)=reshape(intRx_G,Nb_ligne_IMG,Nb_colonne_IMG);
    imgRx(:,:,3)=reshape(intRx_B,Nb_ligne_IMG,Nb_colonne_IMG);
end

figure(5)
subplot 131;
if(data_mode == 'rand_binary_image')
    imagesc(img2send)
elseif(data_mode == 'color_image')
    image(img2send)
end
title('Image emise')

subplot 132;
if(data_mode == 'rand_binary_image')
    imagesc(imgRx)
elseif(data_mode == 'color_image')
    image(imgRx)
end
title('Image recue')

subplot 133;
if(data_mode == 'rand_binary_image')
    imagesc(img2send-imgRx)
elseif(data_mode == 'color_image')
    image(img2send-imgRx)
end
title('diff des images')


%% BER results
disp('--------------------------------------------------------------------')
disp(sprintf('SNR at the receiver side : %d dB',round(Prx_dB-N0dB)))
disp('--------------------------------------------------------------------')

disp(sprintf('BER after BCH : %d',(BER_U)))
disp('--------------------------------------------------------------------')
