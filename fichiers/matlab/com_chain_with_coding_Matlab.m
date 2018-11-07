%% Projet EN321 - Simulation d'une TX/RX sans modulation OFDM - Version double codeur de canal
% La chaine de communication est compos�e de 2 codeurs de canal s�par�s par
% un entrelaceur convolutif. Le premier codeur est un codeur en bloc, le
% second est un codeur en ligne.

clear all;
close all;

instrreset;

%%%%%%%%%%%
%% INITIALIZATION
%%%%%%%%%%%

TB_file_ID = fopen('U.txt','w');

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

nb=2    ; % number of bits per symbol
nb_mac=2
disp('Code MAC binaire correspondant � la modulation num�rique : ')
b0_b2=de2bi(nb,3,'left-msb')

M=2^nb;
M_mac=2^nb_mac;
type_mod='qam';

L=30; % size of the channel
freq_axis = [-1/(2*Te):1/(NFFT*Te):1/(2*Te)-1/(NFFT*Te)];
noise_variance=N0; % noise variance


%% Scrambler parameters
scramb_polynomial=[1 1 1 0 1];
scramb_init_state=[0 0 0 0];
Scrambler_U_obj=comm.Scrambler(2,scramb_polynomial,scramb_init_state); % scrambler creation
Scrambler_U_obj_mac=comm.Scrambler(2,scramb_polynomial,scramb_init_state);

%% BCH parameters
bch_n=7;   % code block-length
bch_k=4;   % code dimension

%% Convolutionnal interleaver
intlvr_line_nb=7; % nb of lines ( = shift-registers) in the interleaver
intlvr_reg_size=1; % number of bits per register

%% Convolutionnal encoder parameters 
trellis = poly2trellis(3,[5 3]); % generator polynomial : (15,13)


%% Data reading/generation

data_mode = 1;   %'rand_binary_image' generation of a random binary image
%data_mode = 0;  %'color_image'

if(data_mode == 1)
    
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
    
elseif(data_mode == 0)
    
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
%%- CHANNEL CODING
%%%---------------------------------------------------------------------%%%


%% Padding for the BCH encoder and the interleaver

full_bch_cwd_nb = floor(U_soft_size/bch_k)
bch_cwd_nb = (full_bch_cwd_nb +1) + (intlvr_line_nb-1);

intlvr_pad_bit_nb = bch_k * intlvr_reg_size * (intlvr_line_nb - 1); % after BCH encoding, there is intlvr_line_nb*intlvr_reg_size*(intlvr_line_nb-1) padding bits for the interleaver
bch_pad_bit_nb = bch_k-(U_soft_size-full_bch_cwd_nb*bch_k);

total_pad_bit_nb = bch_pad_bit_nb + intlvr_pad_bit_nb; 

sprintf('Nb of padding bits for the last BCH codeword + interleaver: %d',total_pad_bit_nb)

padding_bits=rand(1,total_pad_bit_nb) > 0.5;
bch_bit_nb = bch_cwd_nb * bch_n;

V_soft = [U_soft, padding_bits];
V_soft_size = length(V_soft);

%% Calcul padding SOFT
padding = 0;

while(mod(length(V_soft)*bch_n/bch_k*2 + padding,NFFT) ~= 0 || (mod(length(V_soft)*bch_n/bch_k*2 + padding, NFFT*nb) ~= 0))
    %bch_cwd_nb = bch_cwd_nb +1;
    padding = padding + 1;
end

%% [AJOUT] MAC
if type_mod=='psk'
    type_mod_bin = 1;
elseif type_mod=='qam'
    type_mod_bin = 0;
else
    sprintf('Erreur modulation inconnue')
    s=[];
end

nb_bin = de2bi(nb,3,'left-msb');
padd_bin = de2bi(padding,11,'left-msb');
intlvr_pad_bit_nb_bin = de2bi(intlvr_pad_bit_nb,10,'left-msb');


V_mac = [type_mod_bin b0_b2 padd_bin intlvr_pad_bit_nb_bin]
V_mac_size = length(V_mac);

%%Padding MAC
full_bch_cwd_nb_mac = floor(V_mac_size/bch_k)
bch_cwd_nb_mac = (full_bch_cwd_nb_mac +1) + (intlvr_line_nb-1);

intlvr_pad_bit_nb_mac = bch_k * intlvr_reg_size * (intlvr_line_nb - 1); % after BCH encoding, there is intlvr_line_nb*intlvr_reg_size*(intlvr_line_nb-1) padding bits for the interleaver
bch_pad_bit_nb_mac = bch_k-(V_mac_size-full_bch_cwd_nb_mac*bch_k);

total_pad_bit_nb_mac = bch_pad_bit_nb_mac + intlvr_pad_bit_nb_mac; 

sprintf('Nb of padding bits for MAC for the last BCH codeword + interleaver: %d',total_pad_bit_nb_mac)

padding_bits_mac=rand(1,total_pad_bit_nb_mac) > 0.5 ;
bch_bit_nb_mac = bch_cwd_nb_mac * bch_n;

V_mac_padd = [V_mac padding_bits_mac];
V_mac_padd_size = length(V_mac_padd);

%% Calcul padding MAC
padding_mac = 0;

while(mod((length(V_mac_padd)*bch_n/bch_k)*2 + padding_mac,NFFT) ~= 0 || (mod((length(V_mac_padd)*bch_n/bch_k)*2 + padding_mac, NFFT*2) ~= 0))
    padding_mac = padding_mac + 1;
end

% %% Write UART
% s = send_UART(V_mac_padd,V_mac_padd_size)

%% Scrambler MAC
S_mac=step(Scrambler_U_obj_mac,V_mac_padd.');

% BCH Encoder MAC
X_gf_mac = bchenc(gf(reshape(S_mac, bch_k, bch_cwd_nb_mac).',1), bch_n, bch_k); % codeur BCH(bch_n,bch_k)
X_mac = double( X_gf_mac.x );


% Interleaver MAC
P_mac=convintrlv([reshape(X_mac.',1,[])],intlvr_line_nb,intlvr_reg_size);

% Convolutionnal Encoder
C_mac = convenc(P_mac,trellis); 


% % Read UART
% C_hard_mac = recv_UART(s, bch_bit_nb_mac);
% C_hard_mac = reshape(de2bi(C_hard_mac)',1,[]);
% C_mac= C_hard_mac;

% Padding
padding_add_mac = zeros(1,padding_mac);
C_mac_padd = [C_mac padding_add_mac];

%%%--------------------------------------------------------------------%%%%
%%- DIGITAL MODULATION MAC
%%%---------------------------------------------------------------------%%%
XX=bi2de(reshape(C_mac_padd.',length(C_mac_padd)/nb_mac,nb_mac),'left-msb').'; % bit de poids fort � gauche
init_phase=0;
symb_utiles_mac = qammod(XX,M_mac,0,'gray');

% OFDM Modulator 
trame_OFDM_mac=[];
LL=4;
IG = zeros(1,LL-1);

for i = 1:(length(symb_utiles_mac)/NFFT)
    symb_OFDM_mac = ifft(symb_utiles_mac(1,(i-1)*NFFT+1:i*NFFT),NFFT);
    trame_OFDM_mac=[trame_OFDM_mac symb_OFDM_mac];
    %trame_OFDM=[trame_OFDM symb_OFDM(1,end-length(IG)+1:end) symb_OFDM];
end




% %% Write UART
% s = send_UART(V_soft,V_soft_size)

%% Scrambler
S_soft=step(Scrambler_U_obj,V_soft.');

%% BCH Encoder
X_gf_soft = bchenc(gf(reshape(S_soft, bch_k, bch_cwd_nb).',1), bch_n, bch_k); % codeur BCH(bch_n,bch_k)
X_soft = double( X_gf_soft.x );

%% Interleaver
P_soft=convintrlv([reshape(X_soft.',1,[])],intlvr_line_nb,intlvr_reg_size);

%% Convolutionnal Encoder
C_soft = convenc(P_soft,trellis);

% % Read UART
% C_hard = recv_UART(s, bch_bit_nb);
% C_hard = reshape(de2bi(C_hard)',1,[]);
% C_soft= C_hard;

%%Padding
padding_add = rand(1,padding) > 0.5;
C_soft_padd = [C_soft padding_add];





%%%--------------------------------------------------------------------%%%%
%%- DIGITAL MODULATION
%%%---------------------------------------------------------------------%%%
X=bi2de(reshape(C_soft_padd.',length(C_soft_padd)/nb,nb),'left-msb').'; % bit de poids fort � gauche
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

%% OFDM Modulator 
trame_OFDM=[];
LL=40;
IG = zeros(1,LL-1);

for i = 1:(length(symb_utiles)/NFFT)
    symb_OFDM = ifft(symb_utiles(1,(i-1)*NFFT+1:i*NFFT),NFFT);
    trame_OFDM=[trame_OFDM symb_OFDM];
    %trame_OFDM=[trame_OFDM symb_OFDM(1,end-length(IG)+1:end) symb_OFDM];
end

%% Concatenation MAC + trame
trame_OFDM_final = [trame_OFDM_mac trame_OFDM];

%%
%%%--------------------------------------------------------------------%%%%
%%- CHANNEL (normalized channel : average power)
%%%---------------------------------------------------------------------%%%
h = 1; % discrete channel without multi-path
%h=sqrt(1/(2*L))*(randn(1,L)+1i*randn(1,L)); % discrete channel with multi-path
H=fft(h,NFFT);
y = filter(h,1,trame_OFDM_final);
       
%%
%%%--------------------------------------------------------------------%%%%
%% RECEIVER
%%%---------------------------------------------------------------------%%%
noise_variance = 0.001 %0.2
noise = sqrt(noise_variance/(2))*(randn(size(y))+1i*randn(size(y)));
z = y + noise; 
%z=y;

%%Separation MAC/Donnees
z_mac = z(1:128);
z_donnees = z(129:length(z));

%% OFDM Demodulator MAC
symb_estime_mac=[];
for i=1:(length(z_mac)/NFFT)
    temp_mac=fft(z_mac(1,(i-1)*NFFT+1:i*NFFT),NFFT);
    %temp_mac=fft(z_mac(1,(i-1)*NFFT+1:i*NFFT),NFFT)./H;
    %temp_mac=fft(z_mac(1,i*length(IG)+(i-1)*NFFT+1:i*NFFT+i*length(IG)),NFFT);
    %temp_mac=fft(z_mac(1,i*length(IG)+(i-1)*NFFT+1:i*NFFT+i*length(IG)),NFFT)./H;
    symb_estime_mac=[symb_estime_mac temp_mac];
end


%% Demodulation

symb_U_Rx_mac = symb_estime_mac;

init_phase = 0;
s_mac = qamdemod(symb_U_Rx_mac,M_mac,0,'gray');
X_r_mac=de2bi(s_mac,log2(M_mac),'left-msb').'; % bit de poids fort � gauche
C_r_pad_mac=reshape(X_r_mac.',1,[]);

%Deppading 
C_r_mac = C_r_pad_mac(1:length(C_r_pad_mac)-length(padding_add_mac));


%% Viterbi Decoding

trellis_depth=42; % profondeur du trellis

P_r_mac = vitdec(C_r_mac,trellis,trellis_depth,'trunc','hard');

BER_U_A_Viterbi_mac = mean(abs(P_mac-P_r_mac))


%% Deinterleaving

X_r_mac=convdeintrlv(P_r_mac,intlvr_line_nb,intlvr_reg_size);

%% BCH decoding

S_r_mac_gf=bchdec(gf(reshape(X_r_mac,bch_n,bch_cwd_nb_mac).',1),bch_n,bch_k); 
S_r_mac = uint8(S_r_mac_gf.x);

S_r_mac_Depad_temp = reshape(S_r_mac.',1,[]);
S_r_mac_Depad = S_r_mac_Depad_temp(intlvr_pad_bit_nb+1:end);
S_r_mac_Depad = S_r_mac_Depad(1:end-bch_pad_bit_nb_mac);

%% Descrambler

Descrambler_U_obj = comm.Descrambler(2,scramb_polynomial,scramb_init_state);

S_r_mac_Depad=step(Descrambler_U_obj,S_r_mac_Depad.'); % descrambler

%% Depadding
S_r_mac_Depad_final = fliplr(S_r_mac_Depad(1:25));

BER_U_MAC = mean(abs(S_r_mac_Depad_final-uint8(V_mac')));

disp('--------------------------------------------------------------------')

fprintf('BER FINAL MAC : %d \n',BER_U_MAC)

disp('--------------------------------------------------------------------')

%% Decode MAC

type_mod_bin_mac = S_r_mac_Depad_final(1);
b0_b2_mac = bi2de(S_r_mac_Depad_final(2:4)','left-msb');
padding_extrait_mac = double(bi2de(S_r_mac_Depad_final(5:15)','left-msb'));
intlvr_pad_bit_nb_bin_mac = double(bi2de(S_r_mac_Depad_final(16:25)','left-msb'));

if type_mod_bin_mac == 1 
    type_mod_mac='psk';
elseif type_mod_bin_mac == 0
    type_mod_mac='qam';
else
    sprintf('Erreur modulation inconnue')
    s=[];
end




%% OFDM Demodulator 
symb_estime=[];
for i=1:(length(z_donnees)/NFFT)
    temp=fft(z_donnees(1,(i-1)*NFFT+1:i*NFFT),NFFT);
    %temp=fft(z_donnees(1,(i-1)*NFFT+1:i*NFFT),NFFT)./H;
    %temp=fft(y(1,i*length(IG)+(i-1)*NFFT+1:i*NFFT+i*length(IG)),NFFT);
    %temp=fft(y(1,i*length(IG)+(i-1)*NFFT+1:i*NFFT+i*length(IG)),NFFT)./H;
    symb_estime=[symb_estime temp];
end

figure(1)
hold on;
plot(symb_estime,'ok');
plot(symb_utiles,'*r');


%% Channel equalizer
% No channel equalization

%% Demodulation

symb_U_Rx = symb_estime;

init_phase = 0;
if type_mod_mac=='psk'
    if nb==2
        init_phase=pi/4;
    end
       s = pskdemod(symb_U_Rx,M,init_phase,'gray');
       X=de2bi(s,log2(M),'left-msb').'; % bit de poids fort � gauche   
       
else
       s = qamdemod(symb_U_Rx,M,0,'gray');
       X=de2bi(s,log2(M),'left-msb').'; % bit de poids fort � gauche
       
end

C_r_soft_pad=reshape(X.',1,[]);

%%Deppading 

C_r_soft = C_r_soft_pad(1:length(C_r_soft_pad)-padding_extrait_mac);
%C_r_soft = C_r_soft_pad(1:length(C_r_soft_pad)-102);

%% Viterbi Decoding

trellis_depth=42; % profondeur du trellis

P_r_soft = vitdec(C_r_soft,trellis,trellis_depth,'trunc','hard');

BER_U_A_Viterbi = mean(abs(P_soft-P_r_soft))


%% Deinterleaving

X_r_soft=convdeintrlv(P_r_soft,intlvr_line_nb,intlvr_reg_size);

%% BCH decoding

S_r_soft_gf=bchdec(gf(reshape(X_r_soft,bch_n,bch_cwd_nb).',1),bch_n,bch_k); 
S_r_soft = uint8(S_r_soft_gf.x);

S_r_soft_Depad_temp = reshape(S_r_soft.',1,[]);
S_r_soft_Depad = S_r_soft_Depad_temp(intlvr_pad_bit_nb_bin_mac+1:end);
S_r_soft_Depad = S_r_soft_Depad(1:end-bch_pad_bit_nb);


%% Descrambler

Descrambler_U_obj = comm.Descrambler(2,scramb_polynomial,scramb_init_state);

S_r_soft_Depad=step(Descrambler_U_obj,S_r_soft_Depad.'); % descrambler


BER_U = mean(abs(S_r_soft_Depad-uint8(U_soft')))



%% Image reconstruction

if(data_mode == 1)
    imgRx=reshape(S_r_soft_Depad,Nb_ligne_IMG,Nb_colonne_IMG);
elseif(data_mode == 0)
    bitsRx=reshape(S_r_soft_Depad,[],3);
    intRx_R=uint8(bi2de(reshape(bitsRx(:,1),8,[]).','left-msb'));
    intRx_G=uint8(bi2de(reshape(bitsRx(:,2),8,[]).','left-msb'));
    intRx_B=uint8(bi2de(reshape(bitsRx(:,3),8,[]).','left-msb'));

    imgRx(:,:,1)=reshape(intRx_R,Nb_ligne_IMG,Nb_colonne_IMG);
    imgRx(:,:,2)=reshape(intRx_G,Nb_ligne_IMG,Nb_colonne_IMG);
    imgRx(:,:,3)=reshape(intRx_B,Nb_ligne_IMG,Nb_colonne_IMG);
end

figure(5)
subplot 131;
if(data_mode == 1)
    imagesc(img2send)
elseif(data_mode == 0)
    image(img2send)
end
title('Image emise')

subplot 132;
if(data_mode == 1)
    imagesc(imgRx)
elseif(data_mode == 0)
    image(imgRx)
end
title('Image recue')
subplot 133;
if(data_mode == 1)
    imagesc(uint8(img2send)-imgRx)
elseif(data_mode == 0)
    image(uint8(img2send)-imgRx)
end
title('diff des images')


%% BER results
disp('--------------------------------------------------------------------')
fprintf('SNR at the receiver side : %d dB\n',round(Prx_dB-N0dB))
disp('--------------------------------------------------------------------')

fprintf('BER after Viterbi decoding: %d\n',(BER_U_A_Viterbi))
fprintf('BER after BCH : %d\n',(BER_U))
disp('--------------------------------------------------------------------')

