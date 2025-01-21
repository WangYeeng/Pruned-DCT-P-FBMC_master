%==========================================================================
% Author: Ying Wang
% Function: BER Simulation for MIMO Transmission in Pruned DCT-p-FBMC System.
% NOTE: This code simulates both "Flat Channel" and " Doubly Selective Channel".
% E-mail address: wangyingstu@163.com
%==========================================================================
clear; clc; close all;
addpath("Channel\","Coding\","MultiCarrierModulation\");rng(0);
%% Parameters
NrRepetitions       = 3000;                   % Monte Carlo Repetitions 
M_SNR_dB_OFDM       = -10:5:30;               % SNR (dB)
Channel_Infor       = 'Flat';                 % "Flat" simulates a "flat channel"/"Doub" simulates a "doubly selective channel".
% FBMC and OFDM Parameters Configuration
QAM_ModulationOrder = 4;                      % Modulation order, 4, 16, 64,...
SubcarrierSpacing   = 15e3;                   % Subcarrier spacing (15kHz)
CarrierFrequency    = 2.5e9;
L                   = 16;                     % Number of subcarriers
NrBlocks            = 2;                      % Must be even for Alamouti OFDM. 
K_FBMC              = 30;
K_OFDM              = 15;
NrGuardSubcarriers  = 1; 
CP_Length           = 1/SubcarrierSpacing/14; % LTE CP Length (seconds)
SamplingRate        = 15e3*L*1*21;  
% Channel parameters
switch Channel_Infor
    case 'Flat'% flat channel
        Velocity_kmh        = 0;
        PowerDelayProfile   = 'Flat';          
    case 'Doub'% doubly selective channel 
        Velocity_kmh        = 500;
        PowerDelayProfile   = 'TDL-A_120ns';    
end

%% FBMC
FBMC = FBMC(...
    (L+NrGuardSubcarriers)*NrBlocks,...   % Number of subcarriers
    K_FBMC,...                            % FBMC No. of symbols
    SubcarrierSpacing,...                 % Subcarrier spacing (Hz)
    SamplingRate,...                      % Sampling rate (Samples/s)
    0,...                                 % Intermediate frequency of the first subcarrier (Hz).
    false,...                             % Transmit real valued signal
    'PHYDYAS-OQAM',...                    % Prototype filter and OQAM or QAM. 
    4, ...                                % Overlapping factor
    pi, ...                               % Initial phase shift
    true ...                              % Efficient IFFT implementation, true or false
    );

%% OFDM
ZeroGuardTimeLength = ((FBMC.Nr.SamplesTotal-(round(SamplingRate/SubcarrierSpacing)+CP_Length*SamplingRate)*K_OFDM)/2)/SamplingRate;
OFDM = OFDM(...
    (L+NrGuardSubcarriers)*NrBlocks,...   % Number of subcarriers
    K_OFDM,...                            % OFDM No. of symbols
    SubcarrierSpacing,...                 % Subcarrier spacing (Hz)
    SamplingRate,...                      % Sampling rate (Samples/s)
    0,...                                 % Intermediate frequency of the first subcarrier (Hz).
    false,...                             % Transmit real valued signal
    CP_Length, ...                        % Cycle prefix length (s) 1/SubcarrierSpacing/(K/2-1)
    ZeroGuardTimeLength ...               % Length of the guard time (s)
    );

%% Check the number of samples
if  OFDM.Nr.SamplesTotal~=FBMC.Nr.SamplesTotal
   error('Total number of samples must be the same for OFDM and FBMC.');
end
N = OFDM.Nr.SamplesTotal;

%% Alamouti
Alamouti = SpaceCoding(...
    'Alamouti2x1',...                     % Coding method
    0 ...                                 
    );

%% Symbol constellations
QAM = SignalConstellation(QAM_ModulationOrder,'QAM');

%% Channel
ChannelModel = FastFading(...
                SamplingRate,...                               % Sampling rate (Samples/s)
                PowerDelayProfile,...                          % Power delay profile
                FBMC.Nr.SamplesTotal,...                       % Number of total samples
                Velocity_kmh/3.6*CarrierFrequency/2.998e8,...  % Maximum Doppler shift: Velocity_kmh/3.6*CarrierFrequency/2.998e8
                'Jakes',...                                    % Which Doppler model: 'Jakes', 'Uniform', 'Discrete-Jakes', 'Discrete-Uniform'. For "Discrete-", we assume a discrete Doppler spectrum to improve the simulation time. This only works accuratly if the number of samples and the velocity is sufficiently large
                200,...                                        % Number of paths for the WSSUS process
                2,...                                          % Number of transmit antennas
                2,...                                          % Number of receive antennas
                0 ...                                          % Gives a warning if the predefined delay taps of the channel do not fit the sampling rate. This is usually not much of a problem if they are approximatly the same.
                );

%% Calculation
% ML Detection
ML_MapIndex1 = reshape(repmat((1:QAM.ModulationOrder),QAM.ModulationOrder,1),1,QAM.ModulationOrder^2);
ML_MapIndex2 = reshape(repmat((1:QAM.ModulationOrder).',1,QAM.ModulationOrder),1,QAM.ModulationOrder^2);
ML_Mapping = QAM.SymbolMapping([ML_MapIndex1;ML_MapIndex2]);

% Total number of data symbols
NrDataSymbols_OFDM = OFDM.Nr.MCSymbols*OFDM.Nr.Subcarriers;
NrDataSymbols_FBMC = L/2*NrBlocks*FBMC.Nr.MCSymbols;

% Obtain the transmitted and received matrices
TXMatrix_OFDM = OFDM.GetTXMatrix;
RXMatrix_OFDM = OFDM.GetRXMatrix;
TXMatrix_FBMC = FBMC.GetTXMatrix;
RXMatrix_FBMC = FBMC.GetRXMatrix;

% Reduce storage with sparse matrices
TXMatrix_OFDM(abs(TXMatrix_OFDM)<10^-10)=0;
TXMatrix_OFDM = sparse(TXMatrix_OFDM);
RXMatrix_OFDM(abs(RXMatrix_OFDM)<10^-10)=0;
RXMatrix_OFDM = sparse(RXMatrix_OFDM);

TXMatrix_FBMC(abs(TXMatrix_FBMC)<10^-10)=0;
TXMatrix_FBMC = sparse(TXMatrix_FBMC);
RXMatrix_FBMC(abs(RXMatrix_FBMC)<10^-10)=0;
RXMatrix_FBMC = sparse(RXMatrix_FBMC);

%% Precoding matrix
NrSubcarriers = FBMC.Nr.Subcarriers;
NrSymbols     = FBMC.Nr.MCSymbols;
Domain        = 'Time';     % Setup Domain
Mathod        = 'DCT';      % Adopting the DCT method
Kall          = K_FBMC;     % No. of symbols
K             = Kall;       % FBMC No. of symbols
Savdecision   = 0;          % Whether to save data
phi = GetCodingMatrix(...
    Domain,...          % Determine the domain ('Frequency', 'Time', 'DeEigen')
    K,...               % FBMC No. of symbols
    L,...               % Number of subcarriers
    Mathod,...          % Precoding matrix form
    FBMC...             % The functions are embedded in the FBMC system.
    );

FBMC.SetNrSubcarriers(NrSubcarriers);
FBMC.SetNrMCSymbols(NrSymbols);
% Mapping precoding Matrix to large Matrix describing the whole system.
Phi_OneFBMCsymbol = kron(sparse(eye(NrBlocks)),[phi;zeros(1,size(phi,2))]); 
Phi_Totak = kron(eye(K_FBMC),Phi_OneFBMCsymbol);
% Power normalization ensures that OFDM and FBMC are identical.
FBMCPowerNormalization = sqrt(2*(L+NrGuardSubcarriers)/L);
MeanChannelMatrix = abs(Phi_Totak').^2;

%% Pre-configuration
BER_OFDM_Alamouti = nan(length(M_SNR_dB_OFDM),NrRepetitions);
BER_FBMC_Alamouti = nan(length(M_SNR_dB_OFDM),NrRepetitions);
BER_OFDM_SM_ZeroForcing = nan(length(M_SNR_dB_OFDM),NrRepetitions);
BER_FBMC_SM_ZeroForcing = nan(length(M_SNR_dB_OFDM),NrRepetitions);
BER_OFDM_SM_ML = nan(length(M_SNR_dB_OFDM),NrRepetitions);
BER_FBMC_SM_ML = nan(length(M_SNR_dB_OFDM),NrRepetitions);

%% Start Simulation
tic;
for i_rep = 1:NrRepetitions   
%% Binary data streams
BinaryDataStream_Alamouti_OFDM = randi([0 1],NrDataSymbols_OFDM*log2(QAM.ModulationOrder),1);
BinaryDataStream_Alamouti_FBMC = randi([0 1],NrDataSymbols_FBMC*log2(QAM.ModulationOrder),1);

% SM = Spatial Multiplexing
BinaryDataStream_SMAntenna1_OFDM = randi([0 1],NrDataSymbols_OFDM*log2(QAM.ModulationOrder),1); %Spatial Multiplexing
BinaryDataStream_SMAntenna2_OFDM = randi([0 1],NrDataSymbols_OFDM*log2(QAM.ModulationOrder),1); %Spatial Multiplexing
BinaryDataStream_SMAntenna1_FBMC = randi([0 1],NrDataSymbols_FBMC*log2(QAM.ModulationOrder),1); %Spatial Multiplexing
BinaryDataStream_SMAntenna2_FBMC = randi([0 1],NrDataSymbols_FBMC*log2(QAM.ModulationOrder),1); %Spatial Multiplexing

BinaryDataStream_SM_OFDM = [BinaryDataStream_SMAntenna1_OFDM;BinaryDataStream_SMAntenna2_OFDM];
BinaryDataStream_SM_FBMC = [BinaryDataStream_SMAntenna1_FBMC;BinaryDataStream_SMAntenna2_FBMC];


%% Transmitted symbols
% Alamouti
x_Alamouti_OFDM = reshape(QAM.Bit2Symbol(BinaryDataStream_Alamouti_OFDM),OFDM.Nr.Subcarriers,OFDM.Nr.MCSymbols);
x_Alamouti_Coded_OFDM = Alamouti.Encoder(x_Alamouti_OFDM);
x_Alamouti_Coded_Antenna1_OFDM = x_Alamouti_Coded_OFDM(:,:,1);
x_Alamouti_Coded_Antenna2_OFDM = x_Alamouti_Coded_OFDM(:,:,2);

x_Alamouti_FBMC = reshape(QAM.Bit2Symbol(BinaryDataStream_Alamouti_FBMC),L/2*NrBlocks,FBMC.Nr.MCSymbols);
x_Alamouti_Coded_FBMC = Alamouti.Encoder(x_Alamouti_FBMC);
x_Alamouti_Coded_Antenna1_FBMC = x_Alamouti_Coded_FBMC(:,:,1);
x_Alamouti_Coded_Antenna2_FBMC = x_Alamouti_Coded_FBMC(:,:,2);

% Spatial Multiplexing
x_SM_Antenna1_OFDM =  reshape(QAM.Bit2Symbol(BinaryDataStream_SMAntenna1_OFDM),OFDM.Nr.Subcarriers,OFDM.Nr.MCSymbols);
x_SM_Antenna2_OFDM =  reshape(QAM.Bit2Symbol(BinaryDataStream_SMAntenna2_OFDM),OFDM.Nr.Subcarriers,OFDM.Nr.MCSymbols);

x_SM_Antenna1_FBMC =  reshape(QAM.Bit2Symbol(BinaryDataStream_SMAntenna1_FBMC),L/2*NrBlocks,FBMC.Nr.MCSymbols);
x_SM_Antenna2_FBMC =  reshape(QAM.Bit2Symbol(BinaryDataStream_SMAntenna2_FBMC),L/2*NrBlocks,FBMC.Nr.MCSymbols);

%% Time-domain transmitted signals
s_OFDM_Alamouti_Antenna1 = 1/sqrt(2)*OFDM.Modulation(x_Alamouti_Coded_Antenna1_OFDM);
s_OFDM_Alamouti_Antenna2 = 1/sqrt(2)*OFDM.Modulation(x_Alamouti_Coded_Antenna2_OFDM);

TransmittedSymbols_FBMC_Alamouti_Antenna1 = reshape(Phi_Totak*x_Alamouti_Coded_Antenna1_FBMC(:),FBMC.Nr.Subcarriers,FBMC.Nr.MCSymbols);
TransmittedSymbols_FBMC_Alamouti_Antenna2 = reshape(Phi_Totak*x_Alamouti_Coded_Antenna2_FBMC(:),FBMC.Nr.Subcarriers,FBMC.Nr.MCSymbols);
s_FBMC_Alamouti_Antenna1 = FBMCPowerNormalization/sqrt(2)*FBMC.Modulation(TransmittedSymbols_FBMC_Alamouti_Antenna1);
s_FBMC_Alamouti_Antenna2 = FBMCPowerNormalization/sqrt(2)*FBMC.Modulation(TransmittedSymbols_FBMC_Alamouti_Antenna2);

TransmittedSymbols_FBMC_dft_Alamouti_Antenna1 = reshape(Phi_Totak*x_Alamouti_Coded_Antenna1_FBMC(:),FBMC.Nr.Subcarriers,FBMC.Nr.MCSymbols);
TransmittedSymbols_FBMC_dft_Alamouti_Antenna2 = reshape(Phi_Totak*x_Alamouti_Coded_Antenna2_FBMC(:),FBMC.Nr.Subcarriers,FBMC.Nr.MCSymbols);
s_FBMC_DFT_Alamouti_Antenna1 = FBMC.Modulation(TransmittedSymbols_FBMC_dft_Alamouti_Antenna1);
s_FBMC_DFT_Alamouti_Antenna2 = FBMC.Modulation(TransmittedSymbols_FBMC_dft_Alamouti_Antenna2);

s_OFDM_SM_Antenna1 = 1/sqrt(2)*OFDM.Modulation(x_SM_Antenna1_OFDM);
s_OFDM_SM_Antenna2 = 1/sqrt(2)*OFDM.Modulation(x_SM_Antenna2_OFDM);

TransmittedSymbols_FBMC_SM_Antenna1 = reshape(Phi_Totak*x_SM_Antenna1_FBMC(:),FBMC.Nr.Subcarriers,FBMC.Nr.MCSymbols);
TransmittedSymbols_FBMC_SM_Antenna2 = reshape(Phi_Totak*x_SM_Antenna2_FBMC(:),FBMC.Nr.Subcarriers,FBMC.Nr.MCSymbols);
s_FBMC_SM_Antenna1 = FBMCPowerNormalization/sqrt(2)*FBMC.Modulation(TransmittedSymbols_FBMC_SM_Antenna1);
s_FBMC_SM_Antenna2 = FBMCPowerNormalization/sqrt(2)*FBMC.Modulation(TransmittedSymbols_FBMC_SM_Antenna2);

H11 = ChannelModel.GetConvolutionMatrix{1,1};
ChannelModel.NewRealization;
H12 = ChannelModel.GetConvolutionMatrix{2,1};
ChannelModel.NewRealization;
H21 = ChannelModel.GetConvolutionMatrix{1,2};
ChannelModel.NewRealization;
H22 = ChannelModel.GetConvolutionMatrix{2,2};
ChannelModel.NewRealization;

%% Received Noise-Free Signals. We use the time-varying convolution matrix
r_OFDM_Alamouti_Antenna1_NoNoise = H11*s_OFDM_Alamouti_Antenna1+H12*s_OFDM_Alamouti_Antenna2;
r_FBMC_Alamouti_Antenna1_NoNoise = H11*s_FBMC_Alamouti_Antenna1+H12*s_FBMC_Alamouti_Antenna2;

r_OFDM_SM_Antenna1_NoNoise = H11*s_OFDM_SM_Antenna1+H12*s_OFDM_SM_Antenna2;
r_OFDM_SM_Antenna2_NoNoise = H21*s_OFDM_SM_Antenna1+H22*s_OFDM_SM_Antenna2;

r_FBMC_SM_Antenna1_NoNoise = H11*s_FBMC_SM_Antenna1+H12*s_FBMC_SM_Antenna2;
r_FBMC_SM_Antenna2_NoNoise = H21*s_FBMC_SM_Antenna1+H22*s_FBMC_SM_Antenna2;

H_OFDM(:,:,1,1) = reshape(full(sum((RXMatrix_OFDM*H11).*TXMatrix_OFDM.',2)/sqrt(2)),OFDM.Nr.Subcarriers,OFDM.Nr.MCSymbols);
H_OFDM(:,:,1,2) = reshape(full(sum((RXMatrix_OFDM*H12).*TXMatrix_OFDM.',2)/sqrt(2)),OFDM.Nr.Subcarriers,OFDM.Nr.MCSymbols);
H_OFDM(:,:,2,1) = reshape(full(sum((RXMatrix_OFDM*H21).*TXMatrix_OFDM.',2)/sqrt(2)),OFDM.Nr.Subcarriers,OFDM.Nr.MCSymbols);
H_OFDM(:,:,2,2) = reshape(full(sum((RXMatrix_OFDM*H22).*TXMatrix_OFDM.',2)/sqrt(2)),OFDM.Nr.Subcarriers,OFDM.Nr.MCSymbols);

H_FBMC(:,:,1,1) = reshape(full(sum((RXMatrix_FBMC*H11).*TXMatrix_FBMC.',2)/sqrt(2)*FBMCPowerNormalization),FBMC.Nr.Subcarriers,FBMC.Nr.MCSymbols);
H_FBMC(:,:,1,2) = reshape(full(sum((RXMatrix_FBMC*H12).*TXMatrix_FBMC.',2)/sqrt(2)*FBMCPowerNormalization),FBMC.Nr.Subcarriers,FBMC.Nr.MCSymbols);
H_FBMC(:,:,2,1) = reshape(full(sum((RXMatrix_FBMC*H21).*TXMatrix_FBMC.',2)/sqrt(2)*FBMCPowerNormalization),FBMC.Nr.Subcarriers,FBMC.Nr.MCSymbols);
H_FBMC(:,:,2,2) = reshape(full(sum((RXMatrix_FBMC*H22).*TXMatrix_FBMC.',2)/sqrt(2)*FBMCPowerNormalization),FBMC.Nr.Subcarriers,FBMC.Nr.MCSymbols);

H_FBMC_PostCoding(:,:,1,1) = reshape(MeanChannelMatrix*reshape(H_FBMC(:,:,1,1),[],1),L/2*NrBlocks,FBMC.Nr.MCSymbols);
H_FBMC_PostCoding(:,:,1,2) = reshape(MeanChannelMatrix*reshape(H_FBMC(:,:,1,2),[],1),L/2*NrBlocks,FBMC.Nr.MCSymbols);
H_FBMC_PostCoding(:,:,2,1) = reshape(MeanChannelMatrix*reshape(H_FBMC(:,:,2,1),[],1),L/2*NrBlocks,FBMC.Nr.MCSymbols);
H_FBMC_PostCoding(:,:,2,2) = reshape(MeanChannelMatrix*reshape(H_FBMC(:,:,2,2),[],1),L/2*NrBlocks,FBMC.Nr.MCSymbols);

%% Different noise levels
for i_SNR = 1:length(M_SNR_dB_OFDM)
SNR_dB_OFDM = M_SNR_dB_OFDM(i_SNR); 

Pn = OFDM.Implementation.FFTSize/OFDM.Nr.Subcarriers*10^(-SNR_dB_OFDM/10); 
noise_Antenna1 = sqrt(Pn/2)*(randn(N,1)+1j*randn(N,1));
noise_Antenna2 = sqrt(Pn/2)*(randn(N,1)+1j*randn(N,1));

%% Received signal
r_OFDM_Alamouti_Antenna1 = r_OFDM_Alamouti_Antenna1_NoNoise+noise_Antenna1;
r_FBMC_Alamouti_Antenna1 = r_FBMC_Alamouti_Antenna1_NoNoise+noise_Antenna2;

r_OFDM_SM_Antenna1 = r_OFDM_SM_Antenna1_NoNoise+noise_Antenna1;
r_OFDM_SM_Antenna2 = r_OFDM_SM_Antenna2_NoNoise+noise_Antenna2;

r_FBMC_SM_Antenna1 = r_FBMC_SM_Antenna1_NoNoise+noise_Antenna1;
r_FBMC_SM_Antenna2 = r_FBMC_SM_Antenna2_NoNoise+noise_Antenna2;

%% Demodulation
y_OFDM_Alamouti = OFDM.Demodulation(r_OFDM_Alamouti_Antenna1);
y_FBMC_Alamouti = FBMC.Demodulation(r_FBMC_Alamouti_Antenna1);
y_FBMC_Alamouti_PostCoding = reshape(Phi_Totak'*y_FBMC_Alamouti(:),L/2*NrBlocks,FBMC.Nr.MCSymbols);

y_OFDM_SM_Antenna1 = OFDM.Demodulation(r_OFDM_SM_Antenna1);
y_OFDM_SM_Antenna2 = OFDM.Demodulation(r_OFDM_SM_Antenna2);

y_FBMC_SM_Antenna1 = FBMC.Demodulation(r_FBMC_SM_Antenna1);
y_FBMC_SM_Antenna2 = FBMC.Demodulation(r_FBMC_SM_Antenna2);
y_FBMC_SM_Antenna1_PostCoding = reshape(Phi_Totak'*y_FBMC_SM_Antenna1(:),L/2*NrBlocks,FBMC.Nr.MCSymbols);
y_FBMC_SM_Antenna2_PostCoding = reshape(Phi_Totak'*y_FBMC_SM_Antenna2(:),L/2*NrBlocks,FBMC.Nr.MCSymbols);

%% Data Detection
x_est_OFDM_Alamouti = Alamouti.Decoder(y_OFDM_Alamouti,H_OFDM(:,:,1,1:2)*sqrt(2));
x_est_FBMC_Alamouti = Alamouti.Decoder(y_FBMC_Alamouti_PostCoding,H_FBMC_PostCoding(:,:,1,1:2)*sqrt(2));

% OFDM 
x_est_OFDM_SM_ZeroForcing = nan(OFDM.Nr.Subcarriers,OFDM.Nr.MCSymbols,2);
x_est_OFDM_SM_ML = nan(OFDM.Nr.Subcarriers,OFDM.Nr.MCSymbols,2);
for l_OFDM = 1:OFDM.Nr.Subcarriers
   for k_OFDM = 1:OFDM.Nr.MCSymbols
      y_temp = [y_OFDM_SM_Antenna1(l_OFDM,k_OFDM);y_OFDM_SM_Antenna2(l_OFDM,k_OFDM)];
      H_temp = permute(H_OFDM(l_OFDM,k_OFDM,1:2,1:2),[3 4 1 2]);
      
      % ZF Detection
      x_est_OFDM_SM_ZeroForcing(l_OFDM,k_OFDM,1:2) =  H_temp^-1*y_temp;
      
      % ML Detection
      [~,indexMin] = min(sum(abs(repmat(y_temp,1,QAM.ModulationOrder^2)-H_temp*ML_Mapping).^2,1),[],2);
      x_est_OFDM_SM_ML(l_OFDM,k_OFDM,1:2) = ML_Mapping(:,indexMin);
   end
end 

% FBMC
x_est_FBMC_SM_ZeroForcing = nan(L/2*NrBlocks,FBMC.Nr.MCSymbols,2);
x_est_FBMC_SM_ML = nan(L/2*NrBlocks,FBMC.Nr.MCSymbols,2);
for l_FBMC = 1:L/2*NrBlocks
   for k_FBMC = 1:FBMC.Nr.MCSymbols
      y_temp = [y_FBMC_SM_Antenna1_PostCoding(l_FBMC,k_FBMC);y_FBMC_SM_Antenna2_PostCoding(l_FBMC,k_FBMC)];
      H_temp = permute(H_FBMC_PostCoding(l_FBMC,k_FBMC,:,:),[3 4 1 2]);
      
      % ZF Detection
      x_est_FBMC_SM_ZeroForcing(l_FBMC,k_FBMC,1:2) =  H_temp^-1*y_temp;
      
      % ML Detection
      [~,indexMin] = min(sum(abs(repmat(y_temp,1,QAM.ModulationOrder^2)-H_temp*ML_Mapping).^2,1),[],2);
      x_est_FBMC_SM_ML(l_FBMC,k_FBMC,1:2) = ML_Mapping(:,indexMin);
   end
end 

% symbol-to-bit 
DetectedBitStream_OFDM_Alamouti = QAM.Symbol2Bit(x_est_OFDM_Alamouti);
DetectedBitStream_FBMC_Alamouti = QAM.Symbol2Bit(x_est_FBMC_Alamouti);

DetectedBitStream_OFDM_SM_ZeroForcing = QAM.Symbol2Bit(x_est_OFDM_SM_ZeroForcing);
DetectedBitStream_FBMC_SM_ZeroForcing = QAM.Symbol2Bit(x_est_FBMC_SM_ZeroForcing);

DetectedBitStream_OFDM_SM_ML = QAM.Symbol2Bit(x_est_OFDM_SM_ML);
DetectedBitStream_FBMC_SM_ML = QAM.Symbol2Bit(x_est_FBMC_SM_ML);

% BER
BER_OFDM_Alamouti(i_SNR,i_rep) = mean(BinaryDataStream_Alamouti_OFDM~=DetectedBitStream_OFDM_Alamouti);
BER_FBMC_Alamouti(i_SNR,i_rep) = mean(BinaryDataStream_Alamouti_FBMC~=DetectedBitStream_FBMC_Alamouti);

BER_OFDM_SM_ZeroForcing(i_SNR,i_rep) = mean(BinaryDataStream_SM_OFDM~=DetectedBitStream_OFDM_SM_ZeroForcing);
BER_FBMC_SM_ZeroForcing(i_SNR,i_rep) = mean(BinaryDataStream_SM_FBMC~=DetectedBitStream_FBMC_SM_ZeroForcing);

BER_OFDM_SM_ML(i_SNR,i_rep) = mean(BinaryDataStream_SM_OFDM~=DetectedBitStream_OFDM_SM_ML);
BER_FBMC_SM_ML(i_SNR,i_rep) = mean(BinaryDataStream_SM_FBMC~=DetectedBitStream_FBMC_SM_ML);

end

TimePassed = toc;
if mod(i_rep,1)==0
disp(['Realization ' int2str(i_rep) ' of ' int2str(NrRepetitions) '. Time left: ' int2str(TimePassed/i_rep*(NrRepetitions-i_rep)/60) 'minutes']);
end
end

%% Plot
LineWidth  = 1;
markersize = 8;
figure(1);
semilogy(M_SNR_dB_OFDM,mean(BER_FBMC_SM_ZeroForcing,2),'-x blue' ,'Markersize',markersize,'LineWidth',LineWidth);
hold on;box on;grid on;
semilogy(M_SNR_dB_OFDM,mean(BER_OFDM_SM_ZeroForcing,2),'-o blue' ,'Markersize',markersize,'LineWidth',LineWidth);
semilogy(M_SNR_dB_OFDM,mean(BER_FBMC_SM_ML,2)         ,'-+ black','Markersize',markersize,'LineWidth',LineWidth);
semilogy(M_SNR_dB_OFDM,mean(BER_OFDM_SM_ML,2)         ,'-d black','Markersize',markersize,'LineWidth',LineWidth);
semilogy(M_SNR_dB_OFDM,mean(BER_FBMC_Alamouti,2)      ,'--* red' ,'Markersize',markersize,'LineWidth',LineWidth);
semilogy(M_SNR_dB_OFDM,mean(BER_OFDM_Alamouti,2)      ,'--s red' ,'Markersize',markersize,'LineWidth',LineWidth);
xlabel('SNR(dB)');
ylabel('BER');
legend('Pruned DCT-p-FBMC (ZF)','CP-OFDM (ZF)','Pruned DCT-p-FBMC (ML)','CP-OFDM (ML)','Pruned DCT-p-FBMC (Alamouti Code)','CP-OFDM (Alamouti Code)','Location','SouthWest');
set(gca,'FontName','Times New Roman','FontSize',12);

%% It's not easy to simulate once. Save the data.
if Savdecision
    BER_OFDM_ZF           = mean(BER_OFDM_SM_ZeroForcing,2);
    BER_FBMC_ZF           = mean(BER_FBMC_SM_ZeroForcing,2);
    BER_OFDM_ML           = mean(BER_OFDM_SM_ML,2);
    BER_FBMC_ML           = mean(BER_FBMC_SM_ML,2);
    BER_OFDM_Alamouti     = mean(BER_OFDM_Alamouti,2);
    BER_FBMC_Alamouti     = mean(BER_FBMC_Alamouti,2);
    save('.\Results\BER_flat\BER_OFDM_ZF.mat',      'BER_OFDM_ZF'); 
    save('.\Results\BER_flat\BER_FBMC_ZF.mat',      'BER_FBMC_ZF'); 
    save('.\Results\BER_flat\BER_OFDM_ML.mat',      'BER_OFDM_ML'); 
    save('.\Results\BER_flat\BER_FBMC_ML.mat',      'BER_FBMC_ML'); 
    save('.\Results\BER_flat\BER_OFDM_Alamouti.mat','BER_OFDM_Alamouti'); 
    save('.\Results\BER_flat\BER_FBMC_Alamouti.mat','BER_FBMC_Alamouti'); 
end
