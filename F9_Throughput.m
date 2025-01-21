% =========================================================================   
% Author: Ying Wang
% Function: to simulate the throughput of multicarrier systems in SISO.
% e-mail address: wangyingstu@163.com
% =========================================================================   
clear; clc;close all;
addpath("Channel\","Coding\","MultiCarrierModulation\");
%% Parameters
M_SNR_dB            = -5:1.25:30;                   % SNR (dB)
NrRepetitions       = 1000;                         % Monte Carlo Repetitions                           
L                   = 64;                           % Number of subcarriers
SubcarrierSpacing   = 15e3;                         % Subcarrier spacing (15kHz)
CarrierFrequency    = 2.5e9;                        % Carrier Frequency (Hz)
K_FBMC              = 30;                           % FBMC No. of symbols
K_OFDMnoCP          = 15;                           % OFDM No. of symbols (no CP)
K_OFDM              = 14;                           % OFDM No. of symbols 
NrBlocks            = 4;                            % Frame Number of Transmission
NrGuardSubcarriers  = 0;                            % Number of guard subcarriers
CP_Length           = 1/SubcarrierSpacing/14;       % LTE CP Length (seconds)
SamplingRate        = 15e3*14*5*NrBlocks;           % Sampling rate
PowerDelayProfile   = 'TDL-A_300ns';                % Power Delay Profile 
Velocity_kmh        = 200;                          % Velocity (km/h)
Savdecision         = 0;                            % Whether to save data

%% Adaptive Modulation and Coding (CQI List)
% The first column denotes the modulation order. 4, 16, 64, 256, 1024...
% The second column denotes the code rate (must be between 0 and 1)
M_CQI = [4  ,  78/1024;...
         4  , 120/1024;...
         4  , 193/1024;...
         4  , 308/1024;...
         4  , 449/1024;...
         4  , 602/1024;...
         16 , 378/1024;...
         16 , 490/1024;...
         16 , 616/1024;...
         64 , 466/1024;...
         64 , 567/1024;...
         64 , 666/1024;...
         64 , 772/1024;...
         64 , 873/1024;...
         64 , 948/1024]; % http://www.etsi.org/deliver/etsi_ts/136200_136299/136213/08.08.00_60/ts_136213v080800p.

if not(strcmp(mexext,'mexw64'))  
    IndexCodeRateSmallerOneThird =  find(M_CQI(:,2)<1/3);
    if  numel(IndexCodeRateSmallerOneThird)   >0
        M_CQI(IndexCodeRateSmallerOneThird,:) = [];
        warning('A code rate smaller than 1/3 is only supported in Windows 64-bit => CQI values which contain a code rate smaller than 1/3 are discarded!');
    end    
end
    
%% FBMC
FBMC_DCT = FBMC(...
    (L+NrGuardSubcarriers)*NrBlocks,...             % Number of subcarriers
    K_FBMC,...                                      % FBMC No. of symbols
    SubcarrierSpacing,...                           % Subcarrier spacing (Hz)
    SamplingRate,...                                % Sampling rate (Samples/s)
    0,...                                           % Intermediate frequency of the first subcarrier (Hz).
    false,...                                       % Transmit real valued signal
    'PHYDYAS-OQAM',...                              % Prototype filter and OQAM or QAM.
    4, ...                                          % Overlapping factor
    0, ...                                          % Initial phase shift
    true ...                                        % Efficient IFFT implementation, true or false
    );

ConFBMC = FBMC(...
    (L+NrGuardSubcarriers)*NrBlocks,...             % Number of subcarriers
    K_FBMC,...                                      % FBMC No. of symbols
    SubcarrierSpacing,...                           % Subcarrier spacing (Hz)
    SamplingRate,...                                % Sampling rate (Samples/s)
    0,...                                           % Intermediate frequency of the first subcarrier (Hz).
    false,...                                       % Transmit real valued signal
    'PHYDYAS-OQAM',...                              % Prototype filter and OQAM or QAM. 
    4, ...                                          % Overlapping factor
    0, ...                                          % Initial phase shift
    true ...                                        % Efficient IFFT implementation, true or false
    );


%% OFDM 
ZeroGuardTimeLength = ((ConFBMC.Nr.SamplesTotal-(round(SamplingRate/SubcarrierSpacing)+0*SamplingRate)*K_OFDMnoCP)/2)/SamplingRate;
OFDMnoCP = OFDM(...
    (L+NrGuardSubcarriers)*NrBlocks,...             % Number of subcarriers
    K_OFDMnoCP,...                                  % OFDM No. of symbols
    SubcarrierSpacing,...                           % Subcarrier spacing (Hz)
    SamplingRate,...                                % Sampling rate (Samples/s)
    0,...                                           % Intermediate frequency of the first subcarrier (Hz).
    false,...                                       % Transmit real valued signal
    0, ...                                          % Cycle prefix length (s) 1/SubcarrierSpacing/(K/2-1)
    ZeroGuardTimeLength ...                         % Length of the guard time (s)
    );
ZeroGuardTimeLength = ((ConFBMC.Nr.SamplesTotal-(round(SamplingRate/SubcarrierSpacing)+CP_Length*SamplingRate)*K_OFDM)/2)/SamplingRate;
CPOFDM = OFDM(...
   (L+NrGuardSubcarriers)*NrBlocks,...              % Number of subcarriers
    K_OFDM,...                                      % OFDM No. of symbols
    SubcarrierSpacing,...                           % Subcarrier spacing (Hz)
    SamplingRate,...                                % Sampling rate (Samples/s)
    0,...                                           % Intermediate frequency of the first subcarrier (Hz).
    false,...                                       % Transmit real valued signal
    CP_Length, ...                                  % Cycle prefix length (s) 1/SubcarrierSpacing/(K/2-1)
    ZeroGuardTimeLength ...                         % Length of the guard time (s)
    );

%% Check the number of samples
if  CPOFDM.Nr.SamplesTotal~=ConFBMC.Nr.SamplesTotal || OFDMnoCP.Nr.SamplesTotal~=ConFBMC.Nr.SamplesTotal
   error('Total number of samples must be the same for OFDM and FBMC.');
end
N = CPOFDM.Nr.SamplesTotal;

%% Channel
ChannelModel = FastFading(...
    SamplingRate,...                                                    % Sampling rate (Samples/s)
    PowerDelayProfile,...                                               % Power delay profile 
    N,...                                                               % Number of total samples
    Velocity_kmh/3.6*CarrierFrequency/2.998e8,...                       % Maximum Doppler shift: Velocity_kmh/3.6*CarrierFrequency/2.998e8  
    'Jakes',...                                                         % Which Doppler model: 'Jakes', 'Uniform', 'Discrete-Jakes', 'Discrete-Uniform'.                                      
    200, ...                                                            % Number of paths for the WSSUS process.                          
    1,...                                                               % Number of transmit antennas
    1,...                                                               % Number of receive antennas
    0 ...                                                               % Gives a warning if the predefined delay taps of the channel do not fit the sampling rate.
    );

%% Precoding matrix
NrSubcarriers = ConFBMC.Nr.Subcarriers; % Total number of subcarriers==>NrSubcarriers = (L+NrGuardSubcarriers)*NrBlocks£©
NrSymbols     = ConFBMC.Nr.MCSymbols;
Domain        = 'Time';       % Setup Domain
Mathod        = 'DCT';        % Adopting the DCT method
Kall          = K_FBMC;       % No. of symbols
K             = Kall;         % FBMC No. of symbols

phi = GetCodingMatrix(...
    Domain,...          % Determine the domain ('Frequency', 'Time', 'DeEigen')
    K,...               % FBMC No. of symbols
    L,...               % Number of subcarriers
    Mathod,...          % Precoding matrix form
    ConFBMC...          % The functions are embedded in the FBMC system.
    );

ConFBMC.SetNrSubcarriers(NrSubcarriers);
ConFBMC.SetNrMCSymbols(NrSymbols);
% Mapping precoding Matrix to large Matrix describing the whole system.
Phi = kron(sparse(eye(NrBlocks)),phi); 

%% Generate DFT Matrix
DFTMatrix               = fft(eye(NrSubcarriers))/sqrt(NrSubcarriers);

%% Pre-Initialized CQI: Turbo Encoder with QAM
% Space pre-configuration
QAM                  = cell(1,size(M_CQI,1));
PAM                  = cell(1,size(M_CQI,1));
TurboCoding_OFDM     = cell(1,size(M_CQI,1));
TurboCoding_OFDMnoCP = cell(1,size(M_CQI,1));
TurboCoding_FBMC     = cell(1,size(M_CQI,1));
TurboCoding_FBMC_DFT = cell(1,size(M_CQI,1));
TurboCoding_DCT_FBMC  = cell(1,size(M_CQI,1));

for i_cqi = 1:size(M_CQI,1)
    QAMModulationOrder          = M_CQI(i_cqi,1);
    PAMModulationOrder          = sqrt(QAMModulationOrder);
    CodeRate                    = M_CQI(i_cqi,2);
    
    QAM{i_cqi}                  = SignalConstellation(QAMModulationOrder,'QAM');
    PAM{i_cqi}                  = SignalConstellation(PAMModulationOrder,'PAM');

    NrTransmittedBits_OFDM      = NrSubcarriers*K_OFDM*log2(QAMModulationOrder);
    NrTransmittedBits_OFDMnoCP  = NrSubcarriers*K_OFDMnoCP*log2(QAMModulationOrder);
    NrTransmittedBits_FBMC      = NrSubcarriers*K_FBMC*log2(PAMModulationOrder);

    NrTransmittedBits_DCT_FBMC   = (NrSubcarriers-NrGuardSubcarriers*NrBlocks)/2*K_FBMC*log2(QAMModulationOrder);
  
    TurboCoding_OFDM{i_cqi}     = TurboCoding( NrTransmittedBits_OFDM     , round(CodeRate*NrTransmittedBits_OFDM));
    TurboCoding_OFDMnoCP{i_cqi} = TurboCoding( NrTransmittedBits_OFDMnoCP , round(CodeRate*NrTransmittedBits_OFDMnoCP));
    TurboCoding_FBMC{i_cqi}     = TurboCoding( NrTransmittedBits_FBMC     , round(CodeRate*NrTransmittedBits_FBMC));
    TurboCoding_DCT_FBMC{i_cqi} = TurboCoding( NrTransmittedBits_DCT_FBMC  , round(CodeRate*NrTransmittedBits_DCT_FBMC));

end

%% Obtain the transmitted and received matrices of OFDM and FBMC.
GTX_OFDM     = sparse(CPOFDM.GetTXMatrix);
GRX_OFDM     = sparse(CPOFDM.GetRXMatrix');
GTX_OFDMnoCP = sparse(OFDMnoCP.GetTXMatrix);
GRX_OFDMnoCP = sparse(OFDMnoCP.GetRXMatrix');
G_FBMC       = sparse(ConFBMC.GetTXMatrix);
G_DCT_FBMC   = sparse(FBMC_DCT.GetTXMatrix);

%% Normalized OFDM and FBMC
NormalizationOFDM     = sqrt((GRX_OFDM(:,1)'*GRX_OFDM(:,1)));
NormalizationOFDMnoCP = sqrt((GRX_OFDMnoCP(:,1)'*GRX_OFDMnoCP(:,1)));
NormalizationFBMC     = 1/sqrt((G_FBMC(:,1)'*G_FBMC(:,1)));
NormalizationDCT_FBMC  = 1/sqrt((G_DCT_FBMC(:,1)'*G_DCT_FBMC(:,1)));

GTX_OFDM     = GTX_OFDM*NormalizationOFDM;
GRX_OFDM     = GRX_OFDM/NormalizationOFDM;
GTX_OFDMnoCP = GTX_OFDMnoCP*NormalizationOFDMnoCP;
GRX_OFDMnoCP = GRX_OFDMnoCP/NormalizationOFDMnoCP;
G_FBMC       = G_FBMC*NormalizationFBMC;
G_DCT_FBMC  = G_DCT_FBMC*NormalizationDCT_FBMC;

%% Space pre-configuration
M_Througput_OFDM            = nan( length(M_SNR_dB) , NrRepetitions , size(M_CQI,1) );
M_Througput_OFDMnoCP        = nan( length(M_SNR_dB) , NrRepetitions , size(M_CQI,1) );
M_Througput_FBMC            = nan( length(M_SNR_dB) , NrRepetitions , size(M_CQI,1) );
M_Througput_DFT_OFDM        = nan( length(M_SNR_dB) , NrRepetitions , size(M_CQI,1) );
M_Througput_DCT_FBMC       = nan( length(M_SNR_dB) , NrRepetitions , size(M_CQI,1) );

disp('The simulation may take a while ... ');
%% Start Simulation  
tic;
for i_Rep = 1:NrRepetitions                           
    % Channel 
    ChannelModel.NewRealization;
    H = ChannelModel.GetConvolutionMatrix{1};
    noise_unitPower = sqrt(1/2)*( randn(N,1) + 1j * randn(N,1) );
      
    % Calculate one-tap channels
    % h == sum((G'*H).*G.',2)
    h_OFDM      = full( reshape( sum((GRX_OFDM'*H).*GTX_OFDM.',2), NrSubcarriers, [] ) ); 
    h_OFDMnoCP  = full( reshape( sum((GRX_OFDMnoCP'*H).*GTX_OFDMnoCP.',2), NrSubcarriers, [] ) ); 
    h_FBMC      = full( reshape( sum((G_FBMC'*H).*G_FBMC.',2), NrSubcarriers, [] ) ); 
    h_DCT_FBMC = full( reshape( sum((G_DCT_FBMC'*H).*G_DCT_FBMC.',2), NrSubcarriers, [] ) );
    
    M_Througput_OFDM_OneRealization       = nan( length(M_SNR_dB) , size(M_CQI,1) );
    M_Througput_OFDMnoCP_OneRealization   = nan( length(M_SNR_dB) , size(M_CQI,1) );
    M_Througput_FBMC_OneRealization       = nan( length(M_SNR_dB) , size(M_CQI,1) );
    M_Througput_DFT_OFDM_OneRealization     = nan( length(M_SNR_dB) , size(M_CQI,1) );
    M_Througput_DCT_FBMC_OneRealization  = nan( length(M_SNR_dB) , size(M_CQI,1) );  
   
    % Simulation for different modulation orders and code rates
    for i_cqi = 1:size(M_CQI,1)
        % Generate Bitstream
        BinaryDataStream_OFDM     = randi( [0 1] , TurboCoding_OFDM{i_cqi}.NrDataBits , 1 );
        BinaryDataStream_OFDMnoCP = randi( [0 1] , TurboCoding_OFDMnoCP{i_cqi}.NrDataBits , 1 );
        BinaryDataStream_FBMC     = randi( [0 1] , TurboCoding_FBMC{i_cqi}.NrDataBits , 1 );
        BinaryDataStream_DCT_FBMC  = randi( [0 1] , TurboCoding_DCT_FBMC{i_cqi}.NrDataBits , 1 );
        
        % Updating the intertwining of Turbo encoders
        TurboCoding_OFDM{i_cqi}.UpdateInterleaving;
        TurboCoding_OFDMnoCP{i_cqi}.UpdateInterleaving;        
        TurboCoding_FBMC{i_cqi}.UpdateInterleaving;
        TurboCoding_DCT_FBMC{i_cqi}.UpdateInterleaving;
        
        % Turbo Encoding of Data Bits
        CodedBits_OFDM      = TurboCoding_OFDM{i_cqi}.TurboEncoder( BinaryDataStream_OFDM );
        CodedBits_OFDMnoCP  = TurboCoding_OFDMnoCP{i_cqi}.TurboEncoder( BinaryDataStream_OFDMnoCP );
        CodedBits_FBMC      = TurboCoding_FBMC{i_cqi}.TurboEncoder( BinaryDataStream_FBMC );
        CodedBits_DCT_FBMC   = TurboCoding_DCT_FBMC{i_cqi}.TurboEncoder( BinaryDataStream_DCT_FBMC );
        
        % bit intertwined
        BitInterleaving_OFDM     = randperm( TurboCoding_OFDM{i_cqi}.NrCodedBits );
        BitInterleaving_OFDMnoCP = randperm( TurboCoding_OFDMnoCP{i_cqi}.NrCodedBits );
        BitInterleaving_FBMC     = randperm( TurboCoding_FBMC{i_cqi}.NrCodedBits );
        BitInterleaving_DCT_FBMC  = randperm( TurboCoding_DCT_FBMC{i_cqi}.NrCodedBits );
        
        CodedBits_OFDM     = CodedBits_OFDM(     BitInterleaving_OFDM );
        CodedBits_OFDMnoCP = CodedBits_OFDMnoCP( BitInterleaving_OFDMnoCP );
        CodedBits_FBMC     = CodedBits_FBMC(     BitInterleaving_FBMC );
        CodedBits_DCT_FBMC  = CodedBits_DCT_FBMC(  BitInterleaving_DCT_FBMC );
             
        % Mapping bitstreams to data symbols
        x_OFDM     = reshape(QAM{i_cqi}.Bit2Symbol(CodedBits_OFDM)    , NrSubcarriers,K_OFDM);
        x_OFDMnoCP = reshape(QAM{i_cqi}.Bit2Symbol(CodedBits_OFDMnoCP), NrSubcarriers,K_OFDMnoCP);
        x_FBMC     = reshape(PAM{i_cqi}.Bit2Symbol(CodedBits_FBMC)    , NrSubcarriers,K_FBMC)/sqrt(2);  % 1/sqrt(2) => same TX power for OFDM and FBMC
        x_DCT_FBMC  = reshape(QAM{i_cqi}.Bit2Symbol(CodedBits_DCT_FBMC) ,(NrSubcarriers-NrGuardSubcarriers*NrBlocks)/2,K_FBMC);
        
        % Generation of time-domain transmitted signals
        s_FBMC         = ConFBMC.Modulation(x_FBMC)*NormalizationFBMC;
        s_DCT_FBMC      = ConFBMC.Modulation(Phi*x_DCT_FBMC)*NormalizationFBMC;
        s_DFT_OFDM     = CPOFDM.Modulation(DFTMatrix*x_OFDM)*NormalizationOFDM;
        s_OFDM         = CPOFDM.Modulation(x_OFDM)*NormalizationOFDM;
        s_OFDMnoCP     = OFDMnoCP.Modulation(x_OFDMnoCP)*NormalizationOFDMnoCP;

        % pass through channel    
        r_FBMC_noNoise         = H * s_FBMC;
        r_DCT_FBMC_noNoise    = H * s_DCT_FBMC;
        r_DFT_OFDM_noNoise     = H * s_DFT_OFDM;
        r_OFDM_noNoise         = H * s_OFDM;
        r_OFDMnoCP_noNoise     = H * s_OFDMnoCP;

        % Simulation with different SNR
        for i_SNR = 1:length(M_SNR_dB)
            SNR_dB  = M_SNR_dB(i_SNR); 
            Pn      = 10^(-SNR_dB/10); 
            % add noise
            noise          = sqrt(Pn) * noise_unitPower;
            r_FBMC         = r_FBMC_noNoise     + noise;
            r_DCT_FBMC    = r_DCT_FBMC_noNoise  + noise;
            r_DFT_OFDM     = r_DFT_OFDM_noNoise     + noise;
            r_OFDM         = r_OFDM_noNoise     + noise;
            r_OFDMnoCP     = r_OFDMnoCP_noNoise + noise;

            % Received Symbol (Demodulation)
            y_FBMC           = ConFBMC.Demodulation(r_FBMC)/NormalizationFBMC;
            y_DCT_FBMC      = FBMC_DCT.Demodulation(r_DCT_FBMC)/NormalizationDCT_FBMC;
            y_DFT_OFDM       = CPOFDM.Demodulation(r_DFT_OFDM)/NormalizationOFDM;
            y_OFDM           = CPOFDM.Demodulation(r_OFDM)/NormalizationOFDM;
            y_OFDMnoCP       = OFDMnoCP.Demodulation(r_OFDMnoCP)/NormalizationOFDMnoCP;

            % FBMC's Zero Force Equalizer
            x_est_FBMC       = real(y_FBMC./h_FBMC)*sqrt(2);
            % MMSE Equalizer 
            Scaling_OFDM     =  repmat(1./(mean(1 ./( 1 + Pn./abs( h_OFDM ).^2 ),1)),     NrSubcarriers,1);
            Scaling_OFDMnoCP =  repmat(1./(mean(1 ./( 1 + Pn./abs( h_OFDMnoCP ).^2 ),1)), NrSubcarriers,1);
            Scaling_DFT_OFDM =  repmat(1./(mean(1 ./( 1 + Pn./abs( h_OFDM ).^2 ),1)), NrSubcarriers,1);
            Scaling_DCT_FBMC  =  repmat(1./(mean(1 ./( 1 + Pn./abs( h_DCT_FBMC ).^2 ),1)),  NrSubcarriers,1);
            
            e_OFDM           = Scaling_OFDM     .*conj(h_OFDM)    ./( abs(h_OFDM).^2     + Pn );
            e_OFDMnoCP       = Scaling_OFDMnoCP .*conj(h_OFDMnoCP)./( abs(h_OFDMnoCP).^2 + Pn );
            e_DFT_OFDM       = Scaling_DFT_OFDM     .*conj(h_OFDM)    ./( abs(h_OFDM).^2     + Pn );
            e_DCT_FBMC      = Scaling_DCT_FBMC  .*conj(h_DCT_FBMC )./( abs(h_DCT_FBMC ).^2 + Pn );

            x_est_OFDM       = y_OFDM     .* e_OFDM;
            x_est_OFDMnoCP   = y_OFDMnoCP .* e_OFDMnoCP;
            x_est_DFT_OFDM   = DFTMatrix'  * (y_DFT_OFDM   .* e_DFT_OFDM);
            x_est_DCT_FBMC  = Phi'   * (y_DCT_FBMC  .* e_DCT_FBMC);

            % Calculate the LLR value
            AWGNequivalentNoise_OFDM         = repmat(1./mean(1./(1+Pn./abs(h_OFDM).^2),1)-1 ,NrSubcarriers,1);
            AWGNequivalentNoise_OFDMnoCP     = repmat(1./mean(1./(1+Pn./abs(h_OFDMnoCP).^2),1)-1,NrSubcarriers,1);
            AWGNequivalentNoise_DCT_FBMC      = repmat(1./mean(1./(1+Pn./abs(h_DCT_FBMC).^2),1)-1 ,NrSubcarriers/2,1);
            AWGNequivalentNoise_DFT_OFDM     = repmat(1./mean(1./(1+Pn./abs(h_OFDM).^2),1)-1 ,NrSubcarriers,1);
                       
            LLR_FBMC         = PAM{i_cqi}.LLR_AWGN( x_est_FBMC(:)     , 2*Pn .* 1./abs(h_FBMC(:)).^2);
            LLR_OFDM         = QAM{i_cqi}.LLR_AWGN( x_est_OFDM(:)     , AWGNequivalentNoise_OFDM(:));
            LLR_OFDMnoCP     = QAM{i_cqi}.LLR_AWGN( x_est_OFDMnoCP(:) , AWGNequivalentNoise_OFDMnoCP(:));
            LLR_DFT_OFDM     = QAM{i_cqi}.LLR_AWGN( x_est_DFT_OFDM(:)     , AWGNequivalentNoise_DFT_OFDM(:));
            LLR_DCT_FBMC      = QAM{i_cqi}.LLR_AWGN( x_est_DCT_FBMC(:)  , AWGNequivalentNoise_DCT_FBMC(:));
                       
            % bitstream deinterleaving
            LLR_FBMC(BitInterleaving_FBMC)         = LLR_FBMC;
            LLR_DCT_FBMC(BitInterleaving_DCT_FBMC)   = LLR_DCT_FBMC;
            LLR_DFT_OFDM(BitInterleaving_OFDM)         = LLR_DFT_OFDM;
            LLR_OFDM(BitInterleaving_OFDM)         = LLR_OFDM;
            LLR_OFDMnoCP(BitInterleaving_OFDMnoCP) = LLR_OFDMnoCP;          

            % Decode Bitstreams            
            DecodedBits_FBMC     = TurboCoding_FBMC{i_cqi}.TurboDecoder( LLR_FBMC );
            DecodedBits_DCT_FBMC  = TurboCoding_DCT_FBMC{i_cqi}.TurboDecoder( LLR_DCT_FBMC );         
            DecodedBits_OFDM     = TurboCoding_OFDM{i_cqi}.TurboDecoder( LLR_OFDM );
            DecodedBits_OFDMnoCP = TurboCoding_OFDMnoCP{i_cqi}.TurboDecoder( LLR_OFDMnoCP );
            DecodedBits_DFT_OFDM     = TurboCoding_OFDM{i_cqi}.TurboDecoder( LLR_DFT_OFDM );          
            
            % Calculate throughput after decoding 
            M_Througput_FBMC_OneRealization(i_SNR,i_cqi)     = all( DecodedBits_FBMC     == BinaryDataStream_FBMC ) * length(BinaryDataStream_FBMC)/(ConFBMC.PHY.TimeSpacing*(ConFBMC.Nr.MCSymbols));
            M_Througput_DCT_FBMC_OneRealization(i_SNR,i_cqi)  = all( DecodedBits_DCT_FBMC  == BinaryDataStream_DCT_FBMC ) * length(BinaryDataStream_DCT_FBMC)/(ConFBMC.PHY.TimeSpacing*(ConFBMC.Nr.MCSymbols));
 
            M_Througput_OFDM_OneRealization(i_SNR,i_cqi)     = all( DecodedBits_OFDM     == BinaryDataStream_OFDM ) * length(BinaryDataStream_OFDM)/(CPOFDM.PHY.TimeSpacing*(CPOFDM.Nr.MCSymbols));
            M_Througput_OFDMnoCP_OneRealization(i_SNR,i_cqi) = all( DecodedBits_OFDMnoCP == BinaryDataStream_OFDMnoCP ) * length(BinaryDataStream_OFDMnoCP)/(OFDMnoCP.PHY.TimeSpacing*(OFDMnoCP.Nr.MCSymbols));
            M_Througput_DFT_OFDM_OneRealization(i_SNR,i_cqi)     = all( DecodedBits_DFT_OFDM  == BinaryDataStream_OFDM ) * length(BinaryDataStream_OFDM)/(CPOFDM.PHY.TimeSpacing*(CPOFDM.Nr.MCSymbols));
            
        end
    end
    M_Througput_FBMC(:,i_Rep,:)      = M_Througput_FBMC_OneRealization;
    M_Througput_DCT_FBMC(:,i_Rep,:)   = M_Througput_DCT_FBMC_OneRealization;
    
    M_Througput_OFDM(:,i_Rep,:)      = M_Througput_OFDM_OneRealization;
    M_Througput_OFDMnoCP(:,i_Rep,:)  = M_Througput_OFDMnoCP_OneRealization; 
    M_Througput_DFT_OFDM(:,i_Rep,:)      = M_Througput_DFT_OFDM_OneRealization;
   
    TimePassed = toc;
    disp(['Realization to ' int2str(i_Rep) ' of ' int2str(NrRepetitions) ' used ' int2str(TimePassed) 's. Time left:'...
          int2str(TimePassed*(NrRepetitions-i_Rep)/i_Rep/60) 'minutes, Approx: '...
          int2str(TimePassed*(NrRepetitions-i_Rep)/i_Rep/3600) 'hours']);
end

%% Maximization under CQI => Perfect Feedback
Througput_FBMC      = max(M_Througput_FBMC,[],3);
Througput_DCT_FBMC   = max(M_Througput_DCT_FBMC,[],3);
Througput_DFT_OFDM     = max(M_Througput_DFT_OFDM,[],3);

Througput_OFDM      = max(M_Througput_OFDM,[],3);
Througput_OFDMnoCP  = max(M_Througput_OFDMnoCP,[],3);
%% Plot throughput
LineWidth1 = 1;
LineWidth  = 0.75;
figure(1);
plot(M_SNR_dB,mean(Througput_OFDM,2)/1e6     ,'-'  ,'Color',0.85*[1,0,0],'LineWidth',LineWidth1); hold on;grid on;
plot(M_SNR_dB,mean(Througput_OFDMnoCP,2)/1e6 ,'-.' ,'Color',0.75*[1,0,1],'LineWidth',LineWidth1); 
plot(M_SNR_dB,mean(Througput_FBMC,2)/1e6     ,'-'  ,'Color',0.85*[0,1,0],'LineWidth',LineWidth1);
plot(M_SNR_dB,mean(Througput_DFT_OFDM,2)/1e6,':'  ,'Color',0.85*[1,0,0],'LineWidth',LineWidth1+0.3);
plot(M_SNR_dB,mean(Througput_DCT_FBMC,2)/1e6  ,'-'  ,'Color',0.85*[0,0,1],'LineWidth',LineWidth1); 
grid on;box on;
xlabel('SNR (dB)');
ylabel('Throughput (Mbit/s)');
title('SISO throughput, Doubly-selective channel');
s1 = plot([nan nan],[nan nan],'-' ,'Color',0.85*[1,0,0],'LineWidth',LineWidth);
s2 = plot([nan nan],[nan nan],'-.','Color',0.75*[1,0,1],'LineWidth',LineWidth);
s3 = plot([nan nan],[nan nan],'-' ,'Color',0.85*[0,1,0],'LineWidth',LineWidth);
s4 = plot([nan nan],[nan nan],':','Color',0.85*[1,0,0],'LineWidth',LineWidth+0.3);
s5 = plot([nan nan],[nan nan],'-' ,'Color',0.85*[0,0,1],'LineWidth',LineWidth);
legend([s1 s2 s3 s4 s5],{'OFDM (with CP)','OFDM (no CP)','Conventional FBMC','SC-FDMA','Pruned DCT-p-FBMC'},'Location','NorthWest');
set(gca,'FontName','Times New Roman','FontSize',12,'GridLineStyle','-.');

%% It's not easy to simulate once. Save the data.
if Savdecision
    meanThrougput_OFDM     = mean(Througput_OFDM,2);
    meanThrougput_OFDMnoCP = mean(Througput_OFDMnoCP,2);
    meanThrougput_FBMC     = mean(Througput_FBMC,2);
    meanThrougput_DFT_OFDM  = mean(Througput_DFT_OFDM,2);
    meanThrougput_DCT_FBMC  = mean(Througput_DCT_FBMC,2);
    save('.\Results\SISO_Throughput\meanThrougput_OFDM.mat',    'meanThrougput_OFDM'); 
    save('.\Results\SISO_Throughput\meanThrougput_OFDMnoCP.mat','meanThrougput_OFDMnoCP'); 
    save('.\Results\SISO_Throughput\meanThrougput_FBMC.mat',    'meanThrougput_FBMC'); 
    save('.\Results\SISO_Throughput\meanThrougput_DFT_OFDM.mat', 'meanThrougput_DFT_OFDM'); 
    save('.\Results\SISO_Throughput\meanThrougput_DCT_FBMC.mat','meanThrougput_DCT_FBMC'); 
end

    