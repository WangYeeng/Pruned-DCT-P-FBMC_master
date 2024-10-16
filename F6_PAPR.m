%==========================================================================
% Author: Ying Wang
% Function: Simulation for different multicarrier modulations with peak-to-average power ratios.
% e-mail address: wangyingstu@163.com
%==========================================================================
clc;clear; close all;
addpath("Channel\","Coding\","MultiCarrierModulation\","Results\");
load('CCDF_FBMC_DFT_LessPoints.mat');
load('CCDF_FBMC_DFT_xAxis_LessPoints.mat');
%% Parameters
NrRepetitions           = 100000;                         % Monte Carlo Repetitions.
NrSubcarriers           = 64;                           % Number of subcarriers
L                       = NrSubcarriers;
QAM_ModulationOrder     = 16;                           % Modulation order, 4, 16, 64,...
SubcarrierSpacing       = 15e3;                         % Subcarrier spacing (15kHz)
CarrierFrequency        = 2.5e9;                        % Carrier Frequency (Hz)
NrBlocks                = 4;                            % Frame Number of Transmission
NrGuardSubcarriers      = 0;                            % Number of guard subcarriers
K_FBMC                  = 30;                           % FBMC No. of symbols
K_OFDMnoCP              = 15;                           % OFDM No. of symbols (no CP)
K_OFDM                  = 14;                           % OFDM No. of symbols
CP_Length               = 1/SubcarrierSpacing/14;       % LTE CP Length (seconds)
SamplingRate            = SubcarrierSpacing*L*21;       % Sampling rate
PseudoOverlappingFactor = 4;                            % Pseudo Overlapping Factor. 

%% FBMC
ConFBMC = FBMC(...
    (L+NrGuardSubcarriers)*NrBlocks,...             % Number of subcarriers
    K_FBMC,...                                      % FBMC No. of symbols
    SubcarrierSpacing,...                           % Subcarrier spacing (Hz)
    SamplingRate,...                                % Sampling rate (Samples/s)
    0,...                                           % Intermediate frequency of the first subcarrier (Hz).
    false,...                                       % Transmit real valued signal
    'PHYDYAS-OQAM',...                              % Prototype filter and OQAM or QAM.
    PseudoOverlappingFactor, ...                    % Overlapping factor
    0, ...                                          % Initial phase shift
    true ...                                        % Efficient IFFT implementation, true or false
    );

FBMC_DCT = FBMC(...          
    (L+NrGuardSubcarriers)*NrBlocks,...             % Number of subcarriers
    K_FBMC,...                                      % FBMC No. of symbols
    SubcarrierSpacing,...                           % Subcarrier spacing (Hz)
    SamplingRate,...                                % Sampling rate (Samples/s)
    0,...                                           % Intermediate frequency of the first subcarrier (Hz).
    false,...                                       % Transmit real valued signal
    'PHYDYAS-OQAM',...                              % Prototype filter and OQAM or QAM.
    PseudoOverlappingFactor, ...                    % Overlapping factor
    0, ...                                          % Initial phase shift
    true ...                                        % Efficient IFFT implementation, true or false
    );


%% OFDM Object
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
    (L+NrGuardSubcarriers)*NrBlocks,...             % Number of subcarriers
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

%% Symbol constellations
QAM = SignalConstellation(QAM_ModulationOrder,'QAM');
PAM = SignalConstellation(sqrt(QAM_ModulationOrder),'PAM');

%% Precoding matrix
NrSubcarriers = ConFBMC.Nr.Subcarriers;
NrSymbols     = ConFBMC.Nr.MCSymbols;
Domain        = 'Time-DFT';     % Setup Domain
Mathod        = 'DCT';      % Adopting the DCT method
Kall          = K_FBMC;     % No. of symbols
K             = Kall;       % FBMC No. of symbols

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
Phi_OneFBMCsymbol = kron(sparse(eye(NrBlocks)),phi); 
Phi_Totak = kron(eye(K_FBMC),Phi_OneFBMCsymbol);
% Power normalization ensures that OFDM and FBMC are identical.
FBMCPowerNormalization = sqrt(2*(NrSubcarriers+NrGuardSubcarriers)/NrSubcarriers);
MeanChannelMatrix = abs(Phi_Totak').^2;

%% Generate DFT Matrix
DFTMatrix               = fft(eye(NrSubcarriers))/sqrt(NrSubcarriers);

%% Pre-configuration
AveragePower_OFDM         = nan(1,NrRepetitions);
AveragePower_OFDMnoCP     = nan(1,NrRepetitions);
AveragePower_FBMC         = nan(1,NrRepetitions);
AveragePower_DFT_OFDM     = nan(1,NrRepetitions);
AveragePower_DCT_FBMC    = nan(1,NrRepetitions);

PeakPower_OFDM          = nan(NrRepetitions,K_OFDM);
PeakPower_OFDMnoCP      = nan(NrRepetitions,K_OFDMnoCP);
PeakPower_FBMC          = nan(NrRepetitions,K_FBMC/2-2);
PeakPower_DFT_OFDM      = nan(NrRepetitions,K_OFDM);
PeakPower_DCT_FBMC     = nan(NrRepetitions,K_FBMC/2);

%% Start Simulation
tic;
for i_rep = 1:NrRepetitions   
    %% Generate Bitstream
    BinaryDataStream              = randi( [0 1] , K_OFDMnoCP*NrSubcarriers*log2(QAM.ModulationOrder) ,1);
    BinaryDataStream_OFDM         = randi( [0 1] , K_OFDM*NrSubcarriers*log2(QAM.ModulationOrder) ,1); 
    BinaryDataStream_DCT_FBMC    = randi( [0 1] , K_FBMC*L/2*NrBlocks*log2(QAM.ModulationOrder),1);

    %% Mapping bitstreams to transmitted symbols
    x_OFDM      = reshape( QAM.Bit2Symbol(BinaryDataStream_OFDM) , NrSubcarriers, K_OFDM );
    x_OFDMnoCP  = reshape( QAM.Bit2Symbol(BinaryDataStream) , NrSubcarriers , K_OFDMnoCP );
    x_FBMC      = reshape( PAM.Bit2Symbol(BinaryDataStream), NrSubcarriers , K_FBMC );
    x_FBMC_DCT_FBMC   = reshape( QAM.Bit2Symbol(BinaryDataStream_DCT_FBMC), L/2*NrBlocks, K_FBMC);

    %% Generation of time-domain transmitted signals
    s_OFDM              = CPOFDM.Modulation( x_OFDM );
    s_OFDMnoCP          = OFDMnoCP.Modulation( x_OFDMnoCP );
    s_FBMC              = ConFBMC.Modulation( x_FBMC ); 
    s_DFT_OFDM          = CPOFDM.Modulation( DFTMatrix*x_OFDM );

    TransmittedSymbols_DCT_FBMC = reshape(Phi_Totak*x_FBMC_DCT_FBMC(:),ConFBMC.Nr.Subcarriers,ConFBMC.Nr.MCSymbols);
    s_DCT_FBMC           = FBMCPowerNormalization*FBMC_DCT.Modulation( TransmittedSymbols_DCT_FBMC );
    
    %% To calculate the peak-to-average ratio, we truncate the signal at time periods of length 1/F.
    s_OFDM_reshaped                 = reshape( s_OFDM((CPOFDM.Implementation.ZeroGuardSamples+1):(end-CPOFDM.Implementation.ZeroGuardSamples)) , CPOFDM.Implementation.FFTSize+CPOFDM.Implementation.CyclicPrefix , CPOFDM.Nr.MCSymbols );
    s_OFDM_reshaped                 = s_OFDM_reshaped(CPOFDM.Implementation.CyclicPrefix+1:end ,:); % Not necessary

    s_OFDMnoCP_reshaped             = reshape( s_OFDMnoCP((OFDMnoCP.Implementation.ZeroGuardSamples+1):(end-OFDMnoCP.Implementation.ZeroGuardSamples)) , OFDMnoCP.Implementation.FFTSize+OFDMnoCP.Implementation.CyclicPrefix , OFDMnoCP.Nr.MCSymbols );
    s_OFDMnoCP_reshaped             = s_OFDMnoCP_reshaped(OFDMnoCP.Implementation.CyclicPrefix+1:end,:); % Not necessary
    
    s_DFT_OFDM_reshaped             = reshape( s_DFT_OFDM((CPOFDM.Implementation.ZeroGuardSamples+1):(end-CPOFDM.Implementation.ZeroGuardSamples)) , CPOFDM.Implementation.FFTSize+CPOFDM.Implementation.CyclicPrefix , CPOFDM.Nr.MCSymbols );
    s_DFT_OFDM_reshaped             = s_DFT_OFDM_reshaped(CPOFDM.Implementation.CyclicPrefix+1:end,:); % Not necessary

    s_FBMC_reshaped                 = reshape( s_FBMC((CPOFDM.Implementation.ZeroGuardSamples+1):(end-CPOFDM.Implementation.ZeroGuardSamples)) , CPOFDM.Implementation.FFTSize , ConFBMC.Nr.MCSymbols/2 );
    s_FBMC_reshaped(:,[1 end])      = [];
    
    
    s_DCT_FBMC_reshaped              = reshape( s_DCT_FBMC((CPOFDM.Implementation.ZeroGuardSamples+1):(end-CPOFDM.Implementation.ZeroGuardSamples)) , CPOFDM.Implementation.FFTSize , ConFBMC.Nr.MCSymbols/2 );
    s_DCT_FBMC_reshaped              = s_DCT_FBMC_reshaped(OFDMnoCP.Implementation.CyclicPrefix+1:end,:); % Not necessary
    
    %% Signal power, which should be 1.
    AveragePower_OFDM(i_rep)         = mean(mean( abs(s_OFDM_reshaped).^2 ));
    AveragePower_OFDMnoCP(i_rep)     = mean(mean( abs(s_OFDMnoCP_reshaped).^2 ));
    AveragePower_FBMC(i_rep)         = mean(mean( abs(s_FBMC_reshaped).^2 ));
    
    AveragePower_DFT_OFDM(i_rep)     = mean(mean( abs(s_DFT_OFDM_reshaped).^2 ));
    AveragePower_DCT_FBMC(i_rep)     = mean(mean( abs(s_DCT_FBMC_reshaped).^2 ));
  
    %% Peak power values at 1/F per time space.
    PeakPower_OFDM(i_rep,:)          = max( abs(s_OFDM_reshaped).^2 ,[],1);
    PeakPower_OFDMnoCP(i_rep,:)      = max( abs(s_OFDMnoCP_reshaped).^2 ,[],1);
    PeakPower_FBMC(i_rep,:)          = max( abs(s_FBMC_reshaped).^2 ,[],1);
    
    PeakPower_DFT_OFDM(i_rep,:)      = max( abs(s_DFT_OFDM_reshaped).^2 ,[],1);
    PeakPower_DCT_FBMC(i_rep,:)     = max( abs(s_DCT_FBMC_reshaped).^2 ,[],1);
    TimePassed = toc;
    if mod(i_rep,100)==0
        disp(['Realization ' int2str(i_rep) ' of ' int2str(NrRepetitions) '. Time left: ' ...
            int2str(TimePassed/i_rep*(NrRepetitions-i_rep)/60) 'min. Approx:' int2str(TimePassed/i_rep*(NrRepetitions-i_rep)/3600) 'h']);
    end
end


%% Plotting the peak-to-average power ratio
PAPR_OFDM         = 10*log10(PeakPower_OFDM(:));
PAPR_OFDMnoCP     = 10*log10(PeakPower_OFDMnoCP(:));
PAPR_FBMC         = 10*log10(PeakPower_FBMC(:));

PAPR_DFT_OFDM     = 10*log10(PeakPower_DFT_OFDM(:));
PAPR_DCT_FBMC      = 10*log10(PeakPower_DCT_FBMC(:));

[CCDF_OFDM , CCDF_OFDM_xAxis]                 = ecdf(PAPR_OFDM); CCDF_OFDM=1-CCDF_OFDM;
[CCDF_OFDMnoCP , CCDF_OFDMnoCP_xAxis]         = ecdf(PAPR_OFDMnoCP); CCDF_OFDMnoCP=1-CCDF_OFDMnoCP;
[CCDF_FBMC , CCDF_FBMC_xAxis]                 = ecdf(PAPR_FBMC); CCDF_FBMC=1-CCDF_FBMC;
[CCDF_DFT_OFDM , CCDF_DFT_OFDM_xAxis]         = ecdf(PAPR_DFT_OFDM); CCDF_DFT_OFDM=1-CCDF_DFT_OFDM;
[CCDF_DCT_FBMC , CCDF_DCT_FBMC_xAxis]       = ecdf(PAPR_DCT_FBMC); CCDF_DCT_FBMC=1-CCDF_DCT_FBMC;

Dim       = 300;
MaxMinAll = max([min(CCDF_OFDM_xAxis),min(CCDF_OFDMnoCP_xAxis),min(CCDF_FBMC_xAxis),min(CCDF_DCT_FBMC_xAxis)]);
MinMaxAll = min([max(CCDF_OFDM_xAxis),max(CCDF_OFDMnoCP_xAxis),max(CCDF_FBMC_xAxis),max(CCDF_DCT_FBMC_xAxis)]);

CCDF_OFDM_xAxis_LessPoints = linspace(MaxMinAll,MinMaxAll,Dim);
[x, index] = unique(CCDF_OFDM_xAxis);
CCDF_OFDM_LessPoints = 10.^interp1(x,log10(CCDF_OFDM(index)),CCDF_OFDM_xAxis_LessPoints);

CCDF_OFDMnoCP_xAxis_LessPoints = linspace(MaxMinAll,MinMaxAll,Dim);
[x, index] = unique(CCDF_OFDMnoCP_xAxis);
CCDF_OFDMnoCP_LessPoints = 10.^interp1(x,log10(CCDF_OFDMnoCP(index)),CCDF_OFDMnoCP_xAxis_LessPoints);

CCDF_FBMC_xAxis_LessPoints = linspace(MaxMinAll,MinMaxAll,Dim);
[x, index] = unique(CCDF_FBMC_xAxis);
CCDF_FBMC_LessPoints = 10.^interp1(x,log10(CCDF_FBMC(index)),CCDF_FBMC_xAxis_LessPoints);

CCDF_DCT_FBMC_xAxis_LessPoints = linspace(MaxMinAll,MinMaxAll,Dim);
[x, index] = unique(CCDF_DCT_FBMC_xAxis);
CCDF_DCT_FBMC_LessPoints = 10.^interp1(x,log10(CCDF_DCT_FBMC(index)),CCDF_DCT_FBMC_xAxis_LessPoints);

MaxMinAll_DFT = max([min(CCDF_DFT_OFDM_xAxis)]);
MinMaxAll_DFT = min([max(CCDF_DFT_OFDM_xAxis)]);

CCDF_DFT_OFDM_xAxis_LessPoints = linspace(MaxMinAll_DFT,MinMaxAll_DFT,100);
[x, index] = unique(CCDF_DFT_OFDM_xAxis);
CCDF_DFT_OFDM_LessPoints = 10.^interp1(x,log10(CCDF_DFT_OFDM(index)),CCDF_DFT_OFDM_xAxis_LessPoints);


%% Plot
figure(1);
markersize = 9;
LineWidth  = 1;
LineWidth1 = 0.7;
Markersize1 = 8;
semilogy(CCDF_OFDM_xAxis_LessPoints(10:45:Dim)    ,CCDF_OFDM_LessPoints(10:45:Dim),    's','Color',0.75*[0,1,0],'Markersize',markersize,'LineWidth',LineWidth);hold on; 
semilogy(CCDF_OFDMnoCP_xAxis_LessPoints(25:45:Dim),CCDF_OFDMnoCP_LessPoints(25:45:Dim),'x','Color',0.75*[0,0,1],'Markersize',markersize,'LineWidth',LineWidth);
semilogy(CCDF_FBMC_xAxis_LessPoints(40:45:Dim)    ,CCDF_FBMC_LessPoints(40:45:Dim),    'o','Color',0.75*[1,0,1],'Markersize',markersize,'LineWidth',LineWidth);
semilogy(CCDF_DFT_OFDM_xAxis_LessPoints(30:10:end),CCDF_DFT_OFDM_LessPoints(30:10:end),'d','Color',0.50*[0,0,1],'Markersize',markersize,'LineWidth',LineWidth+0.2); 
semilogy(CCDF_DCT_FBMC_xAxis_LessPoints(30:20:Dim) ,CCDF_DCT_FBMC_LessPoints(30:20:Dim), '*','Color',0.75*[1,0,0],'Markersize',markersize,'LineWidth',LineWidth);
semilogy(CCDF_FBMC_DFT_xAxis_LessPoints(30:20:100) ,CCDF_FBMC_DFT_LessPoints(30:20:100), '^','Color',0.75*[0,0,0],'Markersize',markersize,'LineWidth',LineWidth);
grid on; box on;
semilogy(CCDF_OFDM_xAxis_LessPoints    , CCDF_OFDM_LessPoints,      'Color',0.75*[0,1,0],'LineWidth',LineWidth);
semilogy(CCDF_OFDMnoCP_xAxis_LessPoints, CCDF_OFDMnoCP_LessPoints,  'Color',0.75*[0,0,1],'LineWidth',LineWidth);
semilogy(CCDF_FBMC_xAxis_LessPoints    , CCDF_FBMC_LessPoints,      'Color',0.75*[1,0,1],'LineWidth',LineWidth);
semilogy(CCDF_DFT_OFDM_xAxis_LessPoints, CCDF_DFT_OFDM_LessPoints,  'Color',0.60*[0,0,1],'LineWidth',LineWidth+0.2); 
semilogy(CCDF_DCT_FBMC_xAxis_LessPoints,CCDF_DCT_FBMC_LessPoints, 'Color',0.75*[1,0,0],'LineWidth',LineWidth);
semilogy(CCDF_FBMC_DFT_xAxis_LessPoints,CCDF_FBMC_DFT_LessPoints, 'Color',0.75*[0,0,0],'LineWidth',LineWidth);
ylim([10e-4 1]);
xlim([4 12]);
ylabel('CCDF');
xlabel('Peak to Average Power Ratio (dB)');
s1 = plot([nan nan],[nan nan],'s-','Color',0.85*[0,1,0],'Markersize',Markersize1,'LineWidth',LineWidth1);
s2 = plot([nan nan],[nan nan],'x-','Color',0.85*[0,0,1],'Markersize',Markersize1,'LineWidth',LineWidth1);
s3 = plot([nan nan],[nan nan],'o-','Color',0.85*[1,0,1],'Markersize',Markersize1,'LineWidth',LineWidth1);
s4 = plot([nan nan],[nan nan],'d-','Color',0.60*[0,0,1],'Markersize',Markersize1,'LineWidth',LineWidth1);
s5 = plot([nan nan],[nan nan],'*-','Color',0.85*[1,0,0],'Markersize',Markersize1,'LineWidth',LineWidth1);
s6 = plot([nan nan],[nan nan],'^-','Color',0.85*[0,0,0],'Markersize',Markersize1,'LineWidth',LineWidth1);
legend([s1 s2 s3 s4 s5 s6],{'OFDM (with CP)','OFDM (no CP)','Conventional FBMC','SC-FDMA','Pruned DCT-p-FBMC','Pruned DFT-s-FBMC'},'Location','Southwest');
title('16-QAM, $L=64$','Interpreter','latex');
set(gca,'FontName','Times New Roman','FontSize',12);

