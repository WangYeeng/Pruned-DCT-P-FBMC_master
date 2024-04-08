%==========================================================================
% Author: Ying Wang
% Function: Transmitted power and power spectral density of different multicarrier modulation systems.
% e-mail address: wangyingstu@163.com
%========================================================================== 
clc; clear;close all;
addpath("Channel\","Coding\","MultiCarrierModulation\");
%% Parameters
NrRepetitions       = 10000;                       % Monte Carlo Repetitions
L                   = 64;                          % Number of subcarriers
QAM_ModulationOrder = 16;                          % Modulation order, 4, 16, 64,...
SubcarrierSpacing   = 15e3;                        % Subcarrier spacing (15kHz)
K                   = 30;                          % FBMC No. of symbols
K_OFDM              = K/2-1;                       % OFDM No. of symbols (same as in LTE)
CP_Length           = 1/SubcarrierSpacing/14;      % LTE CP Length £¨seconds£©
SamplingRate        = 15e3*L*21;                   % Sampling rate 

%% FBMC 
FBMC = FBMC(...
    L,...                                           % Number of subcarriers
    K,...                                           % FBMC No. of symbols
    SubcarrierSpacing,...                           % Subcarrier spacing (Hz)
    SamplingRate,...                                % Sampling rate (Samples/s)
    0,...                                           % Intermediate frequency of the first subcarrier (Hz).
    false,...                                       % Transmit real valued signal
    'PHYDYAS-OQAM',...                              % Prototype filter and OQAM or QAM. 
    4, ...                                          % Overlapping factor
    0, ...                                          % Initial phase shift
    true ...                                        % Efficient IFFT implementation, true or false
    );
G = FBMC.GetTXMatrix;
%% OFDM 
ZeroGuardTimeLength = ((FBMC.Nr.SamplesTotal-(round(SamplingRate/SubcarrierSpacing)+CP_Length*SamplingRate)*K_OFDM)/2)/SamplingRate;
OFDM = OFDM(...
    L,...                                           % Number of subcarriers
    K_OFDM,...                                      % OFDM No. of symbols
    SubcarrierSpacing,...                           % Subcarrier spacing (Hz)
    SamplingRate,...                                % Sampling rate (Samples/s)
    0,...                                           % Intermediate frequency of the first subcarrier (Hz).
    false,...                                       % Transmit real valued signal
    CP_Length, ...                                  % Cycle prefix length (s) 1/SubcarrierSpacing/(K/2-1)
    ZeroGuardTimeLength ...                         % Length of the guard time (s)
    );

%% Check the number of samples
if  OFDM.Nr.SamplesTotal~=FBMC.Nr.SamplesTotal
   error('Total number of samples must be the same for OFDM and FBMC.');
end
N = OFDM.Nr.SamplesTotal;

%% Symbol constellations
QAM = SignalConstellation(QAM_ModulationOrder,'QAM');
PAM = SignalConstellation(sqrt(QAM_ModulationOrder),'PAM');
%% Precoding matrix
NrSubcarriers = FBMC.Nr.Subcarriers;
Domain        = 'Time';  % Setup Domain
Mathod        = 'DCT';   % Adopting the DCT method
Kall          = K;       % No. of symbols
K             = Kall;    % FBMC No. of symbols
L             = NrSubcarriers;

phi = GetCodingMatrix(...
    Domain,...          % Determine the domain ('Frequency', 'Time', 'DeEigen')
    K,...               % FBMC No. of symbols
    L,...               % Number of subcarriers
    Mathod,...          % Precoding matrix form
    FBMC...             % The functions are embedded in the FBMC system.
    );

%% Generate DFT Matrix
DFTMatrix               = fft(eye(NrSubcarriers))/sqrt(NrSubcarriers);
if size(phi,2)~=(NrSubcarriers)/2
   error('Size should be "(Number of Subcarriers)/2"');
end
%% Calculating Transmitted Power
Pt_OFDM = zeros(N,1);
Pt_FBMC = zeros(N,1);
Pt_Prund_DCT_FBMC = zeros(N,1);
Pt_SC_FDMA        = zeros(N,1);

Pf_OFDM = zeros(N,1);
Pf_FBMC = zeros(N,1);
Pf_Prund_DCT_FBMC = zeros(N,1);
Pf_SC_FDMA        = zeros(N,1);
tic;
for i_rep = 1:NrRepetitions
    %% Generate Bitstream
    BinaryDataStream           = randi([0 1],(K_OFDM+1)*L*log2(QAM.ModulationOrder),1);
    BinaryDataStream_FBMC      = randi([0 1],K*L*log2(PAM.ModulationOrder),1);
    BinaryDataStream_OFDM      = randi([0 1],K_OFDM*L*log2(QAM.ModulationOrder),1);          
    % Bit-to-symbol mapping
    x_OFDM           = reshape(QAM.Bit2Symbol(BinaryDataStream_OFDM),L,K_OFDM);
    x_FBMC           = reshape(PAM.Bit2Symbol(BinaryDataStream_FBMC),L,K);
    x_Prund_DCT_FBMC = reshape(QAM.Bit2Symbol(BinaryDataStream),L/2,K);
    % Generating transmitted signals in the time domain
    s_OFDM           = OFDM.Modulation(x_OFDM);
    s_FBMC           = G*reshape(x_FBMC,L*K,1);
    s_Prund_DCT_FBMC = G*reshape(phi*x_Prund_DCT_FBMC,L*K,1)*sqrt(2);
    s_SC_FDMA        = OFDM.Modulation(DFTMatrix*x_OFDM);
    %% Calculate power
    Pt_OFDM           = Pt_OFDM + abs(s_OFDM).^2;
    Pt_FBMC           = Pt_FBMC + abs(s_FBMC).^2;
    Pt_Prund_DCT_FBMC = Pt_Prund_DCT_FBMC + abs(s_Prund_DCT_FBMC).^2;
    Pt_SC_FDMA        = Pt_SC_FDMA + abs(s_SC_FDMA).^2;
    %% Calculate the power spectral density
    Pf_OFDM           = Pf_OFDM + abs(fft((full(s_OFDM)))).^2;
    Pf_FBMC           = Pf_FBMC + abs(fft(s_FBMC)).^2;
    Pf_Prund_DCT_FBMC = Pf_Prund_DCT_FBMC + abs(fft((full(s_Prund_DCT_FBMC)))).^2;
    Pf_SC_FDMA        = Pf_SC_FDMA + abs(fft(s_SC_FDMA)).^2;
    
    TimePassed = toc;
    if mod(i_rep,1000)==0
        disp(['Realization ' int2str(i_rep) ' of ' int2str(NrRepetitions) '. Time left: ' int2str(TimePassed/i_rep*(NrRepetitions-i_rep)/60) 'minutes = ' int2str(TimePassed/i_rep*(NrRepetitions-i_rep)/60/60) 'hours']);
    end
end
%% Calculate the power and power spectral density
% % OFDM
P_OFDM     = Pt_OFDM/NrRepetitions;
PSD_OFDM   = Pf_OFDM/max(Pf_OFDM);
% % FBMC
P_FBMC     = Pt_FBMC/NrRepetitions;
PSD_FBMC   = Pf_FBMC/max(Pf_FBMC);
% % Prund_DCT_FBMC
P_Prund_DCT_FBMC     = Pt_Prund_DCT_FBMC/NrRepetitions;
PSD_Prund_DCT_FBMC   = Pf_Prund_DCT_FBMC/max(Pf_Prund_DCT_FBMC);
% % SC-FDMA
P_SC_FDMA     = Pt_SC_FDMA/NrRepetitions;
PSD_SC_FDMA   = Pf_SC_FDMA /max(Pf_SC_FDMA);

%% dB
PSD_OFDM             = 10*log10(PSD_OFDM);
PSD_FBMC             = 10*log10(PSD_FBMC);
PSD_Prund_DCT_FBMC   = 10*log10(PSD_Prund_DCT_FBMC);
PSD_SC_FDMA          = 10*log10(PSD_SC_FDMA);
%% plot
LineWidth = 1.5;
f = (0:OFDM.Nr.SamplesTotal-1)*OFDM.PHY.SamplingRate/(OFDM.Nr.SamplesTotal);
f_Normalized = f/FBMC.PHY.SubcarrierSpacing - L+1;
t = (0:OFDM.Nr.SamplesTotal-1)*OFDM.PHY.dt;
t_Normalized = t*FBMC.PHY.SubcarrierSpacing;

figure(1)% power
S1 = plot(t_Normalized,P_Prund_DCT_FBMC ,'Color',0.85*[1 0 0],  'LineWidth',LineWidth);hold on;grid on;
S2 = plot(t_Normalized,P_OFDM ,          'Color',0.85*[1 0 1],  'LineWidth',LineWidth);
S3 = plot(t_Normalized,P_SC_FDMA,        'Color',0.85*[0 0 1],  'LineWidth',LineWidth);
S4 = plot(t_Normalized,P_FBMC ,          'Color',0.85*[0 1 0],  'LineWidth',LineWidth);
ylim([0 1.2]);
xlabel('Normalized time $tF$','Interpreter','latex');
ylabel('Power');
legend([S2,S3,S4,S1],{'CP-OFDM', 'SC-FDMA','Conventional FBMC','Prund DCT-p-FBMC'},'Location','southeast')      
set(gca,'YTick',[0 1],'FontName','Times New Roman','FontSize',12,'GridLineStyle','-.');

figure(2)% power spectral density
plot(f_Normalized,PSD_OFDM ,               'Color',0.85*[1 0 1], 'LineWidth',LineWidth);hold on;grid on;
plot(f_Normalized,PSD_SC_FDMA,        '--','Color',0.85*[0 0 1], 'LineWidth',LineWidth);
plot(f_Normalized,PSD_FBMC ,               'Color',0.85*[0 1 0], 'LineWidth',LineWidth);
plot(f_Normalized,PSD_Prund_DCT_FBMC ,'--','Color',0.85*[1 0 0], 'LineWidth',LineWidth);
ylim([-150 5]);
xlim([-5 10]);
xlabel('Normalized frequency $f/F$','Interpreter','latex');
ylabel('Power spectral density (dB)');
legend({'OFDM (with CP)','SC-FDMA','Conventional FBMC','Pruned DCT-p-FBMC'},'Location','SouthWest')      
set(gca,'YTick',[-200 -150 -100 -50 -0],'FontName','Times New Roman','FontSize',12,'GridLineStyle','-.');
