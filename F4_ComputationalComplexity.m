% =========================================================================   
% Author: Ying Wang
% Function: Calculate the computational complexity of the Pruned DCT-P-FBMC system with the number of multiplications.
% e-mail address: wangyingstu@163.com
% =========================================================================   
clc; clear; close all;
% Parameters
N_FFT =  0:16:4056;                     % FFT size for OFDM and FBMC
NoverL = [2 4 32 128];                  % Factor: "N_FFT/L"                      
O = 4;                                  % Overlapping factor
% Start calculation
for i_NoverL = 1:length(NoverL)
    L = N_FFT/NoverL(i_NoverL);
    Pruned_DCT_FBMC = 2*(L/2+L.*log2(L/2)+N_FFT.*log2(N_FFT)+O*N_FFT);
    SC_FDMA = L.*log2(L)+N_FFT.*log2(N_FFT);
    
    Complexity_OFDM(:,i_NoverL)      = SC_FDMA;   
    Complexity_HOLD_FBMC(:,i_NoverL) = Pruned_DCT_FBMC;
    ComplexityIncrease(:,i_NoverL)   = Pruned_DCT_FBMC./SC_FDMA;
   
end

% Plot the computational complexity relative to the reference
LineWidth = 1;
figure(1);
plot(N_FFT,ComplexityIncrease(:,1),'Color',0.75*[1 0 0],'LineWidth',LineWidth);hold on;grid on;
plot(N_FFT,ComplexityIncrease(:,2),'Color',0.75*[0 1 0],'LineWidth',LineWidth);hold on;grid on;
plot(N_FFT,ComplexityIncrease(:,3),'Color',0.75*[0 0 1],'LineWidth',LineWidth);hold on;grid on;
plot(N_FFT,ComplexityIncrease(:,4),'Color',0.75*[1 0 1],'LineWidth',LineWidth);hold on;grid on;
xlabel('$N_{FFT}$ for Pruned DCT-p-FBMC','interpreter','latex'); 
ylabel('Relative Complexity');
title('Pruned DCT-p-FBMC vs SC-FDMA');
legend([repmat('$N_{FFT}/L = $',length(NoverL),1) int2str(NoverL')],'interpreter','latex');
xlim([min(N_FFT) max(N_FFT)]);
ylim([2.35 3.15]);
set(gca,'Ytick',[2.35 2.55 2.75 2.95 3.15],'FontName','Times New Roman','FontSize',12,'GridLineStyle','--');

