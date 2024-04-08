function C = GetCodingMatrix(...
    Domain,...          % Determine the domain ('Frequency', 'Time', 'DeEigen')
    K,...               % FBMC No. of symbols
    L,...               % Number of subcarriers
    Mathod,...          % Precoding matrix form
    FBMC...
    )
% Function: Generate the initial matrix of the precoding matrix.
  switch Domain 
      case 'Frequency'
     FBMC.SetNrMCSymbols(K); 
     C = FBMC.GetPrecodingMatrixForQAMinOQAM(Domain,2,Mathod);           
      case 'Time'
     FBMC.SetNrSubcarriers(L);
     FBMC.SetNrMCSymbols(1);
     C = FBMC.GetPrecodingMatrixForQAMinOQAM(Domain,1,Mathod);
      case 'DeEigen'
     FBMC.SetNrSubcarriers(L);
     FBMC.SetNrMCSymbols(1);  
     D = FBMC.GetFBMCMatrix;
     [U,V] = eig(D);
     C = U*sqrt(V(:,1:(FBMC.Nr.Subcarriers*FBMC.Nr.MCSymbols)/2));
      case 'Time-DFT'
     FBMC.SetNrMCSymbols(1);
     FBMC.SetNrSubcarriers(L);
     DFTMatrix   = fft(eye(L))/sqrt(L);
     D           = FBMC.GetFBMCMatrix;
     a           = abs(diag(DFTMatrix'*D*DFTMatrix));
     a           = a+randn(size(a))*10e-12;
     a_Tilde     = sort(a,'descend');
     alpha       = a_Tilde((L)/2);
     Index_Tilde = (a>=alpha);
     DFTMatrix_Tilde   = DFTMatrix(:,Index_Tilde) ;
     b_Tilde     = sqrt(2./(a(Index_Tilde)));
     C           = DFTMatrix_Tilde*diag(b_Tilde);
  end
end

