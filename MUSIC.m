function [degrees,delay] = MUSIC(Rsig,M,N,d,lam,subcarr)
    Rsig = reshape(Rsig,[],1);
%   MUSIC
    %   Auto-correlation matrix
     Rxx = Rsig*Rsig'/subcarr/2;
     
     % eigen decomposition technique #1   
     [eigVector, eigValue] = eig(Rxx);
     eigValue = (abs(diag(eigValue)));
     % signal and noise eigenvalue threshold
     eigNoiseThr = 0.01;
     % signal subspace
     eigVectSignal = eigVector(:, eigValue > eigNoiseThr);
     MI=eigVectSignal*inv(eigVectSignal'*eigVectSignal)*eigVectSignal';
     % noise subspace
     MI=eye(size(MI))-MI;

    % eigen decomposition technique #1
    % can't get accurate results in 2D 
%   Auto-correlation matrix
%     Rxx = Rsig*Rsig'/312;
%   eigenvalue decomposition
%     [Vi,Li] = eig(Rxx);
%   sort in descending order
%     [L,I] = sort(diag(Li),'descend');
%     V = Vi(:,I);
% %   Signal Subspace
%     Ps = V(:,1:M)*V(:,1:M)';
% %   Noise Subspace
%     Pn = V(:,1+M:N)*V(:,1+M:N)';
 
    theta1=[0:1:90];
    tau=[0:1e-7:5e-6];
    deltaF = 30e3;
    % The MUSIC spectrum

    for i=1:length(theta1)
        phi1=2*pi*(d/lam)*sin(theta1(i)*pi/180);
        B =zeros([N 1]);
        for k=1:N
            B(k,1)= (exp((k-1)*1i*phi1));
        end
        for j=1:length(tau)
           for l = 1:subcarr/2
                A(:,l) = exp(2*pi*-1i*l*deltaF*tau(j));
           end
           A1 = kron(A.',B);
           PMUSIC(i,j)= (subcarr/2)/(sum(abs(A1'*MI*A1)));
       end
    end
    
