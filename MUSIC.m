function [degrees,delay] = MUSIC(Rsig,M,N,d,lam,subcarr)
    Rsig = Rsig.';
    Rsig = Rsig(:);
%   MUSIC
    %   Auto-correlation matrix
     Rxx = Rsig*Rsig';
     
%      % eigen decomposition technique #1   
%      [eigVector, eigValue] = eig(Rxx);
%      eigValue = ((diag(eigValue)));
%      % signal and noise eigenvalue threshold
%      eigNoiseThr = 1e-6;
%      % signal subspace
%      eigVectSignal = eigVector(:, eigValue > eigNoiseThr);
%      MI=eigVectSignal*inv(eigVectSignal'*eigVectSignal)*eigVectSignal';
%      % noise subspace
%      MI=eye(size(MI))-MI;


[E, D] = eig(Rxx); %eigenvalues and eigenvectors
D = diag(D);
[~,I] = sort(D,'descend');

E = E(:,I); 
En = E(:,M+1:end);% eigenvectors of noise subspace

Pn = En * inv(En'*En) * En';
   
% eigen decomposition technique #1
    % can't get accurate results in 2D 
%  Auto-correlation matrix
%     Rxx = Rsig*Rsig'/312;
% %  eigenvalue decomposition
%     [Vi,Li] = eig(Rxx);
% %  sort in descending order
%     [L,I] = sort(diag(Li),'descend');
%     V = Vi(:,I);
% %   Signal Subspace
%     Ps = V(:,1:M)*V(:,1:M)';
% %   Noise Subspace
%     Pn = V(:,1+M:N)*V(:,1+M:N)';

    theta1=[0:1:90];
    tau=[0:1e-9:5e-8];
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
                A(:,l) = exp(2*pi*-1i*(l*deltaF)*tau(j));
           end
           A1 = kron(B,A.');
           PMUSIC(i,j)= 1/(sum(abs(A1'*Pn*A1)));
       end
    end
 
    [C,I] = max(PMUSIC(:));
    [I1,I2] = ind2sub(size(PMUSIC),I);
    deg = I1-1;
    del = I2-1;
%     figure(1);
%     plot(theta1,10*log10(PMUSIC));
%     title('MUSIC spectrum');
%     xlabel('Angle [degrees]');
%     ylabel('PMUSIC');   
    figure(1);
    [X,Y] = meshgrid(theta1,tau);
    surf(X,Y,10*log10(PMUSIC)')
    shading interp 
    colorbar    
    title('MUSIC spectrum');
    xlabel('Angle [degrees]');
    ylabel('Delay (s)');
    zlabel('PMUSIC');
    
    degrees = deg;
    delay = del;
end
