clc;
clear all;

numFrames = 1; % 10 ms frames
%snr = 1; % SNR in dB
fc = 3e9;


% Create UE/carrier configuration
ue = nrCarrierConfig;
ue.NSizeGrid = 52;          % Bandwidth in number of resource blocks (52RBs at 15kHz SCS for 10MHz BW)
ue.SubcarrierSpacing = 15;  % 15, 30, 60, 120, 240 (kHz)
ue.CyclicPrefix = 'Normal'; % 'Normal' or 'Extended'

nTxAnts = 1;  % Number of transmit antennas (1,2,4)
nRxAnts = 1;  % Number of receive antennas
nLayers = min(nTxAnts,nRxAnts);

srs = nrSRSConfig;
srs.NumSRSSymbols = 1;          % Number of OFDM symbols allocated per slot (1,2,4)
srs.SymbolStart = 8;            % Starting OFDM symbol within a slot
srs.NumSRSPorts = 1;            % Number of SRS antenna ports (1,2,4).
srs.FrequencyStart = 0;         % Frequency position of the SRS in BWP in RBs
srs.NRRC = 0;                   % Additional offset from FreqStart specified in blocks of 4 PRBs (0...67)
srs.CSRS = 14;                  % Bandwidth configuration C_SRS (0...63). It controls the allocated bandwidth to the SRS
srs.BSRS = 0;                   % Bandwidth configuration B_SRS (0...3). It controls the allocated bandwidth to the SRS
srs.BHop = 0;                   % Frequency hopping configuration (0...3). Set BHop < BSRS to enable frequency hopping
srs.KTC = 2;                    % Comb number (2,4). Frequency density in subcarriers
srs.Repetition = 2;             % Repetition (1,2,4). It disables frequency hopping in blocks of |Repetition| symbols
srs.SRSPeriod = [2 0];          % Periodicity and offset in slots. SRSPeriod(2) must be < SRSPeriod(1)
srs.ResourceType = 'aperiodic';  % Resource type ('periodic', 'semi-persistent','aperiodic'). Use 'aperiodic' to disable inter-slot frequency hopping
srs.SRSPositioning = 1;
srs.NSRSID = 0;

srs2 = nrSRSConfig;
srs2.NumSRSSymbols = 1;          % Number of OFDM symbols allocated per slot (1,2,4)
srs2.SymbolStart = 8;            % Starting OFDM symbol within a slot
srs2.NumSRSPorts = 1;            % Number of SRS antenna ports (1,2,4).
srs2.FrequencyStart = 0;         % Frequency position of the SRS in BWP in RBs
srs2.NRRC = 0;                   % Additional offset from FreqStart specified in blocks of 4 PRBs (0...67)
srs2.CSRS = 14;                  % Bandwidth configuration C_SRS (0...63). It controls the allocated bandwidth to the SRS
srs2.BSRS = 0;                   % Bandwidth configuration B_SRS (0...3). It controls the allocated bandwidth to the SRS
srs2.BHop = 0;                   % Frequency hopping configuration (0...3). Set BHop < BSRS to enable frequency hopping
srs2.KTC = 2;                    % Comb number (2,4). Frequency density in subcarriers
srs2.Repetition = 2;             % Repetition (1,2,4). It disables frequency hopping in blocks of |Repetition| symbols
srs2.SRSPeriod = [2 0];          % Periodicity and offset in slots. SRSPeriod(2) must be < SRSPeriod(1)
srs2.ResourceType = 'aperiodic';  % Resource type ('periodic', 'semi-persistent','aperiodic'). Use 'aperiodic' to disable inter-slot frequency hopping
srs2.SRSPositioning = 1;
srs2.NSRSID = 56;

channel = nrTDLChannel;
channel.DelayProfile = 'Custom';
channel.FadingDistribution = 'Rician';
channel.KFactorFirstTap = 13.3;
channel.PathDelays = [0];
channel.AveragePathGains = [0];
channel.NumTransmitAntennas = 1;
channel.NumReceiveAntennas = 1;
channel.SampleRate = 15360000;
channel.NormalizePathGains = false;

channel2 = nrTDLChannel;
channel2.DelayProfile = 'Custom';
channel2.FadingDistribution = 'Rayleigh';
channel2.PathDelays = [0];
channel2.AveragePathGains = [-3];
channel2.NumTransmitAntennas = 1;
channel2.NumReceiveAntennas = 1;
channel2.SampleRate = 15360000;
channel2.NormalizePathGains = false;

% Number of slots to simulate
numSlots = numFrames*ue.SlotsPerFrame;

% Total number of subcarriers and symbols per slot
K = ue.NSizeGrid * 12;
L = ue.SymbolsPerSlot;

%             [srsIndices2,srsIndInfo2] = nrSRSIndices(ue,srs2);
%             srsSymbols2 = nrSRS(ue,srs2);
% 
%             % Create a slot-wise resource grid empty grid and map SRS symbols
%             txGrid2 = nrResourceGrid(ue,nTxAnts);
%             txGrid2(srsIndices2) = srsSymbols2;
% 
%             % OFDM Modulation
%             [txWaveform2,waveformInfo1] = nrOFDMModulate(ue,txGrid2);

for nSlot = 0:numSlots-1
            % Update slot counter
            ue.NSlot = nSlot;

            % Generate SRS and map to slot grid
            [srsIndices,srsIndInfo] = nrSRSIndices(ue,srs);
            srsSymbols = nrSRS(ue,srs);

            % Create a slot-wise resource grid empty grid and map SRS symbols
            txGrid = nrResourceGrid(ue,nTxAnts);
            txGrid(srsIndices) = srsSymbols;

            % Determine if the slot contains SRS
            isSRSSlot= ~isempty(srsSymbols);
            % OFDM Modulation
            [txWaveform,waveformInfo] = nrOFDMModulate(ue,txGrid);
            
            c = physconst('LightSpeed');
         %  lambda (λ)
            wavelength = c/fc; 
         %  Element Spacing
            d = 0.5*wavelength;
         %  Number of Elements   
            N = 12;
            
            theta = [20;60];
         %  Number of sources
            M = length(theta);

            for k = 1:M
                SteeringVector(:, k) = exp(1i*2*pi*(d/wavelength)*sind(theta(k))*[0:N-1]'); 
            end
            
%              Transmission through channel
            [rxWaveform1,pathGains] = channel(txWaveform);
            [rxWaveform2,pathGains2] = channel2(txWaveform);
            chanFilt = comm.ChannelFilter('SampleRate', 15360000,'PathDelays',[1e-6]);
            rxWaveform1(:,2) = chanFilt(rxWaveform2,pathGains2);

            x1 = SteeringVector*rxWaveform1.';
            x = x1;
            x = x.';
 
%           x = awgn(x,10,'measured');
          % OFDM Demodulation
            rxGrid = nrOFDMDemodulate(ue,x);
            subcarr = length(rxGrid(:,1));
          % Check Normal/Extended Cyclic Prefix
            if strcmp(ue.CyclicPrefix,'normal')
                %   Repeat for SRS Number in a single slot
                for i=9:9
                    rxGridEst=zeros(K,N);
                    C = permute(rxGrid,[1 3 2]);
                    rxGridEst(:,:) = C(:,:,i);
                    zer = zeros(1,subcarr/2);
                    srsSym = srsSymbols(1:subcarr/2).';
%                     srscsi = [srsSym; zer];
%                     srscsi = reshape(srscsi,1,[]);
                    rxGridEst = rxGridEst.';
%                      for i = 1:N
%                          CSI(:,i)= rxGridEst(i,:)./srscsi;
%                      end
%                      for i = 1:312
%                          CSI1(i,:) = CSI(2*i-1,:);
%                      end
                    conjSRS = conj(srsSym);
                    rxGridEst1(1:N,1:subcarr/2) = zeros;
                    for ii = 1:subcarr/2
                        rxGridEst1(1:N,ii*2-1) = conjSRS(ii)*rxGridEst(1:N,ii*2-1);
                    end
                    for iii = 1:subcarr/2
                        rxGridEst2(:,iii) = rxGridEst1(:,2*iii-1);
                    end
                    [degrees] = MUSIC(rxGridEst2,M,N,d,wavelength,subcarr);
                 end 
            elseif strcmp(ue.CyclicPrefix,'extended')
                for i=6:6+srs.NumSRSSymbols-1
                    rxGridEst=zeros(K,N);
                    C = permute(rxGrid,[1 3 2]);
                    rxGridEst(:,:) = C(:,:,i);
                    rxGridEst = rxGridEst;
                    rxGridEst = rxGridEst.';
                    [degrees] = MUSIC(rxGridEst,M,N,d,wavelength);
                end
            else 
                disp('Unexpected Cyclic Prefix Value')
            end
end
