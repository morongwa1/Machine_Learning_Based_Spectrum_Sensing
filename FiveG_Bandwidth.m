% clear;
% clc;
addpath('functions')

ifWrite2File = 1;                                                           % 1 if data is to be written to file, 0 otherwise
numfiles = 2;                                                               % 1 file for ML trainig, 1 for testing

numUsers1 = 5;                                                              % number of users/number of different applications for allocation
numUsers2 = 5;   
numUsers3 = 5;
numUsers4 = 5;  
numUsers5 = 5;
numUsers6 = 5;
numUsers7 = 5;


polar = 0;                                                                  % bits polarity
modType = 'QPSK';                                                           % type of modulation used
numBitsPerSymb = NumberBits(modType);                                       % number of bits per modulation symbol

% OFDM variables
FFTLength = 4096;                                                           % numberc of fft samples used in ofdm
numSubcarriers = 12;                                                        % number of subcarriers per LTE resource block
numCPSamples = 72;                                                          % number of samples of OFDM symbol used for cyclic prefix
timeSlotsNum = 80;                                                          % total number of time slots used in program
resBlocksNum = 150;                                                          % number of resource blocks for chosen band (for example resBlocksNum = 50 for 10 MHz LTE channel)
ofdmSymbPerTimeSlot = 1;                                                    % number of OFDM sybmols per one time slot

resourceAlloc = zeros(resBlocksNum+1, timeSlotsNum+1);                      % table of resources assigned to users (for visualization purposes)
resourceAllocShort1 = zeros(numUsers1,4);                                     % table of resources assigned to users (for chcecking spectrum occupancy)
                                                                            %   columns: first freq, last freq, first time slot, last time slot
                                                                            %   rows: users
                                                                            
resourceAllocShort2 = zeros(numUsers2,4);
resourceAllocShort3 = zeros(numUsers3,4);
resourceAllocShort4 = zeros(numUsers4,4);
resourceAllocShort5 = zeros(numUsers5,4);
resourceAllocShort6 = zeros(numUsers6,4);
resourceAllocShort7 = zeros(numUsers7,4);

detectedRes = zeros(resBlocksNum*numSubcarriers, timeSlotsNum);
allocPossible1 = zeros(1, numUsers1);                                          % determines if resources allocation is possible at a given time slot

allocPossible2 = zeros(1, numUsers2); 
allocPossible3 = zeros(1, numUsers3); 
allocPossible4 = zeros(1, numUsers4); 
allocPossible5 = zeros(1, numUsers5); 
allocPossible6 = zeros(1, numUsers6);
allocPossible7 = zeros(1, numUsers7); 



LTETransmData = zeros(1,timeSlotsNum*ofdmSymbPerTimeSlot...
    *(FFTLength+numCPSamples));                                              % LTE data transmitted in time domain

                                                           

Pf = 0.1;                                                                    % Probability of false alarm - used in ED threshold calculation

%% channel parameters

% SNR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SNR = -40:1:20;
SNR = -40;

% Fading %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fs = 14.336e6;                                                      % Sample rate (Hz)
pathDelays = [0 30e-9 70e-9 90e-9 110e-9 190e-9 410e-9];            % Path delays (s)
pathPower = [0 -1 -2 -3 -8 -17.2 -20.8];                            % Path power (dB)
fD = 0;                                                             % Maximum Doppler shift (Hz)
numSamples = length(LTETransmData);                                 % Number of samples


%% signal generation
for num_files = 1:1
disp(['      number_of_iteration = ' num2str(num_files)])
        for fileIdx=1:numfiles
            LTETransmData = zeros(1,timeSlotsNum*ofdmSymbPerTimeSlot*...
                (FFTLength+numCPSamples));
            resourceAlloc = zeros(resBlocksNum+1, timeSlotsNum+1);
            
            resourceAllocShort1 = zeros(numUsers1,4);
            allocPossible1 = zeros(1, numUsers1);
            
            resourceAllocShort2 = zeros(numUsers2,4);
            allocPossible2 = zeros(1, numUsers2);
            
            resourceAllocShort3 = zeros(numUsers3,4);
            allocPossible3 = zeros(1, numUsers3);
            
            resourceAllocShort4 = zeros(numUsers4,4);
            allocPossible4 = zeros(1, numUsers4);

            resourceAllocShort5 = zeros(numUsers5,4);
            allocPossible5 = zeros(1, numUsers5);
            
            
            resourceAllocShort6 = zeros(numUsers6,4);
            allocPossible6 = zeros(1, numUsers6);  
            
            resourceAllocShort7 = zeros(numUsers7,4);
            allocPossible7 = zeros(1, numUsers7);   
            
            detect2 = zeros(1, length(SNR));
            falseAl2 = zeros(1, length(SNR));

            % fading channel impulse response generation
            rchan = comm.RayleighChannel('SampleRate',fs, ...
                'PathDelays',pathDelays,'AveragePathGains',pathPower, ...
                'MaximumDopplerShift',fD,'FadingTechnique','Sum of sinusoids');

            disp('   resource allocation');

            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % LTE-RB generation
            for i = 1:numUsers1
                disp(['Numerology 0 User: ', num2str(i)]);
                while allocPossible1(i)==0
                    numTSlots = randomParam(16, 0, true, 1, timeSlotsNum);                      % number of time slots (time domain)
                    firstTimeSlot = randomParam(32, 32, true, 1, timeSlotsNum-numTSlots+1);     % first time slot index
                    numResourceBlocks = randomParam(4, 0, true, 1, 50);               % number of resource blocks (frequency domain)
                    firstRB = randomParam(25, 20, true, 1, 50-numResourceBlocks+1);    % first resource block index
                    lastTimeSlot = firstTimeSlot + numTSlots - 1;                               % last time slot index
                    lastRB = firstRB + numResourceBlocks - 1;                                    % last resource block index

                    numBits = ofdmSymbPerTimeSlot*numSubcarriers*numBitsPerSymb*...
                        numTSlots*numResourceBlocks;                                            % number of data bits to generate in the beginning of program

                    [resourceAllocShort1, allocPossible1(i)] = CheckResourcesAllocation(firstRB, ...
                        lastRB, firstTimeSlot, lastTimeSlot, resourceAllocShort1, ...
                        i, timeSlotsNum, 50);
                end

                if allocPossible1(i)
                    resourceAlloc(firstRB:lastRB, firstTimeSlot:lastTimeSlot) = ...
                        i*ones(numResourceBlocks, numTSlots);

                    % calculating frequency subcarriers indexes for resource blocks allocated for given user
                    if firstRB > lastRB
                        disp([firstRB lastRB])
                    end
                    [firstSubcarrier, lastSubcarrier] = SubcarrierAllocation([firstRB lastRB], ...
                        FFTLength, numSubcarriers, 50);

                    %transmiting ofdm data
                    ofdmData = OneLTEUserTransm(numBits, polar, modType, ...
                        numTSlots, FFTLength, numCPSamples, firstSubcarrier,...
                        ofdmSymbPerTimeSlot);
                    ofdmVectL = numTSlots*ofdmSymbPerTimeSlot*length(ofdmData(:,1));
                    ofdmVect = reshape(ofdmData, [1,ofdmVectL]);
                    LTETransmData = transmissionLTE(ofdmVect, firstTimeSlot, numTSlots, ...
                        LTETransmData, FFTLength, numCPSamples, ofdmSymbPerTimeSlot);
                end

            end
            
            
            
            

            for i = 1:numUsers2
                disp(['Numerology 1 User: ', num2str(i)]);
                while allocPossible2(i)==0
                    numTSlots = randomParam(4, 0, true, 1, timeSlotsNum);                      % number of time slots (time domain)
                    firstTimeSlot = randomParam(29, 32, true, 1, timeSlotsNum-numTSlots+1);     % first time slot index
                    numResourceBlocks = randomParam(8, 0, true, 1, 100);               % number of resource blocks (frequency domain)
                    firstRB = randomParam(50+50+25, 15, true, 100, 150-numResourceBlocks+1);    % first resource block index
                    lastTimeSlot = firstTimeSlot + numTSlots - 1;                               % last time slot index
                    lastRB = firstRB + numResourceBlocks - 1;                                    % last resource block index

                    numBits = ofdmSymbPerTimeSlot*numSubcarriers*numBitsPerSymb*...
                        numTSlots*numResourceBlocks;                                            % number of data bits to generate in the beginning of program

                    [resourceAllocShort2, allocPossible2(i)] = CheckResourcesAllocation(firstRB, ...
                        lastRB, firstTimeSlot, lastTimeSlot, resourceAllocShort2, ...
                        i, timeSlotsNum, 150);
                end

                if allocPossible2(i)
                    resourceAlloc(firstRB:lastRB, firstTimeSlot:lastTimeSlot) = ...
                        i*ones(numResourceBlocks, numTSlots);

                    % calculating frequency subcarriers indexes for resource blocks allocated for given user
                    if firstRB > lastRB
                        disp([firstRB lastRB])
                    end
                    [firstSubcarrier, lastSubcarrier] = SubcarrierAllocation([firstRB lastRB], ...
                        FFTLength, numSubcarriers, 150);

                    %transmiting ofdm data
                    ofdmData = OneLTEUserTransm(numBits, polar, modType, ...
                        numTSlots, FFTLength, numCPSamples, firstSubcarrier,...
                        ofdmSymbPerTimeSlot);
                    ofdmVectL = numTSlots*ofdmSymbPerTimeSlot*length(ofdmData(:,1));
                    ofdmVect = reshape(ofdmData, [1,ofdmVectL]);
                    LTETransmData = transmissionLTE(ofdmVect, firstTimeSlot, numTSlots, ...
                        LTETransmData, FFTLength, numCPSamples, ofdmSymbPerTimeSlot);
                end

            end
            


            for i = 1:numUsers3
                disp(['Numerology 2 User: ', num2str(i)]);
                while allocPossible3(i)==0
                    numTSlots = randomParam(1, 0, true, 1, timeSlotsNum);                      % number of time slots (time domain)
                    firstTimeSlot = randomParam(5,20, true, 1, timeSlotsNum-numTSlots+1);     % first time slot index
                    numResourceBlocks = randomParam(16, 0, true, 1, 150);               % number of resource blocks (frequency domain)
                    firstRB = randomParam(50+50+25, 8, true, 100, 150-numResourceBlocks+1);    % first resource block index
                    lastTimeSlot = firstTimeSlot + numTSlots - 1;                               % last time slot index
                    lastRB = firstRB + numResourceBlocks - 1;                                    % last resource block index

                    numBits = ofdmSymbPerTimeSlot*numSubcarriers*numBitsPerSymb*...
                        numTSlots*numResourceBlocks;                                            % number of data bits to generate in the beginning of program

                    [resourceAllocShort3, allocPossible3(i)] = CheckResourcesAllocation(firstRB, ...
                        lastRB, firstTimeSlot, lastTimeSlot, resourceAllocShort3, ...
                        i, timeSlotsNum, 150);
                end

                if allocPossible3(i)
                    resourceAlloc(firstRB:lastRB, firstTimeSlot:lastTimeSlot) = ...
                        i*ones(numResourceBlocks, numTSlots);

                    % calculating frequency subcarriers indexes for resource blocks allocated for given user
                    if firstRB > lastRB
                        disp([firstRB lastRB])
                    end
                    [firstSubcarrier, lastSubcarrier] = SubcarrierAllocation([firstRB lastRB], ...
                        FFTLength, numSubcarriers, 150);

                    %transmiting ofdm data
                    ofdmData = OneLTEUserTransm(numBits, polar, modType, ...
                        numTSlots, FFTLength, numCPSamples, firstSubcarrier,...
                        ofdmSymbPerTimeSlot);
                    ofdmVectL = numTSlots*ofdmSymbPerTimeSlot*length(ofdmData(:,1));
                    ofdmVect = reshape(ofdmData, [1,ofdmVectL]);
                    LTETransmData = transmissionLTE(ofdmVect, firstTimeSlot, numTSlots, ...
                        LTETransmData, FFTLength, numCPSamples, ofdmSymbPerTimeSlot);
                end

            end
             
                                
            for i = 1:numUsers4
                disp(['Numerology 1 User: ', num2str(i)]);
                while allocPossible4(i)==0
                    numTSlots = randomParam(4, 0, true, 1, timeSlotsNum);                      % number of time slots (time domain)
                    firstTimeSlot = randomParam(75, 20, true, 70, timeSlotsNum-numTSlots+1);     % first time slot index
                    numResourceBlocks = randomParam(8, 0, true, 1, 100);               % number of resource blocks (frequency domain)
                    firstRB = randomParam(50+50+25, 10, true, 100, 150-numResourceBlocks+1);    % first resource block index
                    lastTimeSlot = firstTimeSlot + numTSlots - 1;                               % last time slot index
                    lastRB = firstRB + numResourceBlocks - 1;                                    % last resource block index

                    numBits = ofdmSymbPerTimeSlot*numSubcarriers*numBitsPerSymb*...
                        numTSlots*numResourceBlocks;                                            % number of data bits to generate in the beginning of program

                    [resourceAllocShort4, allocPossible4(i)] = CheckResourcesAllocation(firstRB, ...
                        lastRB, firstTimeSlot, lastTimeSlot, resourceAllocShort4, ...
                        i, timeSlotsNum, 150);
                end

                if allocPossible4(i)
                    resourceAlloc(firstRB:lastRB, firstTimeSlot:lastTimeSlot) = ...
                        i*ones(numResourceBlocks, numTSlots);

                    % calculating frequency subcarriers indexes for resource blocks allocated for given user
                    if firstRB > lastRB
                        disp([firstRB lastRB])
                    end
                    [firstSubcarrier, lastSubcarrier] = SubcarrierAllocation([firstRB lastRB], ...
                        FFTLength, numSubcarriers, 150);

                    %transmiting ofdm data
                    ofdmData = OneLTEUserTransm(numBits, polar, modType, ...
                        numTSlots, FFTLength, numCPSamples, firstSubcarrier,...
                        ofdmSymbPerTimeSlot);
                    ofdmVectL = numTSlots*ofdmSymbPerTimeSlot*length(ofdmData(:,1));
                    ofdmVect = reshape(ofdmData, [1,ofdmVectL]);
                    LTETransmData = transmissionLTE(ofdmVect, firstTimeSlot, numTSlots, ...
                        LTETransmData, FFTLength, numCPSamples, ofdmSymbPerTimeSlot);
                end

            end
            
        
             for i = 1:numUsers5
                disp(['Numerology 0 User: ', num2str(i)]);
                while allocPossible5(i)==0
                    numTSlots = randomParam(16, 0, true, 1, timeSlotsNum);                      % number of time slots (time domain)
                    firstTimeSlot = randomParam(55, 25, true, 40, timeSlotsNum-numTSlots+1);     % first time slot index
                    numResourceBlocks = randomParam(4, 0, true, 1, 100);               % number of resource blocks (frequency domain)
                    firstRB = randomParam(60+15, 10, true, 60, 100-numResourceBlocks+1);    % first resource block index
                    lastTimeSlot = firstTimeSlot + numTSlots - 1;                               % last time slot index
                    lastRB = firstRB + numResourceBlocks - 1;                                    % last resource block index

                    numBits = ofdmSymbPerTimeSlot*numSubcarriers*numBitsPerSymb*...
                        numTSlots*numResourceBlocks;                                            % number of data bits to generate in the beginning of program

                    [resourceAllocShort5, allocPossible5(i)] = CheckResourcesAllocation(firstRB, ...
                        lastRB, firstTimeSlot, lastTimeSlot, resourceAllocShort5, ...
                        i, timeSlotsNum, 100);
                end

                if allocPossible5(i)
                    resourceAlloc(firstRB:lastRB, firstTimeSlot:lastTimeSlot) = ...
                        i*ones(numResourceBlocks, numTSlots);

                    % calculating frequency subcarriers indexes for resource blocks allocated for given user
                    if firstRB > lastRB
                        disp([firstRB lastRB])
                    end
                    [firstSubcarrier, lastSubcarrier] = SubcarrierAllocation([firstRB lastRB], ...
                        FFTLength, numSubcarriers, 100);

                    %transmiting ofdm data
                    ofdmData = OneLTEUserTransm(numBits, polar, modType, ...
                        numTSlots, FFTLength, numCPSamples, firstSubcarrier,...
                        ofdmSymbPerTimeSlot);
                    ofdmVectL = numTSlots*ofdmSymbPerTimeSlot*length(ofdmData(:,1));
                    ofdmVect = reshape(ofdmData, [1,ofdmVectL]);
                    LTETransmData = transmissionLTE(ofdmVect, firstTimeSlot, numTSlots, ...
                        LTETransmData, FFTLength, numCPSamples, ofdmSymbPerTimeSlot);
                end

            end
            
              for i = 1:numUsers6
                disp(['Numerology 2 User: ', num2str(i)]);
                while allocPossible6(i)==0
                    numTSlots = randomParam(1, 0, true, 1, timeSlotsNum);                      % number of time slots (time domain)
                    firstTimeSlot = randomParam(5, 15, true, 1, timeSlotsNum-numTSlots+1);     % first time slot index
                    numResourceBlocks = randomParam(16, 0, true, 1, 100);               % number of resource blocks (frequency domain)
                    firstRB = randomParam(25, 10, true, 1, 50-numResourceBlocks+1);    % first resource block index
                    lastTimeSlot = firstTimeSlot + numTSlots - 1;                               % last time slot index
                    lastRB = firstRB + numResourceBlocks - 1;                                    % last resource block index

                    numBits = ofdmSymbPerTimeSlot*numSubcarriers*numBitsPerSymb*...
                        numTSlots*numResourceBlocks;                                            % number of data bits to generate in the beginning of program

                    [resourceAllocShort6, allocPossible6(i)] = CheckResourcesAllocation(firstRB, ...
                        lastRB, firstTimeSlot, lastTimeSlot, resourceAllocShort6, ...
                        i, timeSlotsNum, 50);
                end

                if allocPossible6(i)
                    resourceAlloc(firstRB:lastRB, firstTimeSlot:lastTimeSlot) = ...
                        i*ones(numResourceBlocks, numTSlots);

                    % calculating frequency subcarriers indexes for resource blocks allocated for given user
                    if firstRB > lastRB
                        disp([firstRB lastRB])
                    end
                    [firstSubcarrier, lastSubcarrier] = SubcarrierAllocation([firstRB lastRB], ...
                        FFTLength, numSubcarriers, 50);

                    %transmiting ofdm data
                    ofdmData = OneLTEUserTransm(numBits, polar, modType, ...
                        numTSlots, FFTLength, numCPSamples, firstSubcarrier,...
                        ofdmSymbPerTimeSlot);
                    ofdmVectL = numTSlots*ofdmSymbPerTimeSlot*length(ofdmData(:,1));
                    ofdmVect = reshape(ofdmData, [1,ofdmVectL]);
                    LTETransmData = transmissionLTE(ofdmVect, firstTimeSlot, numTSlots, ...
                        LTETransmData, FFTLength, numCPSamples, ofdmSymbPerTimeSlot);
                end

              end
            

%            ALLOCATED RESOURCES VISUALISATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    figure('Renderer', 'painters', 'Position', [200 200 320 200])
                    s = surf(sign(resourceAlloc));
                    view(0, 90)
                    s.EdgeColor = 'none';
                    xlabel('time slot')
                    ylabel('frequency')
                    xlim([1,timeSlotsNum+1])
                    ylim([1,resBlocksNum+1])
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
        
        
            % even out power
            [LTETransmData, powers] = rescaleLTEData(LTETransmData, resourceAlloc,...
                FFTLength, numCPSamples, ofdmSymbPerTimeSlot);
            
            %% AWGN channel

            for SNRind = 1:length(SNR)
                disp(['      SNR = ' num2str(SNR(SNRind))])
                variance = 0.002;           

                %% transmition through channel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                % fading
                LTETransmDataChannel = rchan(LTETransmData.').';
                [C, LTETransmDataChannel] = changeSigPow(variance, ...
                   SNR(SNRind), LTETransmDataChannel);
               
                
                noise =  sqrt(variance/2)*(randn(size(LTETransmDataChannel)) ...
                    + 1i*randn(size(LTETransmDataChannel)));
                noisySignal = LTETransmDataChannel + noise;
                
                %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                [detectedRes, energyMtrx] = energyDetection(FFTLength, ...
                    numSubcarriers, resBlocksNum, ...
                    timeSlotsNum, ofdmSymbPerTimeSlot, numCPSamples, Pf, noisySignal);
                decisionTable = detectedRes;
                [detectProc, falseAlarmProc] = ...
                    detectionCheck(sign(resourceAlloc(1:resBlocksNum, 1:timeSlotsNum)), ...
                    detectedRes);

                detect2(SNRind) = detect2(SNRind) + detectProc;
                falseAl2(SNRind) = falseAl2(SNRind) + falseAlarmProc;


                %%  DETECTED RESOURCES VISUALISATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             figure('Renderer', 'painters', 'Position', [200 200 320 200])
%                             s2 = surf(decisionTable);
%                             s2.EdgeColor = 'none';
%                             view(0, 90)
%                             xlabel('time')
%                             ylabel('frequency')
%                             zlabel('user')
%                             xlim([1,timeSlotsNum+1])
%                             ylim([1,resBlocksNum+1])
%                 %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %                 %% Write to file
                if ifWrite2File == 1
                    if mod(fileIdx,2) == 1
                        fileEnd = '.txt';
                    else
                        fileEnd = '_test.txt';
                    end


 
                    if SNR(SNRind)<0
                        fileName = ['data_' num2str(abs(SNR(SNRind))) fileEnd];
                    else
                        fileName = ['data' num2str(abs(SNR(SNRind))) fileEnd];
                    end
                    fileID = fopen(fileName, 'a+');

                    rB = 1;
                    tS = 1;
                    [maxRB, maxTS] = size(decisionTable);

                    forget_Factor = 0.9;

                    for rB = 1:resBlocksNum

                        serialNeighbours = 0;
                        ff_value = 0;
                        ff_valueEn = 0;

                        for tS = 1:timeSlotsNum
                            if (resourceAlloc(rB,tS) == 1 || resourceAlloc(rB,tS) == 2 || resourceAlloc(rB,tS) == 3 || resourceAlloc(rB,tS) == 4 || resourceAlloc(rB,tS) == 5)
                            


                                fprintf(fileID, '%d %d %d\n', ...
                                    rB, tS, sign(resourceAlloc(rB,tS)));
                            end
                        end
                    end

                    fclose(fileID);
                 end
              end
% %     
% % %% plot energy detecion probability of detection and probability of false alarm
%     figure
%     plot(SNR, detect2, 'b', 'LineWidth', 1)
%     hold on
%     plot(SNR, falseAl2, 'r', 'LineWidth', 1)
%     xlabel('SNR [dB]')
%     ylabel('Probability [%]')
%     legend('probability of detection', 'probability of false alarm');
% % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        end
end