function [detectedRes, energyMatrix, Pn, abs_vv] = energyDetection(FFTLength, ...
    numSubcarriers, resBlocksNum, timeSlotsNum, ofdmSymbPerTimeSlot, numCPSamples, Pf, data, variance)
%ENERGYDETECTION spectrum sensing using energy detection for each
%subcarrier in resource block
%   Parameters:
%       FFTLength - length of OFDM symbol without cyclic prefix
%       numSubcarriers - number of subcarriers in resource block
%       resBlockNum - number of resource blocks in band
%       timeSlotsNum - number of time slots of generated data
%       numCPSamples - number of samples in cyclic prefix
%       Pf - probability of false alarm
%       data - signal for sensing
%       variance - noise variance for checking correctness of noise
%       estimation (not used)
%   Output:
%       detectedRes - detected resources (subcarriers)
%       Pn - detection threshold
%       vv - estimated noise variance

%The grid is made again for every user that is why the indices of the start and end are defined again


L = length(data);

detectedRes = zeros(resBlocksNum, timeSlotsNum);
energyMatrix = zeros(resBlocksNum, timeSlotsNum);

data = cyclicPrefixRemoval(data, FFTLength, numCPSamples);

timeSamps = FFTLength;  

tInd1 = 1;
tInd2 = tInd1 + timeSamps*ofdmSymbPerTimeSlot - 1;    
nSampsT = tInd2 - tInd1 + 1;     

usableSamples = numSubcarriers*resBlocksNum*ofdmSymbPerTimeSlot;   
marginSamples = nSampsT - usableSamples;   

firstFreqIdx = round(0.5*marginSamples+1);    
lastFreqIdx = round(0.5*marginSamples + usableSamples);   

fInd1 = firstFreqIdx;  
fInd2 = fInd1 + numSubcarriers*ofdmSymbPerTimeSlot - 1; 
nSampsF = fInd2 - fInd1+1;   

for m = 1:timeSlotsNum   
    t1 = tInd1;   
    t2 = t1 + FFTLength - 1;  
    
    Pxx = zeros(size(data(tInd1:tInd2)));
    for nS = 1:ofdmSymbPerTimeSlot
        Pxx(nS:ofdmSymbPerTimeSlot:length(Pxx)) = fft(data(t1:t2), FFTLength);
        t1 = t2+1;
        t2 = t1 + FFTLength - 1;
    end
    
    z = [Pxx(1:firstFreqIdx-1) Pxx(lastFreqIdx+1:length(Pxx))];           %collecting the data of margin samples
    vv1 = sum(abs(z).^2)/(length(z));
    vv = vv1;
    abs_vv = abs(vv1-vv)/vv;   %noise power to be used in Q function
    
    % Threshold:
    Pn = vv*(qfuncinv(Pf)*sqrt(1/nSampsF) + 1 );

    
    vect = zeros(1,resBlocksNum);
    
    for n = 1:resBlocksNum   
        vect(n) = sum(abs(Pxx(fInd1:fInd2)).^2)/nSampsF;
        energyMatrix(n,m) = sum(abs(Pxx(fInd1:fInd2)).^2)/nSampsF;

        if energyMatrix(n,m) > Pn
            detectedRes(n,m) = 1;
        else
            detectedRes(n,m) = 0;
        end
        
        if (fInd2 - fInd1 + 1) ~= nSampsF
            disp('ERROR!!!!')
            disp(fInd2 - fInd1 + 1)
            disp(nSampsF)
        end
        
        fInd1 = fInd2 + 1;
        fInd2 = fInd1 + numSubcarriers*ofdmSymbPerTimeSlot - 1;
    end
    
    tInd1 = tInd2 + 1;
    tInd2 = tInd1 + timeSamps*ofdmSymbPerTimeSlot - 1;
    fInd1 = firstFreqIdx;
    fInd2 = fInd1 + numSubcarriers*ofdmSymbPerTimeSlot - 1;
end

end

