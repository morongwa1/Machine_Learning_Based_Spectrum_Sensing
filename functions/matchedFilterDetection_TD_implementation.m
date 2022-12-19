function [detectedRes, energyMatrix, Pn] = matchedFilterDetection_TD_implementation(FFTLength, ...
    numSubcarriers, resBlocksNum, timeSlotsNum, ofdmSymbPerTimeSlot, numCPSamples, Pf, data, templatedata)
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


detectedRes = zeros(resBlocksNum, timeSlotsNum);
energyMatrix = zeros(resBlocksNum, timeSlotsNum);
testMatrix = zeros(resBlocksNum, timeSlotsNum);

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

%Every time a new 4G grid is incoming having the relevant user, we are
%doing this to extract the data of the relevant user.

for m = 1:timeSlotsNum
    t1 = tInd1;
    t2 = t1 + FFTLength - 1;
    
    Pxx = zeros(size(data(tInd1:tInd2)));
    
    for nS = 1:ofdmSymbPerTimeSlot
        Pxx(nS:ofdmSymbPerTimeSlot:length(Pxx)) = fft(data(t1:t2), FFTLength);
        
        time_domain_pilot_samples = templatedata(t1:t2);
        
        template = fliplr(conj(time_domain_pilot_samples));           %Pilot Signal in time domain ------ conjugated and time reversed

        convolution = conv(template, data(t1:t2));
            
        t1 = t2+1;
        t2 = t1 + FFTLength - 1;
    end
    f_domain_conv = fft(convolution, FFTLength);

    
    z = [Pxx(1:firstFreqIdx-1) Pxx(lastFreqIdx+1:length(Pxx))];     % Data which does not belong to user (margin samples)
    vv1 = sum(abs(z).^2)/(length(z));                               % Noise variance (average value)
    vv = vv1;
    
    noise = imresize(z,size(Pxx));
    threshold = sum(abs(noise.*fft(conj(time_domain_pilot_samples),FFTLength)))/length(noise);
    
    % Threshold:
    Pn = threshold*(qfuncinv(Pf)*sqrt(1/nSampsF) + 1);


    for n = 1:resBlocksNum        
        testMatrix(n,m) = sum(abs(f_domain_conv(fInd1:fInd2)))/nSampsF;
        
        if testMatrix(n,m) > Pn
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