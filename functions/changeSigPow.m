function [C, signalChanged] = changeSigPow(Pn, SNR, signal)
%changeSigPow - changes signal power depending on the noise power and
%signal to noise power ratio.
%   Parameters:
%       - Pn - noise power (sigma^2)
%       - SNR - signal to noise ratio (in dB)
%       - signal - signal samples (vector)
%   Output:
%       - C - rescaling factor
%       - signalChanged - rescaled signal samples

Ps2 = Pn * 10^(SNR/10);
sig = nonzeros(signal)';
Ps1 = sum(abs(sig).^2)/length(sig);

C = sqrt(Ps2/Ps1);

signalChanged = C*signal;
end

