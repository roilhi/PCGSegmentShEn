function [ShEn] = ShannonSmoothEnvelope(x,th,Fs)
%===================================================================
%  ____  _
% / ___|| |__   __ _ _ __  _ __   ___  _ __  ___
% \___ \| '_ \ / _` | '_ \| '_ \ / _ \| '_ \/ __|
%  ___) | | | | (_| | | | | | | | (_) | | | \__ \
% |____/|_| |_|\__,_|_| |_|_| |_|\___/|_| |_|___/
% 
%  _____                _
% | ____|_ ____   _____| | ___  _ __   ___
% |  _| | '_ \ \ / / _ \ |/ _ \| '_ \ / _ \
% | |___| | | \ V /  __/ | (_) | |_) |  __/
% |_____|_| |_|\_/ \___|_|\___/| .__/ \___|
%                              |_|
%   ____      _            _       _   _
%  / ___|__ _| | ___ _   _| | __ _| |_(_) ___  _ __
% | |   / _` | |/ __| | | | |/ _` | __| |/ _ \| '_ \
% | |__| (_| | | (__| |_| | | (_| | |_| | (_) | | | |
%  \____\__,_|_|\___|\__,_|_|\__,_|\__|_|\___/|_| |_|
%
% Function based on [1]:
% [1] Varghees, V. N., & Ramachandran, K. I. (2017). 
%     Effective Heart Sound Segmentation and Murmur Classification Using 
%     Empirical Wavelet Transform and Instantaneous Phase for 
%     Electronic Stethoscope. IEEE Sensors Journal, 17(12), 3861-3872.
% Created by: Roilhi Frajo Ibarra Hernandez, May 2018
% frajo@cicese.edu.mx; roilhi-frajo.ibarra-hernandez@irisa.fr

% 1) Adaptive Amplitude Tresholding
N = length(x);
L = 0.03*Fs; % Window length of 30ms
DL = buffer(x,L); %Signal windowing in non-overlaping frames (columns)
muL = mean(DL,2);
% Adaptive noise-level tresholding
eta_n = median(muL)/0.6745;
eta_n = eta_n+0.01;
xTh = zeros(size(x));
for k=1:N
    if abs(x(k))>eta_n
        xTh(k) = abs(x(k));
    else
        xTh(k) = 0;
    end
end
% Shannons entropy calculation
ShPrEn = -xTh.*log(xTh);

ShPrEn(isnan(ShPrEn))=0;

% Zero-phase smoothing filter
% Creating the impulse response (rectangular pulse)
h_k = [zeros(1,0.090*Fs) ones(1,0.060*Fs)];
% Filtering to smooth the entropy and getting the envelope
ShEn = filtfilt(h_k,1,ShPrEn);
ShEn = ShEn/max(ShEn);
%sigma = std(Sm);
sigma = th;
% Envelope tresholding
for r = 1:length(ShEn)
    if ShEn(r) > sigma
        ShEn(r)=ShEn(r);
    else
        ShEn(r)=0;
    end
end



