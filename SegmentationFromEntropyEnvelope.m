% =========================================================================
%
%     ____   ____ ____
%    |  _ \ / ___/ ___|
%    | |_) | |  | |  _
%    |  __/| |__| |_| |
%    |_|    \____\____|
%    
%     ____                                  _        _   _
%    / ___|  ___  __ _ _ __ ___   ___ _ __ | |_ __ _| |_(_) ___  _ __
%    \___ \ / _ \/ _` | '_ ` _ \ / _ \ '_ \| __/ _` | __| |/ _ \| '_ \
%     ___) |  __/ (_| | | | | | |  __/ | | | || (_| | |_| | (_) | | | |
%    |____/ \___|\__, |_| |_| |_|\___|_| |_|\__\__,_|\__|_|\___/|_| |_|
%                |___/
%      __                       _____       _
%     / _|_ __ ___  _ __ ___   | ____|_ __ | |_ _ __ ___  _ __  _   _
%    | |_| '__/ _ \| '_ ` _ \  |  _| | '_ \| __| '__/ _ \| '_ \| | | |
%    |  _| | | (_) | | | | | | | |___| | | | |_| | | (_) | |_) | |_| |
%    |_| |_|  \___/|_| |_| |_| |_____|_| |_|\__|_|  \___/| .__/ \__, |
%                                                        |_|    |___/
%     _____                _
%    | ____|_ ____   _____| | ___  _ __   ___
%    |  _| | '_ \ \ / / _ \ |/ _ \| '_ \ / _ \
%    | |___| | | \ V /  __/ | (_) | |_) |  __/
%    |_____|_| |_|\_/ \___|_|\___/| .__/ \___|
%                                 |_|
% ========================================================================
clear
close all
clc

addpath('/usr/local/mptk/matlab/'); % Path to the MP decomposition

% Dictionary for the MP decomposition
dict = dictread('/usr/local/mptk/reference/dictionary/gabor_32_64_FFT1024.xml');


Fs = 2000;
% -------------------------------------------------
[PCG,segCell] = readPCGsignal('e',210); %Replace this line with the audioread function of MATLAB
% PCG = audioread(Path, bla, bla bla)

%[despiked_signal] = schmidt_spike_removal(PCG, Fs);
% Frequencies for the signal filtering (low and high) in Hz
wn1 = 25/(Fs/2);
wn2 = 600/(Fs/2);
% 6th order Butterworth filtering
[b,a] = butter(6,[wn1 wn2]);
% Signal filtering
PCG = filtfilt(b,a,PCG);
% Signal normalization
PCG = (PCG-mean(PCG))./std(PCG);
PCG = PCG./max(PCG);

N = length(PCG);
NFr = Fs/2; % 0.5 seconds frames
numIter = round(N/NFr)*15; % 15 atoms/0.5 seconds
% MP decomposition of the signal (mptk function)
[book,residual,~] = mpdecomp(PCG,Fs,dict,numIter);
RecSig = mprecons(book,dict);

% Function to detect the Shannon Envelope, we can play with the amplitude 
% threshold (0.5 in this case).
ShEn = ShannonSmoothEnvelope(RecSig,0.5,Fs);

% Derivative of the entropy function to detect the pre-onsets and
% pre-offsets 
Time = (1:N)/Fs;
DerSh = diff(ShEn)'./diff(Time);
DerSh = [0 DerSh]; %Zero padding
% Normalizing the derivative
DerSh = DerSh/max(DerSh);
% Thresholding of the derivative function
for d=1:length(DerSh)
    if DerSh(d)<=0.1 && DerSh(d)>=-0.1
        DerSh(d)=0;
    end
end
DerSh= round(DerSh);

% Position of detected onsets (pre-onsets)
PosOnsets = find(DerSh==1);
%Position of detected offsets (pre-offsets)
PosOffsets = find(DerSh==-1);
% Checking out the distance beweeen onsets and offsets (in miliseconds)
DisOnOff = 1000*((PosOffsets-PosOnsets)/Fs);

% Calculating the difference between possible offsets (in miliseconds) 
OnDif = 1000*(diff(PosOnsets)/Fs);
OffDif = 1000*(diff(PosOffsets)/Fs);
% Checking out wich onsets are too close (<210ms) 1st criteria
CloseOnsets = find(OnDif<=210);
% Checking the offsets wich are too close (<200ms) 2nd criteria
CloseOffsets = find(OffDif<=200);

% Checking the onset-offset intervals of less than 70ms. 3rd criteria
ThCrit = find(DisOnOff<=60);
DisOnOffMax = find(DisOnOff>=150);

% Checking the onsets-offset distance which is too large (more than 750 ms). 
FarOnsets = find(OnDif>=750);

PointerOn = zeros(size(DerSh));
PointerOff = PointerOn;
PointerTh = PointerOff;

IndexOn = PosOnsets(CloseOnsets);
IndexOff = PosOffsets(CloseOffsets+1);
Index3 = PosOffsets(ThCrit);
IndexMaxLen = PosOffsets(DisOnOffMax);
InFarOn = PosOnsets(FarOnsets);


PointerOn(IndexOn) = 1;
PointerOff(IndexOff) = -1;

% The last criteria is to check the amplitude levels of segments detected
% by the other 3 criterias

% a) Detecting the union of two close onsets and offsets
OnUOff = ismember(CloseOnsets,CloseOffsets);
%OnUOff = [CloseOnsets CloseOffsets];
% Location of all onsets and onsets detected too close
OnUOff = CloseOnsets(OnUOff==1);
% Location of all onsets/offsets too close, and too short

OnOffSh = ismember(OnUOff,ThCrit);
OnOffSh = OnUOff(OnOffSh==1);
%AllCrit = sort([OnOffSh OnUOff ThCrit]);
AllCrit = sort([OnOffSh ThCrit]);
CheckAmpOn = PosOnsets(AllCrit);
CheckAmpOff = PosOffsets(AllCrit);
CheckedAmp = zeros(1,length(CheckAmpOn));
for k=1:length(CheckAmpOn)
    x1 = CheckAmpOn(k);
    x2 = CheckAmpOff(k);
    AmpInt = PCG(x1:x2);
    if max(abs(AmpInt))< 0.4
        CheckedAmp(k) =1;
    end
end
% Indexes which have the 4 criteria (offsets to be eliminated)
In4crit = PosOnsets(AllCrit(CheckedAmp==1));

In4crit2 = [In4crit PosOffsets(AllCrit(CheckedAmp==1))];
NewSeg = DerSh;

NewSeg(In4crit2) = 0;
HMM = [segCell{:,1}];
SEG = zeros(1,N);

figure
subplot(211)
plot(Time,PCG), hold on, plot(Time,DerSh),grid, plot(Time,ShEn), plot(Time(IndexOn),PointerOn(IndexOn),'c*'), 
plot(Time(IndexOff),PointerOff(IndexOff),'k*'), plot(Time(Index3),PointerTh(Index3),'r*'),
plot(Time(In4crit),PointerTh(In4crit),'g*'), plot(Time(InFarOn),DerSh(InFarOn),'m*'),
plot(Time(IndexMaxLen),PointerTh(IndexMaxLen),'b*')
subplot(212)
plot(Time,PCG),hold on, plot(Time,NewSeg),plot(Time(HMM),SEG(HMM),'r*'),grid on,
plot(Time(IndexOn),PointerOn(IndexOn),'c*'), plot(Time(IndexOff),PointerOff(IndexOff),'k*')
plot(Time(InFarOn),DerSh(InFarOn),'m*'),
plot(Time(IndexMaxLen),PointerTh(IndexMaxLen),'b*')
figure
subplot(2,1,1)
plot(Time,PCG),hold on, plot(Time,NewSeg),plot(Time(HMM),SEG(HMM),'r*'),grid
 subplot(2,1,2)
 bookplot(book),grid, hold on, plot(Time,NewSeg*200,'r'),plot(Time(HMM),SEG(HMM),'b*')
%  
%  figure, plot(atoms(5,:)./Fs,atoms(2,:),'bx'),grid
% hold on, plot(Time,NewSeg*4,'r'), hold on, plot(Time,PCG),plot(Time(HMM),SEG(HMM),'k*'),plot(Time,ShEn), plot(Time(IndexOn),PointerOn(IndexOn),'c*'), 
% plot(Time(IndexOff),PointerOff(IndexOff),'k*'), plot(Time(Index3),PointerTh(Index3),'r*'),plot(Time(In4crit),PointerTh(In4crit),'g*')
 
 % subplot(2,1,2)
% plot(atoms(5,:)./Fs,atoms(3,:).*Fs,'ko'), hold on, plot(Time,abs(NewSeg)*200,'r'),grid
% figure
% plot(PCG), hold on, plot(DerSh),grid, plot(ShEn), plot(abs(DerSh))

