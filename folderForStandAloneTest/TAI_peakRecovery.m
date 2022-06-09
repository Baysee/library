% T-TAI on a QAM16 signal
% For a gain of 2 on a 10 GS/s signal. Nb of output peaks: 1 or 2
% need to insure that bit rate and Talbot parameters are commensurate with the length of the whole signal to avoid any sliping.
%% To use discrete spectral phase filtering, uncomment line 114, i.e. [phaseGVD]=genDiscPhase(lent,nus,q,df,discSpec,1);



%% Time Frequency Vectors
lent = 2^22;
time_window = 0.1e-12*lent;%second. Actual time array=+-time_window/2
t = linspace(-time_window/2,time_window/2,lent);
dt = t(2)-t(1); Fs = 1/dt; f = linspace(-Fs/2,Fs/2,lent); df = (f(2)-f(1));
fG=  f*10^-9; tps = t*10^12;



%% Talbot Parameters

q = 10; p = 1;
s = generateSparameter(p,q);



% Adjust the phase levels of the Talbot phase to fit with the bit rate
tqAim=1/(10e9);                      % tq is modified such that ts is commensurate with dt

% No need to modify beloe
tsAim=tqAim/q;
ndt_ts=round(tsAim/dt);
ndt_tq=q*ndt_ts;
ts=dt*ndt_ts;
tq=q*ts;
nus=1/tq;
nuq=nus*q;



% Chose SUT; First a data signal
% From the SUT generation, one needs to define:
% SUT: signal
% nPeaksPerBit can be set to a small integer (usually refers to the number
% of TAI peaks per data signal bit)
% sutUsedRange: The range of interest, to get rid of edge effects
%% Data signal generation


nPeaksPerBit=5;

P = 7;                                                  % PRBS order
L = 2^P-1;                                              % PRBS length

digimod_format ='SP-PAM'; 'QAM';                                 % Digital modulation format
M = 2;%16;                                                 % Number of symbols (modulation order)
ndt_bitRate = (ndt_ts*q*nPeaksPerBit);


K = log2(M);                                            % Number of bits per symbol
N = L*K;                                                % Number of bits
% B = idinput([N,2],'prbs',[1,1],[0,1]).'; B = B(1,:);    % Bit sequence
B = randi(2,[1,N])-1;
D = bi2de(reshape(B,L,K)).';                            % Symbol indices


% Symbol sequence % QAM
S = (qammod(D,M)+1)/2;


S_I = real(S);                                              % In-phase channel
S_Q = imag(S);                                              % Quadrature channel
S_M = abs(S);                                               % Magnitude
S_P = angle(S);                                             % Phase

% figure % Constellation
% plot(S,'o')
% axis square

SUT = repelem( S, ndt_bitRate);                              % repeat each element ndt_bitRate time to get to sampling rate
SUT = repmat (SUT, 1, ceil(lent/(ndt_bitRate*L)) );
SUTmids=round(ndt_bitRate/2):round(ndt_bitRate):lent; 
range=-round(ndt_bitRate*0.6/2):round(ndt_bitRate*0.6/2);
sutUsedRange=reshape(SUTmids+range',[numel(range)*numel(SUTmids),1]);

%% Alternative signal: a sinusoid
% 
% 
% nPeaksPerBit=50; % Number of TAI peaks per sinusoid period
% singleSineLen=ndt_tq*nPeaksPerBit; % Length of a single sinusoid period
% singleSine=sin(2*pi/(ndt_tq*nPeaksPerBit)*[0:singleSineLen-1]); % Single sinusoid period
% SUT = repmat (singleSine, 1, ceil(lent/(singleSineLen)) ); %% Repeat single period many times; Need to make SUT a bit too long; will use only what is necessary later
% sutUsedRange=round(lent/20):round(lent-lent/20); % Region of analysis
% 
% 



%% Keep section below


 % Redefine time=frequency vectors for new length to ensure lent is commensurate with bit length
bitRate = 1/(ndt_bitRate*dt);
 lent=ndt_bitRate*ceil(lent/ndt_bitRate); 
t = (1:lent)*dt-round(lent)*dt;
f = linspace(-Fs/2,Fs/2,lent); df = (f(2)-f(1));
fG=  f*10^-9; tps = t*10^12;
CW=ones(1,lent);
SUT = SUT(1:lent);
% SUT=CW;
SUT=SUT(1:lent);
SUT_f = nfft(SUT,dt);

% Filter SUT
SUTraw=SUT;
filtFreqSUT=(bitRate)*40; plotit=0;
[SUT]=(filterSig(1,round(1/df*filtFreqSUT),SUTraw,plotit));
SUT_f = nfft(SUT,dt);

figure; plot(abs(SUT)); 
hold on ; plot(sutUsedRange,abs(SUT(sutUsedRange)),'.')

%% Talbot Phase Manipulations
%%%%%%%%%%%%%%%%%
%%%% GVD %%%%%%%%
%%%%%%%%%%%%%%%%%


phi2 = p*q*ts^2/(2*pi);                                 % GVD From Theoretical parameters (using s parameter)
phi2perKm =   2.1823e-23;  % GVD per KM
nkm = phi2/phi2perKm;

phaseGVD = phi2/2*(2*pi*f).^2;                          % Phase Equation


% Discrete spectral phase filtering filtering
discSpec=wrapTo2Pi(s/q*pi*(0:q-1).^2);   % Values of Discrete filter
% [phaseGVD]=genDiscPhase(lent,nus,q,df,discSpec,1);
% phaseGVD=circshift(phaseGVDini,(fnus-vChange));   % Circshift the



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Temporal Phase Modulation %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

GV=wrapTo2Pi(pi*p/q*(1+q*mod(q,2))*(1:q).^2); % Use this form with a Discrete filter
GV=circshift(GV,floor(q/2)+1);

%% Corrections for AWG

[phaseTemporalRaw]=genDiscPhase(lent,ts,q,dt,GV,0);

filtFreq=(6/ts); plotit=1;                            % Filter Temporal phase mod signal to simulate finite AWG RF-bandwidth
phaseTemporal=phaseTemporalRaw;                         % Record signal without filter if needed
% Uncomment Line below to filter temporal phase signal
   [phaseTemporal]=real(filterSig(3,round(1/df*filtFreq),phaseTemporalRaw,plotit));                    


phaseTemporal_f=nfft(phaseTemporal,dt);




%% Talbot Processing

% Temporal Phase Modulation

PM = SUT.*exp(1j*phaseTemporal);
PM_f = nfft(PM,dt);

% Dispersive Propagation

disp_f = PM_f.*exp(1j*phaseGVD);
disp = nifft( disp_f,Fs);

% Plot basic results

dispPow=abs(disp).^2;

figure;
subplot(4,1,1)
plot(f,abs(disp_f).^2);
yyaxis right; 
plot(f,phaseGVD)
subplot(4,1,2)
plot(f,abs(disp_f).^2);
% yyaxis right; 
% plot(f,phaseGVD)
xlim([-2/tq 2/tq])

% temporal Plots
subplot(4,1,3)
plot(t,abs(disp).^2);
yyaxis right
plot(t,phaseTemporal)
subplot(4,1,4)
plot(t,abs(disp).^2);
yyaxis right
plot(t,phaseTemporal)
xlim([-2*tq 2*tq])




%%%%%%
% Find Peaks using a constant val signal
%%%%%%

testpk=ones(1,lent);

PM_pk = testpk.*exp(1j*phaseTemporal);                              % Talbot on CW
PM_f_pk = nfft(PM_pk,dt);
disp_f_pk = PM_f_pk.*exp(1j*phaseGVD);
disp_pk = nifft( disp_f_pk,Fs);

cwpk_pow=abs((disp_pk).^2);

searchRange=1:(ndt_tq*3);
cwRange=cwpk_pow(searchRange);
[maxVal,maxInd]=max(cwRange);
rise=find(cwRange(1:maxInd)<maxVal/2,1,'last');
fall=maxInd-1+find(cwRange(maxInd:end)<maxVal/2,1);
mid=round((rise+fall)/2);

peakInds=(mid+4*ndt_tq):ndt_tq:(lent-4*ndt_tq);

omitPeaks=logical(ismember(peakInds,sutUsedRange)-1);
peakInds(omitPeaks)=[];





% Take Averages near peak val
nAvs = 2 ;                                                          % Make even number. will check this number of indices each side and take average
avInds = repelem(peakInds,2*nAvs+1,1); pms=-nAvs:nAvs;              % Index over which to average
avIndsAll_Mat = avInds + pms';                                      % Matrix form
avIndsAll=reshape(avIndsAll_Mat,1,numel(avIndsAll_Mat));            % Vector form

avIndsAll_Mat(avIndsAll_Mat<1)=1;
avIndsAll_Mat(avIndsAll_Mat>lent)=lent;
nearPeakVals_Mat = disp(avIndsAll_Mat);                             % Values to average in Mat form
nearPeakVals=reshape(nearPeakVals_Mat,1,numel(nearPeakVals_Mat));   % Values to av in Vector form
meanVals=mean(nearPeakVals_Mat);                                    % Mean


% Find constant phase difference
sampSUT_P=wrapTo2Pi(angle(SUT(avIndsAll)));
sampOut_Pini=wrapTo2Pi(angle(disp(avIndsAll)));
angleDiff=(sampSUT_P-sampOut_Pini);
[negInds]=find(angleDiff<0);                        % Find negative values, then remove 2pi

sampOut_Pini(negInds)=sampOut_Pini(negInds)-2*pi;
angleDiff=(sampSUT_P-sampOut_Pini);

mean_AD=mode(angleDiff);                        % angle difference that appears most often
% SampSUT_P=wrapTo2Pi(sampSUT_Pini);
sampOut_P=(sampOut_Pini+mean_AD);
[negInds]=find(sampOut_P<0); 
sampOut_P(negInds)=sampOut_P(negInds)+2*pi;
% ploting section


% Parameters
ht = t(round(lent/2));
tlim = ([ht ht+20/bitRate])*10^12;
tmid=mean(tlim); tps=tps-tmid; tlim=tlim-tmid;






%% Plotting section

% FOnt Sizes 
AFS=14;

figSize=[0, 0, 14, 8];
fg=figure('rend','painters', 'Units', 'Inches', 'Position', figSize);

subplot(3,5,[3:5, 8:10])
ax=gca;
ax.Position=[0.4217 0.3553 0.4833 0.5697];
maxSUT=max(abs(SUT).^2);

h_TAI=plot(tps,abs(disp).^2/maxSUT,'color',[0.8500    0.3250    0.0980]);
hold on
avIndsAll
plot(tps(peakInds),abs(disp(peakInds)).^2/maxSUT,'o','color',[0.8500    0.3250    0.0980])
plot(tps(avIndsAll),abs(disp(avIndsAll)).^2/maxSUT,'x','color',[0.8500    0.3250    0.0980])
h_SUT=plot(tps,abs(SUT).^2/maxSUT,'color',[ 0    0.4470    0.7410]);
plot(tps(peakInds),abs(SUT(peakInds)).^2/maxSUT,'o','color',[ 0    0.4470    0.7410])
ylabel('Power (n.u.)')

% yyaxis right
% h_SUT_P=plot(tps(avIndsAll),sampSUT_P,'--','color',  [0.4660    0.6740    0.1880]);
% hold on
% h_TAI_P=plot(tps(avIndsAll),sampOut_P,':','color',  [0.4660    0.6740    0.1880]);
% plot(tps,angle(disp))


xlim(tlim); ylim([0 21])
xlabel('Time (ps)')
ylabel('Phase (rad)')
set(gca,'FontSize',AFS);
% legend([h_SUT h_TAI h_SUT_P h_TAI_P],'SUT power','TAI power','SUT phase','TAI phase','Location','northwest')

% subplot(3,5,[13])
% ax=gca;
%  ax.Position=[0.4206 0.1100 0.1587 0.1747];
% plot(tps,abs(disp_pk).^2,'color',[0.4246 0.1100 0.1548 0.1660]);
% hold on
% plot(tps(peakInds),abs(disp_pk(peakInds)).^2,'o','color',[0.8500    0.3250    0.0980])
% plot(tps(avIndsAll),abs(disp_pk(avIndsAll)).^2,'.','color',[0.8500    0.3250    0.0980])
% xlim([-tq tq]*10^12);
% xlabel('Time (ps)')

subplot(3,5,[13:15])
ax=gca;
ax.Position=[0.4216 0.1100 0.4834 0.1730];%[0.6184 0.1100 0.2866 0.1730];
h_TAI=plot(tps,abs(disp).^2/maxSUT,'color',[0.8500    0.3250    0.0980]);
hold on
plot(tps(peakInds),abs(disp(peakInds)).^2/maxSUT,'o','color',[0.8500    0.3250    0.0980])
plot(tps(avIndsAll),abs(disp(avIndsAll)).^2/maxSUT,'.','color',[0.8500    0.3250    0.0980])
% xlim([-2*tq 2*tq]*10^12+2*tq*10^12);
xlim([-2*tq 2*tq]*10^12);


% Constellation
subplot(3,5,[1,2,6,7])
ax=gca;
ax.Position=[0.0506 0.3559 0.3214 0.5691];

% plot(abs(disp(avIndsAll)).*exp(1j*sampOut_P),'o','color',[0.9290    0.6940    0.1250])
% hold on

%% Unprocessed constallation up to a phase shift
if mod(q,2)==0
plot(disp(avIndsAll)*exp(1j*pi/4),'o','color',[0.8500    0.3250    0.0980])
else
 plot(disp(avIndsAll),'o','color',[0.8500    0.3250    0.0980])
end

%% Processed constellation
% if mod(q,2)==0
% plot(disp(avIndsAll)/sqrt(q)*exp(1j*pi/4)-0.065+1j*0.032 -mean(disp(avIndsAll)*exp(1j*pi/4)/sqrt(q))+mean(SUT(avIndsAll)),'o','color',[0.8500    0.3250    0.0980])
% else
%  plot(disp(avIndsAll)/sqrt(q) -mean(disp(avIndsAll)/sqrt(q))+mean(SUT(avIndsAll)),'o','color',[0.8500    0.3250    0.0980])
% end
hold on
plot((SUT(avIndsAll)),'o','color',[ 0    0.4470    0.7410]);
% plot(SUTraw(avIndsAll),'o');
xlabel('In phase')
ylabel('Quadrature')


% Details
subplot(3,5,[11 12])
ax=gca;
dim=[0.0506 0.1100 0.3274 0.1574];
ax.Position=dim;
str=['Bit Rate: ', num2str(bitRate*10^-9,3), ' GHz', newline, ...
    num2str(nPeaksPerBit) ' Peak/bit, q=' num2str(q), newline, ...
    'Peak width: ', num2str(ts*10^12), ' ps', newline,...
    'AWG sampling Rate: ', num2str(1/ts*10^-9,4), ' GHz', newline,...
    'Peak separation: ' , num2str(tq*10^12,4), ' ps', newline,...
    'Dispersion: ', num2str(phi2*10^24,4), ' ps^2 (', num2str(nkm,3), ' km SMF)']
annotation('textbox',dim,'String',str)
xticks(''); yticks('')
set(findall(gcf,'-property','FontSize'),'FontSize',14)



function s=generateSparameter(p,q)

if mod(q,2)==1
    s=2*mulinv(2*p,q);
else
    s= mulinv(p,2*q) ;
end

end


function y=mulinv(x,p)
% 
% if ~isprime(p)
%     disp('The field order is not a prime number');
%     return
% elseif sum(x>=p)
%     disp('All or some of the numbers do not belong to the field');
%     return
% elseif sum(~x)
%     disp('0 does not have a multiplicative inverse');
%     return
% end

k=zeros(size(x));   %set the counter array
m=mod(k*p+1,x);     %find reminders
while sum(m)        %are all reminders zero?
    k=k+sign(m);    %update the counter array
    m=mod(k*p+1,x); %caculate new reminders 
end
y=(k*p+1)./x;       %find the multiplicative inverses of X
end

function ifftout=nifft(sig,varargin)
%%%Need to first ifftshift to swap for ifft (Matlab assumes swaped fft)
%normalize if varargin is given

unnormdifft=fftshift(ifft(ifftshift(sig)));

scale=nargin-1;

if scale==1
    Fs=varargin{1};
ifftout=Fs*unnormdifft;

else
ifftout=1/(max(abs(unnormdifft)))*unnormdifft;

end

end

function fftout=nfft(sig,varargin)
%fft gives the swaped spectrum.
unnormdfft=fftshift(fft(ifftshift(sig)));

scale=nargin-1;

if scale==1
    dt=varargin{1};
fftout=dt*unnormdfft;
        else
        fftout=1/(max(abs(unnormdfft)))*unnormdfft;
end
end

function [output]=filterSig(m,t0,sig,plotit)

inds=1:numel(sig); half=round(numel(sig)/2);

Sig_f=nrmd_fft(sig,1,1);
filt=superGauss(0,t0,m,inds,half);

sigSize=size(sig);
filtSize=size(filt);

if sigSize(1)~=filtSize
    filt=filt'
end
output=nrmd_ifft(Sig_f.*filt,1,1);
if plotit
    
    figure;
    subplot(2,1,1)
    plot(abs(Sig_f));
    yyaxis right
    plot(filt)
    subplot(2,1,2)
    plot(abs(sig))
    hold on
    plot(abs(output))
    plot(real(output))
    legend('input','abs','real')
end


end



function waveform=superGauss(C,t0,m,xs,center)
waveform=exp(-(1+1j*C)/2*((xs-center)/t0).^(2*m));
end

function [phase]=genDiscPhase(lent,level,m,dt,phaseVals,offset)
%level is the time/freq length of each phase level. Offset sets the phase to have the middle of a phaselevel at the center of lent
% 
% 
% if length(phaseVals)~=m
%    'error! not enough phase vals!'
% end

m=numel(phaseVals);

npLevel=round(level/(dt));%Number of points per unit. Unit~Tr contains m pulses, contained within period T1*m. level_t is Tr/m
npUnit=npLevel*m;

nReps=ceil(lent/npUnit);

% unit=zeros(npLevel,m);%Square matrix m by npLevel_t, with all phase values
% for i=1:numel(unit(1,:))
%     unit(:,i)=phaseVals(i);%GaussVals created at begining, from s & m.
%     
% end

unit=repelem(phaseVals,npLevel); %repeat each phase val 



phaseTemporal=repmat(unit,1,nReps); 
phase= phaseTemporal(1:lent);
if offset==1
    
   cent=round(lent/2);
   left=mod(cent,npLevel);
   
   shift=left-round(npLevel/2);
   phase = circshift(phase,shift);
   
end
end
