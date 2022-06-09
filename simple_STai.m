%Perform spectral TAI.

%% ---Define time-frequency vectors

lent=2^20;
time_window=20e-9;%second. Actual time array=+-time_window/2

t=linspace(-time_window/2,time_window/2,lent);
dt=t(2)-t(1);Fs=1/dt;
f=linspace(-Fs/2,Fs/2,lent);
df=(f(2)-f(1));
fG=f*10^-9;tps=t*10^12;%GHz
scale=1;


%% -- Talbot Conditions; GVD and Temporal Phase set -------------- %%%%%%%%

% p spectRal phase Parameter. s Temporap phase parameter
% m multiplicationparameter
q=20;p=1;
s=generateSparameter(p,q);%mod(1+m*mod(m,2),2*m);


%-- Define GVD via sampling parameters for temporal phase modulation -----%
% % 
% nuR=12*18e9;
% tR=1/(nuR);%Single temporal level duration



%% S-TAI
% 
% nuR=1.0e9;
% nus=2e9;
nus=2.98e9;
nuq=nus*q;

ts=1/nuq;
tq=q*ts;

% %% define sig
 BW=400e9;

SUT_f=superGauss(0,BW,4,f,-200e9)+superGauss(0,BW,4,f,200e9);
SUT_f=superGauss(0,BW,4,f,0);
SUT=nifft(SUT_f,Fs);
% 
% bitPeriod=1/400e9;
% nBitPeriod=round(bitPeriod/dt);
% SUT=digiModFunc('SP-PAM',2,16,nBitPeriod,0);
% 
% SUT=SUT(1:lent).*superGauss(0,1/8*time_window,20,t,0);
% SUT_f=nfft(SUT,Fs);
% 
% 

x=t; Tr=1/nus;%0.1e-9;
fwhm=dt*15; center=0;
% SUT=generatePeriodicGauss(x,Tr,fwhm,center).*superGauss(0,30*Tr,8,t,0);
% SUT_f=nfft(SUT,dt);
% % figure;plot(t,abs(in));

%% --  Phase Modulation ---------------------------------- %%%%%%%%

%% SpectRal phase filtering
% 
% % %%% Discrete Spectral Phase
% discSpec=wrapTo2Pi(s/q*pi*(0:q-1).^2); 
% discSpec=wrapTo2Pi(p(1+q*mod(q,2))/q*pi*(0:q-1).^2); 
% [phaseGVD]=genDiscPhase(lent,nus,q,df,discSpec,1);
% 

%%%% Continuous spectral phase
phi2=s/q*2*pi/(2*pi*nus)^2;%1/(m*2*pi*nus^2);%p/m*2*pi/(2*pi*nus)^2;
phi2perKm=   2.1823e-23;
phi2/phi2perKm
phaseGVD=phi2/2*(2*pi*f).^2;

comb=generatePeriodicGauss(f,1/Tr,df,center)/(max(generatePeriodicGauss(f,1/Tr,df,center)));
% figure;plot(comb); hold on; plot(phi2/2*(2*pi*f).^2)

figure;plot(fG,wrapTo2Pi(comb.*(phi2/2*(2*pi*f).^2)))
hold on
plot(fG,circshift(phaseGVD,round(nus/df)))
% plot(abs((SUT_f)/max(abs(SUT_f))))


%% Temporal phase modulation
%S-TAI

%%% Discrete 
GV=wrapTo2Pi(p*(1+q*mod(q,2))/q*pi*(0:q-1).^2);
% GV=wrapTo2Pi(-pi*(m-1)/m*(0:m-1).^2);
[phaseTemporalRaw]=genDiscPhase(lent,ts,q,dt,GV,1);
%  phaseTemporal=phaseTemporalRaw;
 % Filter signal to a 
 phaseTemporal=real(filtSG_tf(phaseTemporalRaw,t,f,round(60e9/df),1,1));

%%% Can also generate continuous temporal phase, but requires very big
%%% excursion for interesting results

figure;plot(t*1e12,phaseTemporal)
xlim([-1.5*tq*1e12 1.5*tq*1e12])
yticks([0:0.25:2]*pi)


%% Talbot Magic - Dispersion then temporal Phase modulation (S-TAI) - %%%%%
%%S-TAI 

dispersed_f=SUT_f.*exp(1j*phaseGVD);
dispersed=nrmd_ifft(dispersed_f,Fs,scale);

temporalRaw=dispersed.*exp(1j*phaseTemporal);
spectrumRaw=nrmd_fft(temporalRaw,dt,scale).*exp(-1j*phaseGVD);
temporalRaw=nifft(spectrumRaw,Fs);

[pks,locs]=findpeaks(abs(spectrumRaw).^2,'MinPeakDistance',round(nuq/df*0.9));


%% PlottingSection --------------------------------------------- %%%%%%%%%%
FS=18;


xlimt=['auto'];
xlimf=['auto'];
figure
subplot(2,1,2)


nrmFac=max(abs(SUT).^2);
plot(tps,abs(SUT).^2/nrmFac,'DisplayName','Input')
hold on
ylabel('Intensity')
plot(tps,(abs(temporalRaw).^2)/nrmFac,'DisplayName','Temporal S-TAI')
hold on
% plot(tps,abs(temporal_OFT_pmOff).^2/nrmFac,'-b','DisplayName','OFT (pm off)')
% hold on
% plot(tps,abs(temporal_OFT_pmOn).^2/nrmFac,'-r','DisplayName',['OFT (pm on, amp: ' num2str(OFTamp)])
% 
% ylabel('Intensity')
% yyaxis right
% 
% plot(tps,phaseTemporal,'DisplayName','Temporal Phase')
xlabel('time (ps)')
xlim(xlimt)
set(gca,'FontSize',FS)
legend('Show')
% -- SpectRal Plot --------------------------------%

subplot(2,1,1)
hold on
plot(fG,abs(SUT_f).^2,'DisplayName','Input')
% plot(fG,10*log10(abs(spectrumRaw).^2),'DisplayName','Output')
% yyaxis right
hold on
plot(fG,(abs(spectrumRaw).^2),'DisplayName','Output')
xlabel('frequency (GHz)')
ylabel('Intensity')
yyaxis right
ca=clean_angle(spectrumRaw);
plot(fG,phaseGVD,'DisplayName','GVD phase')
% plot(fG,wrapTo2Pi(ca))
% alpha(0.3)
% ylabel('Phase (Rad)')

xlim(xlimf)
set(gca,'FontSize',FS)
legend('Show')





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



function [sigFilt,sig_f,wind,sigFilt_f]=filtSG_tf(sig,t,f,nBW,m,displayPlot)
%,varargin)
% Varargin should

lent=numel(sig);
inds=1:lent;


sig_f=nfft(sig,1);

wind=zeros(1,lent);
wind=superGauss(0,nBW,m,inds,round(lent/2));

sigFilt_f=sig_f.*wind;
sigFilt=nifft(sigFilt_f,1);

if displayPlot==1
    figure;
    subplot(2,1,2)
    plot(f,real(sig_f));
    hold on
    plot(f,real(sigFilt_f));
    yyaxis right
    plot(f,wind);
    subplot(2,1,1)
    plot(t,real(sig))
    hold on
    plot(t,abs(sig))
    plot(t,real(sigFilt))
    plot(t,abs(sigFilt))
end

end


% 
function s=generateSparameter(p,q)

if mod(q,2)==1
    s=2*mulinv(2*p,q);
else
    s= mulinv(p,2*q) ;
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
