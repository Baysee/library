%Perform spectRal TAI. SpectRal phase mod favours continuous ie quadratic
%phase variation. Temporal phase modulation can be Discrete, Quadratic or
%sinusoidal. 

%% ---Define time-frequency vectors

lent=2^20;
time_window=50e-9;%second. Actual time array=+-time_window/2

t=linspace(-time_window/2,time_window/2,lent);
dt=t(2)-t(1);Fs=1/dt;
f=linspace(-Fs/2,Fs/2,lent);
df=(f(2)-f(1));
fG=f*10^-9;tps=t*10^12;%GHz
scale=1;


%% -- Talbot Conditions -------------- %%%%%%%%



%% Talbot Parameters

q = 10; p = 1;
s = generateSparameter(p,q);

bitRateAim =1.1e9;                                      % Aim for bit rate; adjust to commensurate with dt.
nPeakPerBit=5;

% Adjust the phase levels of the Talbot phase to fit with the bit rate
tqAim=1/(nPeakPerBit*bitRateAim);                      % tq is modified such that ts is commensurate with dt
tsAim=tqAim/q;
ndt_ts=round(tsAim/dt);
ndt_tq=q*ndt_ts;
ts=dt*ndt_ts;
tq=q*ts;
nus=1/tq;
nuq=nus*q;


nBitsTot=time_window*2/(tq*nPeakPerBit);

%% define sig

% %% Periodic Talbot Amplificaiton
og_sig=generatePeriodicGauss(t,ts,dt*100,1).*superGauss(0,8*tq,15,t,0); %og_sig_f=circshift(og_sig_f,round(nuR/2/df));
  og_sig=circshift(og_sig,-round(ts/dt/4));
  
 %%% T-TAI
%    og_sig=ones(lent,1)';%genDiscPhase(lent,tq,2,dt,[0 1 ],1);ones(1,lent);
  og_sig=superGauss(0,8*tq,5,t,0);%genDiscPhase(lent,tq,2,dt,[0 1 ],1);ones(1,lent);


og_sig_f=nfft(og_sig,dt);
% % % % 

%% --  Phase Modulation ---------------------------------- %%%%%%%%

%% SpectRal phase filtering

%% T-TAI



% 
% Discrete spectral phase
%  discSpec=wrapTo2Pi(-pi*(m-1)/m*(0:m-1).^2);
% [phaseGVD]=genDiscPhase(lent,nus,m,df,discSpec,1);
% a=1;

% % % % %  discSpec=(p/m*pi*(0:m-1).^2);
 phi2=p*m*(tq/m)^2/(2*pi);
 phaseGVD=phi2/2*(2*pi*f).^2;
 phi2perKm=   2.1823e-23;

%% Temporal phase modulation




% % %T-TAI

 GV=wrapTo2Pi(s/m*pi*(0:m-1).^2); 
[phaseTemporal]=genDiscPhase(lent,ts,m,dt,GV,0); 
phaseTemporal=real(filtSG_tf(phaseTemporal,t,f,round(60e9/df),10,1));

% Sometimes need to shift the temporal signal a bit
phaseTemporal=circshift(phaseTemporal,-round(36.2e-12/dt)+round(tq/2/dt));




%% Talbot Magic - Dispersion then temporal Phase modulation (S-TAI) - %%%%%

% T-TAI

dispersed=og_sig.*exp(1j*phaseTemporal);
dispersed_f=nrmd_fft(dispersed,dt,scale);


spectrumRaw=dispersed_f.*exp(1j*phaseGVD);
temporalRaw=nrmd_ifft(spectrumRaw,Fs,scale);


%%%%%%
%% Find Peaks using a constant val signal
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







%% PlottingSection --------------------------------------------- %%%%%%%%%%
FS=18;

%Plotting parameters
% xlimf=[-2*ogBW 2*ogBW]*2*10^-9;

xlimt=[-200 200];%['auto'];%[-400 400];
xlimf=[-900 0]
xlimf=[-800 800]


xlimt=['auto'];%['auto'];%[-400 400];
xlimf=['auto'];
figure
subplot(2,1,2)


nrmFac=1;%max(abs(og_sig).^2);
plot(tps,abs(og_sig).^2/nrmFac,'DisplayName','Input')
hold on
ylabel('Intensity')
plot(tps,(abs(temporalRaw).^2)/nrmFac,'DisplayName','Temporal S-TAI')
hold on
% plot(tps,abs(temporal_OFT_pmOff).^2/nrmFac,'-b','DisplayName','OFT (pm off)')
% hold on
% plot(tps,abs(temporal_OFT_pmOn).^2/nrmFac,'-r','DisplayName',['OFT (pm on, amp: ' num2str(OFTamp)])
ylabel('Intensity')
yyaxis right

plot(tps,phaseTemporal,'DisplayName','Temporal Phase')
xlabel('time (ps)')
xlim(xlimt)
set(gca,'FontSize',FS)
legend('Show')
% -- SpectRal Plot --------------------------------%

subplot(2,1,1)
hold on
plot(fG,abs(og_sig_f).^2,'DisplayName','Input')
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




function create_AWGfile(s,m,gaussVals)

filename=['s',num2str(s),'m',num2str(m),'_',date];

Pm=gaussVals;

%Pm=repelem(Pm,4)

nrmd_Pm=Pm/max(Pm);
% figure;
% plot(nrmd_Pm,'o')
Signal_AWG7122C=nrmd_Pm;


%  output_file = ['AWG7122C' filename '.txt'];
%  save ( output_file,'Signal_AWG7122C', '-ascii')
%  fid=fopen(output_file,'w');
%  fprintf(fid,'%f\n',Pm);
%  fclose(fid);
%  
 output_file = ['AWG7122C_norm_' filename '.txt'];
  save (  output_file,'Signal_AWG7122C', '-ascii')
 fid=fopen(output_file,'w');
 fprintf(fid,'%f\n',Signal_AWG7122C);
 fclose(fid);

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

function y=generatePeriodicGauss(x,Tr,fwhm,center)


sigma=fwhm/(2*sqrt(2*log(2)));

dx=mean(diff(x));
nSingle=round(Tr/dx);
nReps=ceil(numel(x)/nSingle);
single=singleGauss(sigma,0,linspace(-nSingle/2*dx,nSingle/2*dx,nSingle),0);
y=repmat(single,1,nReps);
y=y(1:numel(x));

if center==1
    lent=numel(x)
   cent=round(lent/2);
   left=mod(cent,nSingle);
   
   shift=left-round(nSingle/2);
   y = circshift(y,shift);
   
end

end

function waveform=superGauss(C,t0,m,xs,center)
waveform=exp(-(1+1j*C)/2*((xs-center)/t0).^(2*m));
end

function y = singleGauss(std,center,xs,prop)


if prop==1
    y=1/(sqrt(2*pi*std^2))*exp(-(xs-center).^2/(2*std^2));
else
    y=exp(-(xs-center).^2/(2*std^2));
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

function [sigFilt,sig_f,wind,sigFilt_f]=filtSG_tf(sig,t,f,nBW,m,displayPlot)
% filtSG_tf(sig,t,f,nBW,m,displayPlot)


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

function cleanedup_angle=clean_angle(sig)
%this function samples the angle to points with a non-zero abs value.
%The angle is then restricted to the 0 to 2*pi domain
bare_angle=angle(sig);
%Sampled_sig returns a signal with a value of one whereever sig>1
tol=0.2;
sampled_sig=abs(sig)>tol*max(abs(sig));
%select only the angle where sig>tol*maxVal
sampled_angle=sampled_sig.*bare_angle;

% if min(sampled_angle<0)
%     sampled_angle=sampled_angle+2*pi;%make all values positive. 
% end
cleanedup_angle=mod(sampled_angle,(2*pi));
end
