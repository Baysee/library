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


%% -- Talbot Conditions; GVD and Temporal Phase set -------------- %%%%%%%%

% p spectRal phase Parameter. s Temporap phase parameter
% m multiplicationparameter
m=12;p=1;%% keep p=1 for now (need to update next line)
% % %  [s,c,GaussVals]=find_Gauss_para(p,m,1);
s=generateSparameter(p,m);%mod(1+m*mod(m,2),2*m);
% s=-(q-1);
% s=m+1;18.4-14.
% s=-(m-1)

%-- Define GVD via sampling parameters for temporal phase modulation -----%
% % 
% nuR=12*18e9;
% tR=1/(nuR);%Single temporal level duration


%% T-TAI
%  nts=round(1/(12*18e9)/dt);

ts=1/120e9;%nts*dt;%tR;
tq=m*ts;

% tq=m*(1/20.02e9);%1e-9;%1/1e9;
% ts=tq/m;


nus=1/tq;
nuq=nus*m;

%% S-TAI
% 
% nuR=1.0e9;
% % nus=2e9;
% nuq=m*50e9;%nus*m;
% nus=nuq/m;
% 
% ts=1/nuq;
% tq=m*ts;

% %% define sig
 BW=0.1e9;
% 
% og_sig_f=superGauss(0,BW,4,f,-200e9)+superGauss(0,BW,4,f,200e9);
% og_sig_f=superGauss(0,BW,4,f,0);
% 
% % 
% og_sig_f=generatePeriodicGauss(f,nuq,df*5,1); %og_sig_f=circshift(og_sig_f,round(nuR/2/df));
%   og_sig=nrmd_ifft(og_sig_f,Fs,scale);

% trunc=superGauss(0,time_window/10,5,t,0);
% % trunc=superGauss(0,tq*5,5,t,0)+superGauss(0,tq*2,1,t,tq*4)+superGauss(0,tq*2,1,t,-tq*5)/2;
% og_sig=trunc;%generatePeriodicGauss(t,tq,dt*2,1);%.*trunc; %og_sig_f=circshift(og_sig_f,round(nuR/2/df));
% og_sig_f=nrmd_fft(og_sig,dt,scale);
% 

% og_sig_fIni=generatePeriodicGauss(f,nus,nus/20,0);
% [pks,locs]=findpeaks(abs(og_sig_fIni));
% fMid=locs(round(numel(locs)/2));
% distFreq=400e9;
% % 
% [S_OS]=digiModFunc('SP-PAM',2,6,1,0);
% [og_sig]=genDiscPhase(lent,tq,2^6-1,dt,S_OS,0); 
% og_sig=real(filtSG(og_sig,round(50e9/df),1,1));
% widthf0=7*nus/2;
% og_sig_f=(superGauss(0,widthf0,20,f,fMid-distFreq)+superGauss(0,widthf0,20,f,fMid+distFreq-nus)).*og_sig_fIni;
% %  og_sig=nrmd_ifft(og_sig_f,Fs,scale);
% % % % 

og_sig=generatePeriodicGauss(t,ts,dt*100,1).*superGauss(0,8*tq,15,t,0); %og_sig_f=circshift(og_sig_f,round(nuR/2/df));

%    og_sig=ones(lent,1)';%genDiscPhase(lent,tq,2,dt,[0 1 ],1);ones(1,lent);
%   og_sig=superGauss(0,5e-9,5,t,0);%genDiscPhase(lent,tq,2,dt,[0 1 ],1);ones(1,lent);
  og_sig=circshift(og_sig,-round(ts/dt/4));
% % % % % 
og_sig_f=nfft(og_sig,dt);
% % % % 

%% --  Phase Modulation ---------------------------------- %%%%%%%%

%% SpectRal phase filtering

% S-TAI
% discSpec=wrapTo2Pi(s/m*pi*(0:m-1).^2); 
%  discSpec=wrapTo2Pi(p*(1+m*mod(m,2))/m*pi*(0:m-1).^2);
% [phaseGVD]=genDiscPhase(lent,nus,m,df,discSpec,0);
% phaseGVD=circshift(phaseGVD,round(-50e9/df));
% %  discSpec=(p/m*pi*(0:m-1).^2);

% nReps=150;
%  discSpec=(p*(1+m*mod(m,2))/m*pi*(0:nReps*m-1).^2);
%  discSpec=[fliplr(discSpec(2:end)),discSpec];
%  [phaseGVD]=genDiscPhase(lent,nuR,nReps*m-1,df,discSpec,1);

 
%  discSpec=circshift(discSpec,-1);
%  [phaseGVD]=genDiscPhase(lent,nuR,m,df,discSpec,1);
% % 
% % % % phaseGVD=circshift(phaseGVD,-round(tR/dt/2));
% %  phi2=p/m*2*pi/(2*pi*nus)^2;
%   phi2perKm=   2.1823e-23;
% %   phi2/phi2perKm
%  
%  
%  phi2= 120*phi2perKm;%2.1823e-22;%8.7292e-22;%2.5465e-21
%  phi2/phi2perKm
%  nus=1/((2*pi)*sqrt(phi2*m/p/(2*pi)))
%  nuq=m*nus;
% %  nuq=119.04e9/2
% % ;
% %  nus=nuq/m;%
% ts=1/nuq;
% tq=m*ts;
% % phi2=00.98*phi2
% % % %  phi2=p/m*2*pi*(1+m*mod(m,2))/(2*pi*nuR)^2;
% TOD=0;%1e-36*(2*pi*f).^3;
% phaseGVD=phi2/2*(2*pi*f).^2+TOD;
%  
% %   [~,mdS]=min(phaseGVD);
%  [~,mphi]=min(phi2/2*(2*pi*f).^2);
%  shift=mphi-mdS-round(nuR/df/2);
% phaseGVD=circshift(phaseGVD,shift);   




% 
% 
%% T-TAI



% 
% discSpec=wrapTo2Pi(p*(1+m*mod(m,2))/m*pi*(0:m-1).^2);
 discSpec=wrapTo2Pi(-pi*(m-1)/m*(0:m-1).^2);
[phaseGVD]=genDiscPhase(lent,nus,m,df,discSpec,1);
a=1;

% % % % %  discSpec=(p/m*pi*(0:m-1).^2);
%  phi2=p*m*(tq/m)^2/(2*pi)
%  phaseGVD=phi2/2*(2*pi*f).^2;
%  phi2perKm=   2.1823e-23;
% %  phi2=phi2perKm*240;
% %   phaseGVD=phi2/2*(2*pi*f).^2;
%  phi2/phi2perKm
%  a=1
 
% % % % 
% disp=-(595e3)*17e-6;%-10.08;%-(120e3)*17.3e-6-1.3;%-9.956;%-10.08;%120e3*19.46e-6;%10.08;%(590e3)*17.3e-6;%19.64e-6 gives optimal value for 120km;
% % % %disp=-10.013;%-1.55;%247.92e-3;%125.5e3*17.3e-6;%Now using FBG equivalent to 120 km of SMF @ 17.3e-6 s/m^2
% % % %Units: 17.3 ps/(nm km)=17.3 (e-12 s)/(e-9 m e3 m)=e-6 s/,
% lambda=1550e-9;%m
%  speedoflight=299792458;%m/s
%  phi2= -lambda^2*disp/(2*pi*speedoflight);
% % disp=-phi2*2*pi*speedoflight/lambda.^2;
% kmFiber=disp/(17e-6)/1e3 %div by 17 ps/(nm*km)=10^-12/(10^-9*10^3)=10^-6
% phaseGVD=phi2/2*(2*pi*f).^2;
% 
% ts=sqrt(abs(2*pi*phi2/(m)));
% tq=ts*m;
% 1/ts
% % [95.6, 100.2]; [107.2 112.2] ; [115.2 120.2]]
% %  ts=1/36.68e9;
% tq=ts*m;
% 
% nus=1/tq; 
% nuq=1/ts
% a=1;


% phaseGVD=(0.99*phi2)/2*(2*pi*f).^2;

%% Temporal phase modulation
%S-TAI
% 
% GV=wrapTo2Pi(s/m*pi*(0:m-1).^2); 
% GV=wrapTo2Pi(p*(1+m*mod(m,2))/m*pi*(0:m-1).^2);
% % 
% [phaseTemporalRaw]=genDiscPhase(lent,ts,m,dt,GV,0);
%  phaseTemporal=phaseTemporalRaw;%^filtSG(phaseTemporalRaw,round(100e9/df/2),1,1)*max(phaseTemporalRaw);
% % phaseTemporal=circshift(phaseTemporal,round(50e9/df));
% phaseTemporal=circshift(phaseTemporal,round(tS/dt/3));
% % phaseTemporal=pi*s/m*(t/tS).^2;
% a=1;



% % %T-TAI
% % % 
% GV=wrapTo2Pi(s/m*pi*(0:m-1).^2); 
% % GV=wrapTo2Pi(p/m*pi*(0:m-1).^2); 
 GV=wrapTo2Pi(s/m*pi*(0:m-1).^2); 
%  GVcorr=GV/max(GV)
%  GVmeas=[0 0.8 0.16 1.16 241/290 181/290 157/290 157/290 192/290 246/290 1 0.1 0.75]
 
%  create_AWGfile(s,m,GVcorr+(GVcorr-GVmeas));
  
%  GV=wrapTo2Pi(s/m*pi*((0:m-1)-round(m/2)).^2); 
% % % % 
% % GV=wrapTo2Pi(-pi*(m-1)/m*(0:m-1).^2);
% GV=circshift(GV,10);
[phaseTemporal]=genDiscPhase(lent,ts,m,dt,GV,0); 
% phIni=phaseTemporal;%
% phaseTemporal=phIni;%real(filtSG(phaseTemporal,round(80e9/df),1,1)*max(phaseTemporal));
% figure;plot(f*1e-9,abs(nfft(phIni)));hold on;plot(f*1e-9,abs(nfft(phaseTemporal)));
% GV=round(GV/(2*pi)*2^8)/2^8*2*pi;; 
% phaseTemporal=real(filtSG(phaseTemporal,round(800e9/df),1,1));%*max(phaseTemporal));
phaseTemporal=circshift(phaseTemporal,-round(36.2e-12/dt)+round(tq/2/dt));
b=real(filtSG_tf(phaseTemporal,t,f,round(60e9/df),10,1));




%% Talbot Magic - Dispersion then temporal Phase modulation (S-TAI) - %%%%%
%%S-TAI 

% dispersed_f=og_sig_f.*exp(1j*phaseGVD);
% dispersed=nrmd_ifft(dispersed_f,Fs,scale);
% 
% temporalRaw=dispersed.*exp(1j*phaseTemporal);
% spectrumRaw=nrmd_fft(temporalRaw,dt,scale);


% dispersed_f=og_sig_f.*exp(1j*phaseGVD);
% dispersed=nrmd_ifft(dispersed_f,Fs,scale);
% 
% temporalRaw_pmoff=dispersed;%.*exp(1j*phaseTemporal);
% spectrumRaw_pmoff=nrmd_fft(temporalRaw_pmoff,dt,scale);


% Frequency to time mapping

% % %  phi2=p/m*2*pi*(1+m*mod(m,2))/(2*pi*nuR)^2;
% phiOFT=(600*phi2perKm);
% phaseGVD_OFT=phiOFT/2*(2*pi*f).^2;
% 
% spectrumRaw_OFT_pmOn=spectrumRaw.*exp(1j*phaseGVD_OFT);
% temporal_OFT_pmOn=nrmd_ifft(spectrumRaw_OFT_pmOn,Fs,scale);
% 
% 
% spectrumRaw_pmoff=dispersed_f.*exp(1j*phaseGVD_OFT);
% temporal_OFT_pmOff=nrmd_ifft(spectrumRaw_pmoff,Fs,scale);
% 
% OFTamp=max(abs(temporal_OFT_pmOn).^2)/max(abs(temporal_OFT_pmOff).^2)
% sum(abs(temporal_OFT_pmOn).^2)/sum(abs(temporal_OFT_pmOff).^2)
% 
% % Verify validity of OFT
% t0=2*abs(t(find(abs(temporalRaw).^2>(max(abs(temporalRaw).^2/2)),1)))
% ratioOFT=abs(phiOFT)/(t0^2/pi)

% %  

% spectrumRaw=spectrumRaw.*(exp(-1j*phaseGVD));

% T-TAI
% 
dispersed=og_sig.*exp(1j*phaseTemporal);
dispersed_f=nrmd_fft(dispersed,dt,scale);


spectrumRaw=dispersed_f.*exp(1j*phaseGVD);
temporalRaw=nrmd_ifft(spectrumRaw,Fs,scale);



% 
% figure;
% phases=zeros(m,numel(phaseGVD));
% modSig=phases;
% for i=1:m
%     
%    phases(i,:)=-discSpec(i)+circshift(phaseGVD,i*round(nus/df));
%    modSig(i,:)=og_sig_f.*exp(1j*phases(i,:));
%    subplot(m+1,1,i)
%    plot(f,abs(og_sig_f));
%    yyaxis right
%    plot(f,phases(i,:))
% end
% 
% subplot(m+1,1,m+1)
% plot(f,abs(sum(modSig)))


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



