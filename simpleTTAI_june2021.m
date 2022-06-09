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

% p spectRal phase Parameter. s Temporap phase parameter
% m multiplicationparameter
m=17;p=2;%% keep p=1 for now (need to update next line)
s=generateSparameter(p,m);%mod(1+m*mod(m,2),2*m);

%  phi2=p*m*(ts)^2/(2*pi);
%  phaseGVD=phi2/2*(2*pi*f).^2;
 phi2perKm=   2.1823e-23;
 
 ts=sqrt((600*phi2perKm)*(2*pi)/(p*m))
% ts=1/12e9;%nts*dt;%tR;
tq=m*ts;

% tq=m*(1/20.02e9);%1e-9;%1/1e9;
% ts=tq/m;


nus=1/tq;
nuq=nus*m;

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
 phi2=p*m*(ts)^2/(2*pi);
 phaseGVD=phi2/2*(2*pi*f).^2;
 phi2perKm=   2.1823e-23;
phi2/phi2perKm
%% Temporal phase modulation




% % %T-TAI

     GV=wrapTo2Pi(s/m*pi*(0:m-1).^2); 
[phaseTemporal]=genDiscPhase(lent,ts,m,dt,GV,0); 
phaseTemporal=real(filtSG_tf(phaseTemporal,t,f,round(20e9/df),10,1));

% Sometimes need to shift the temporal signal a bit
phaseTemporal=circshift(phaseTemporal,-round(36.2e-12/dt)+round(tq/2/dt));




%% Talbot Magic - Dispersion then temporal Phase modulation (S-TAI) - %%%%%

% T-TAI

dispersed=og_sig.*exp(1j*phaseTemporal);
dispersed_f=nrmd_fft(dispersed,dt,scale);


spectrumRaw=dispersed_f.*exp(1j*phaseGVD);
temporalRaw=nrmd_ifft(spectrumRaw,Fs,scale);




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



