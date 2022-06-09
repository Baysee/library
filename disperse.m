% function disperse()% [f,out]=disperse(t,in,scale)
t=linspace(-120e-9,120e-9,2^19);
lent=numel(t);


scale=1;
%t must be in picosecond
dt=t(2)-t(1);
Fs=1/dt;
f=linspace(-Fs/2,Fs/2,numel(t));
x=t; Tr=4e-9; fwhm=dt*3; center=0;
in=generatePeriodicGauss(x,Tr,fwhm,center).*superGauss(0,20*Tr,8,t,0);

envelope=zeros(1,lent);
envi=[1 1 1 4 4 4 1 1 1 6 6 6 1 1 1 4 4 4 1 1 1 6 ].^2;
envelope=repmat(envi,1,ceil(lent/22));
envelope=envelope(1:lent);
env_f=nfft(envelope);
10e-12/dt
%  in_f=superGauss(0,5000e9/2,10,f,00);

% in_f=in_f+(5*(rand(L,1)'+1i*rand(L,1)')).*SuperGauss(0,500e6,30,f,0);;
% % 
% % Tr=1/4.86e9;
% % singleL=linspace(-Tr/2,Tr/2,round(Tr/dt));
% % Gauss=singleGauss(2e-12,0,singleL,0);
% % nreps=ceil((max(t)-min(t))/Tr);
% % in=repmat(Gauss,1,nreps);
% % in=in(1:L);
%   in=nifft(in_f,Fs,scale);
%  delay=50e-12;%picosecond
%  delayind=round(delay/dt);in=in+[in(delayind:end),in(1:delayind-1)];
%  a=1;
 
% in=SuperGauss(0,0.1e-9,30,t,0);

% tsub=-2e-9:dt:2e-9;
% decimatedM=in((L-numel(tsub))/2:(L+numel(tsub))/2-1)/max(in)*10^-6;
% nreps=floor(L/length(decimatedM));
% 
% fulldM=repmat(decimatedM,1,nreps);
% fulldM=[fulldM,zeros(1,L-length(fulldM))];
% 
% in=in/max(in)*sqrt(20*10^-3);
 in_f=nfft(in,dt,scale);


disp=120e3*17.3e-6;%Now using FBG equivalent to 120 km of SMF @ 17.3e-6 s/m^2
%Units: 17.3 ps/(nm km)=17.3 (e-12 s)/(e-9 m e3 m)=e-6 s/,
% lambda=1550e-9;%m
% speedoflight=299792458;%m/s
% phi2=(lambda^2/(2*pi*speedoflight))*disp;%1981e-3/(2*pi*speedoflight)*lambda^2;%*22e-24*120;%-(lambda^2/(2*pi*speedoflight))*disp;%Experimental phi2 in s^2
% %   phi2=22e-24*(235);%-(lambda^2/(2*pi*speedoflight))*disp;%Experimental phi2 in s^2
phi2=Tr^2/(2*pi)*1/40;



%% Dispersion
 phaseMod=exp(1j*phi2*(2*pi*f).^2/2);
dispersed_f=phaseMod.*in_f.*exp(1j*env_f);


%% Time lens
%  phaseMod=exp(1j*phi2*(2*pi*f).^2/2);
% dispersed_f=in_f;



%% Find FWHM


out=nifft(dispersed_f,Fs,scale);%+2*superGauss(0,1e-12,10,t,1.79e-9);

figure;plot(t,abs(out).^2)
figure;plot(f,abs(dispersed_f).^2)
% dispersed_f=nfft(out,dt,scale);
% 
% rise_out=find(abs(out).^2>max(abs(out).^2)/2,1);
% fall_out=find(abs(out).^2>max(abs(out).^2)/2,1,'last');
% FWHM_out=t(fall_out)-t(rise_out);
% 
% rise_in=find(abs(in).^2>max(abs(in).^2)/2,1);
% fall_in=find(abs(in).^2>max(abs(in).^2)/2,1,'last');
% FWHM_in=t(fall_in)-t(rise_in);
% 
% 
% avgin=mean(abs(in).^2);
% 
% avgout=mean(abs(out).^2);
% 



%% Plot
figure
subplot(211)
plot(t,abs(in).^2,'DisplayName',['input waveform (power=',num2str(avgin), '), FWHM ', num2str(FWHM_in*1e12), ' ps'] )
hold on
% legend('Show')
% plot(t,SuperGauss(0,4e-9,30,t,0),'DisplayName','1/250mhz')
 yyaxis right
% subplot(212)
plot(t,abs(abs(out)).^2,'DisplayName',['output waveform  (power=',num2str(avgout), '), FWHM ', num2str(FWHM_out*1e12), ' ps'])
xlabel('time (s)')
legend('Show')
subplot(212)
plot(f,in_f,'DisplayName','Input Spectrum')
plot(f,abs(dispersed_f).^2,'DisplayName','dispersed spectrum')
hold on
plot(f,real(dispersed_f),'DisplayName','dispersed spectrum')
% yyaxis right
%  plot(f,unwrap(clean_angle(dispersed_f)),'DisplayName','real part phase mod')
ylabel('real part')
xlabel('frequency (Hz)')
legend('Show')
hold off
