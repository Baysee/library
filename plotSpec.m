


function plotSpec(SUT_t,tempPhase, nLevels,AnW_Len,AnW_Inc, nFreqs, dt,xf,xt,flim,tlim)
%
%Example call:
% 
% SUT_t=(E_t_GVD); 
% % SUT_t=(E_t); 
% AnW_Len=round(10*tq/dt);
% AnW_Inc=round(AnW_Len*0.7);
% nFreqs=4005*2;
% tempPhase=0;
% nLevels=0;
% % xf=abs(E_f_PM); xt=abs(E_t); flim='auto'; tlim='auto'; 
% xf=f; xt=t; flim=[-130 130]*10^9; tlim=[-13.6 13.6]*10^-9; 
% 
%     plotSpec(SUT_t,tempPhase, nLevels,AnW_Len,AnW_Inc, nFreqs, dt,xf,xt,flim,tlim)

%% Get Spectrogram data
%   [t_Spec,f_Spec,SUT_Spec]=QuickSpectrogram(SUT_t,AnW_Len,AnW_Inc,dt);
 SUT_f=nfft(SUT_t,dt);
  [SUT_Spec,f_Spec,t_Spec]=spectrogram(SUT_t,AnW_Len,AnW_Inc,nFreqs,1/dt,'centered','yaxis');
t_Spec=t_Spec-mean(t_Spec);
  
% spectrogram(y_t(inds2Plot),windowL,overlap,nFreqs,Fs,'centered','yaxis')
SpecPow=abs(SUT_Spec).^2;






%% Main plot

 nsp=5; %number of grids
 
figure;

%% Frequency Plot
subplot(nsp,nsp,1:nsp:(nsp^2-nsp))
% logSut=log(abs(SUT_f).^2); logSpec=log(sum(SpecPow,2));  %To plot in log
logSut=(abs(SUT_f).^2); logSpec=(sum(SpecPow,2));  % Plot linear
% 
%%% If plotting in log, normalize with - 
% plot(logSpec-max(logSpec),(f_Spec),'Marker','.')
% hold on
% plot(logSut-max(logSut),(xf));

% Plot linear
plot(logSpec/max(logSpec),(f_Spec),'Marker','.')
hold on
plot(logSut/max(logSut),(xf));

ylabel('Frequency'); xlabel('Power')
ylim(flim);


subplot(nsp,nsp,reshape((1:nsp:(nsp^2-nsp))+[1:nsp-1]',[1,(nsp-1)^2]))
imgPlot=log(SpecPow);
%% Plot imagesc
% imagesc(t_Spec,f_Spec,log(SpecPow));
cHigh=max(max(imgPlot));
cLow=cHigh-10;
CLIM = [cLow cHigh];
imagesc(t_Spec,f_Spec,imgPlot,CLIM );

set(gca, 'YDir','normal')


xlim(tlim); ylim(flim);

%% Time Plot
subplot(nsp,nsp,(nsp^2-(nsp-2)):nsp^2)

plot(t_Spec,(sum(SpecPow,1))/max(sum(SpecPow,1)),'Marker','.')
% yyaxis right;
hold on
plot(xt,abs(SUT_t).^2/max(abs(SUT_t).^2))
% yyaxis right;
% plot(xt,tempPhase);

% xlim([min(t_Spec) max(t_Spec)]);
xlim(tlim)
xlabel('Time');ylabel('Intensity')


 end