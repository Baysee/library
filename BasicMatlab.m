% % Basic Matlab program


%% Define time and frequency vectors
lent=2^19;
time_window=5000; %picosecond. 
t=linspace(-time_window/2,time_window/2,lent); dt=t(2)-t(1);
Fs=1/dt; f=linspace(-Fs/2,Fs/2,lent); df=(f(2)-f(1));
fG=f*10^3; % GHz
tns=f*1e-3; % ns

%% Generate Signal Under Test (SUT)
SUT_f=superGauss(0,0.5,20,f,0);
SUT_t=nifft(SUT_f,Fs);

%% Plot Signal

figure;
subplot(2,1,1)
plot(t,abs(SUT_t).^2);
xlabel('Time (ns)'); ylabel('Intensity (n.u.)')
xlim([-10 10])
subplot(2,1,2)
plot(fG,abs(SUT_f).^2);
xlabel('Frequency (GHz)'); ylabel('Power (n.u.)')
xlim([-1500 1500])

 


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

function waveform=superGauss(C,t0,m,xs,center)
waveform=exp(-(1+1j*C)/2*((xs-center)/t0).^(2*m));
end