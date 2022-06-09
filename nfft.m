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