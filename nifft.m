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