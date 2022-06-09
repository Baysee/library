function ifftout=nrmd_idct(sig,Fs,scale)
%%%Need to first ifftshift to swap for ifft (Matlab assumes swaped fft)
%normalize according to amplitude
unnormdifft=fftshift(idct(ifftshift(sig)));
if scale==1
ifftout=Fs*unnormdifft;
else if scale==2
    ifftout=1/sqrt(length(sig))*unnormdifft;
else
ifftout=1/(max(abs(unnormdifft)))*unnormdifft;
    end
end
%ifftout=df*unnormdifft;
end