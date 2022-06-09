function fftout=nrmd_fft(sig,dt,scale)
%fft gives the swaped spectrum.
unnormdfft=fftshift(fft(ifftshift(sig)));

if scale==1
fftout=dt*unnormdfft;
    else if scale==2
    fftout=1/sqrt(length(sig))*unnormdfft;
        else
        fftout=1/(max(abs(unnormdfft)))*unnormdfft;
        end
end
%fftout=dt*unnormdfft;
end