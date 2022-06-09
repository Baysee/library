function waveform=superGauss(C,t0,m,xs,center)
waveform=exp(-(1+1j*C)/2*((xs-center)/t0).^(2*m));
end