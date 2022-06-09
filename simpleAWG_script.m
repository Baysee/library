function simpleAWG_script(vals,Srate,m,Os)

nrmd_pm=vals/max(vals);
Signal_AWG7122C=nrmd_pm;

SrateEff=Srate;
Srate=Os*Srate;
    OSopt=['-Os' num2str(Os) '-SrateEff-' num2str(SrateEff,4)];

%   Repeat signal to make it longer  
    sigLen=numel(Signal_AWG7122C);
    sigLen_t=sigLen/(Srate*10^9);
    
    targetLen_t=500e-9;
    nReps=floor(targetLen_t/sigLen_t);
    
    Signal_AWG7122C=repmat(Signal_AWG7122C,1,nReps);


marker=zeros(1,numel(Signal_AWG7122C));marker(1)=1;
Signal_AWG7122C=[Signal_AWG7122C;marker;marker]';

%  output_file = ['AWG7122C' filename '.txt'];
%  save ( output_file,'Signal_AWG7122C', '-ascii')
%  fid=fopen(output_file,'w');
%  fprintf(fid,'%f\n',Pm);
%  fclose(fid);
dt=datetime;
date=[num2str(dt.Day),'-',num2str(dt.Month),'_',num2str(dt.Hour),'h',num2str(dt.Minute)];
filename=[num2str(Srate,4),'Gbs-m-',num2str(m),OSopt,'_',date,'.txt'];
csvwrite(filename,Signal_AWG7122C)

%  output_file = [filename '.txt'];
%   save (  output_file,'Signal_AWG7122C', '-ascii')
%  fid=fopen(output_file,'w');
%  fprintf(fid,'%f\n',Signal_AWG7122C);
%  fclose(fid);
 
 end