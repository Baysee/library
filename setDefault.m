function [AFS,LFS,TFS,DFL,PosFig,col1,col2,col3,col4,col5]=setDefault()

% SETTINGS 
DimEcran = get(0,'ScreenSize');
Dx = 0.60; Dy = 0.8;
PosFig  = [(1-Dx)*DimEcran(3),(0.9-Dy)*DimEcran(4),Dx*DimEcran(3),Dy*DimEcran(4)];

AFS = 12;   LFS = 18;   TFS = 18; DFL=1;  %Axis, label and title font size



col1 = [5,112,0]/256;[7,160,32]/300;%Green
col2 = [0,0,0];%Gray
col3 = [30/256,48/256,160/256];%blue
col4 = [231,190,90]/256;%Yellow
col5 = [0 0 0];%Black
end