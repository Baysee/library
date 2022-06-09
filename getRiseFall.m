function [rise,fall]=getRiseFall(powerf_linNorm,percent)
rise=find(powerf_linNorm>percent,1); fall=find(powerf_linNorm>percent,1,'last');