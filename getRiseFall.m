function [rise,fall]=getRiseFall(powerf_linNorm,percent)
rise=find(powerf_linNorm>percent*max(powerf_linNorm),1); fall=find(powerf_linNorm>percent*max(powerf_linNorm),1,'last');