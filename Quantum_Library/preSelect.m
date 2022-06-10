
function [inds1,data1_tm,inds1_tm,hist1_tm]=preSelect(tm1,hist1_inds,data1,indsDat1)
% tm1 = time_win{1};                                      % select counts within windows
% tm1(tm1<1)=[];
inds1 = ismember(hist1_inds,tm1);                       % Take only the counts there are within the preselections
data1_tm=data1(inds1,:);
inds1_tm=indsDat1(inds1);
[hist1_tm,~]=histcounts(data1_tm(:,3),0:2^15);
end
