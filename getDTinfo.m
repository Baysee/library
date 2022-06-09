function [x,y,inds]=getDTinfo(fig)

datacursormode on
dcm_obj = datacursormode(fig);


info_struct = getCursorInfo(dcm_obj);


nTips=numel(info_struct);
x=zeros(1,nTips); y=x; inds=x;


for i=1:nTips
[dati]=info_struct(i).Position;
x(i)=dati(1); y(i)=dati(2);
inds(i)=info_struct(i).DataIndex;

end
end