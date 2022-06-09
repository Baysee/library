function convertPtuMatBatch()

a=dir('*.ptu')
b=struct2cell(a);
fnames=b(1,:);

for i=1:numel(fnames)
    convert_HH_mPTU2mat(fnames{i})
end
end