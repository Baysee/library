% function convertPtuMatBatch_dirs()


allFiles=dir('**')

root='D:\MATLAB_DataAnalysis\data\Crockett_QTAI20211130\QTAI_20211130\QTAI_exp\quantumData'


% 
% a=dir('*.ptu')
% b=struct2cell(a);
% fnames=b(1,:);

for i=1:numel(allFiles)
    
    folder=allFiles(i).folder;
    fn=allFiles(i).name;
    if numel(fn)<5
        continue
    end
    extension=fn(end-3:end)
    i
    folder
    if strcmp(extension,'.ptu')
    cd(folder)
     convert_HH_mPTU2mat(fn)
     cd(root)
end

end