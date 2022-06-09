function fNames=findFnames(targetDir,ext)

a=dir([targetDir '/*' ext]);
b=struct2cell(a);
fNames=b(1,:);

end
