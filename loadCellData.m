function varargout = loadCellData(fn)

s=load(fn);

fns=fieldnames(s);

for i=1:numel(fns)
    varargout{i}=s.(fns{i});
end