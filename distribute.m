function [varargout]=distribute(varargin)
% Assign inputs to outputs

% if numel(varargout)~=numel(varargin)
%     error('Must have same number of inputs as outputs in distribute function!')
% end

for i=1:numel(varargin)
    varargout{i}=varargin{i};
end
end