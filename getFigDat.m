function [xdata,ydata]=getFigDat(fig)

% fig = gcf;

% xDats=findobj(fig,'-property','XData');
% yDats=findobj(fig,'-property','YData');
lineObjs=findobj(fig,'type','line');
xdata = get(lineObjs, 'XData');
ydata = get(lineObjs, 'YData');
end
