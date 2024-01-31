function [] = saveFigure(folder, name)
% savefigure(folder, name) saves all open figures as a .fig file at the
% path 'folder\name.fig'.

h = findobj('type','figure');
savefig(flip(h),[folder '\' name '.fig'])
end

