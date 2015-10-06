clear all; close all; clc;

X = dlmread('~/Dropbox/Cloud_Waldir_Gustavo/GPeixoto/dat/PCMn-x.dat');
Y = dlmread('~/Dropbox/Cloud_Waldir_Gustavo/GPeixoto/dat/PCMn-y.dat');
Z = dlmread('~/Dropbox/Cloud_Waldir_Gustavo/GPeixoto/dat/PCMn-Z.dat');

samplex = X(:,1); PCMnx = X(:,2);
sampley = Y(:,1); PCMny = Y(:,2);
samplez = Z(:,1); PCMnz = Z(:,2);

CX = gray2RGBColorMap(PCMnx);
cx = colormap(CX);

CY = gray2RGBColorMap(PCMny);
cy = colormap(CY);

CZ = gray2RGBColorMap(PCMnz);
cz = colormap(CZ);

axis off

hx = colorbar;
set(hx,'YTickLabel',samplex)
set(gcf,'color','None');
export_fig(hx,'../figs/mcpnx.pdf');

hy = colorbar;
set(hy,'YTickLabel',sampley)
set(gcf,'color','None');
export_fig(hy,'../figs/mcpny.pdf');

hz = colorbar;
set(hz,'YTickLabel',samplez)
set(gcf,'color','None');
export_fig(hz,'../figs/mcpnz.pdf');
