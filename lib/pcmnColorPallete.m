%% pcmnColorPallete
%{  
    Prints to PDF a color pallete associating 
    a parameter value to their respective 
    sample number.

    author: Gustavo Peixoto
    date: Oct 7, 2015

    required: 
        - gray2RGBColorMap: this function builds a new RGB colormap
                            from the scalar parameter data.
                            (Gibbon code, by Kevin Moerman)

        - export_fig: this function allows to export 
                      the colorbar to figure and produce a
                      .pdf with reduced printed area. 
                      Differs from the standard 'print'.
                      (Available on Mathworks, by Yair Altman)
                      
                      

%}

clear all; close all; clc;

%% Load Files of PCMn

X = dlmread('~/Dropbox/Cloud_Waldir_Gustavo/GPeixoto/dat/PCMn-x.dat');
Y = dlmread('~/Dropbox/Cloud_Waldir_Gustavo/GPeixoto/dat/PCMn-y.dat');
Z = dlmread('~/Dropbox/Cloud_Waldir_Gustavo/GPeixoto/dat/PCMn-Z.dat');

xname = '../figs/mcpnPalleteX.pdf';
yname = '../figs/mcpnPalleteY.pdf';
zname = '../figs/mcpnPalleteZ.pdf';

%% Sample Numbers and Parameter Values

samplex = X(:,1); PCMnx = X(:,2);
sampley = Y(:,1); PCMny = Y(:,2);
samplez = Z(:,1); PCMnz = Z(:,2);

%% Color Pallete X 

CX = gray2RGBColorMap(PCMnx);
cx = colormap(CX);
axis off
set(gcf,'color','None');
sx = num2str(samplex);
sx = cellstr(sx);
hx = colorbar;
set(hx,'YTick',1:length(samplex));
set(hx,'YTickLabel',sx');
export_fig(hx,xname); 


%% Pallete Y 

CY = gray2RGBColorMap(PCMny);
cy = colormap(CY);
axis off
set(gcf,'color','None');
sy = num2str(sampley);
sy = cellstr(sy);
hy = colorbar;
set(hy,'YTick',1:length(sampley));
set(hy,'YTickLabel',sy);
export_fig(hy,yname);

%% Pallete Z 

CZ = gray2RGBColorMap(PCMnz);
cz = colormap(CZ);
axis off
set(gcf,'color','None');
sz = num2str(samplez);
sz = cellstr(sz);
hz = colorbar;
set(hz,'YTick',1:length(samplez));
set(hz,'YTickLabel',sz);
export_fig(hz,zname);



