clear all;
clc;
format long;
formatSpec = '%g';
fid = fopen('inputmesh.dat');
i = 1;
identifier = {'# Elements', '# number of mesh points', '# number of elements'};
%type (strcat(meshName, '.mphtxt'));
% get data from comsol output file

while ~feof(fid)
    readLine = fgetl(fid);
    if strcmp(readLine, identifier(1))
        IEN = [];
        while(~isempty(readLine))
            readLine = fgetl(fid);
            B = str2num(readLine);
            IEN = [IEN; B];
            i = i+1;
        end
    end
    if strfind(readLine, char(identifier(2)))
        k = strfind(cellstr(readLine), ' ');
        numnode = str2num(readLine(1:k{1}(1)-1));
    end
    if strfind(readLine, char(identifier(3)))
        k = strfind(cellstr(readLine), ' ');
        numcell = str2num(readLine(1:k{1}(1)-1));
    end
end
fclose(fid);
IEN = IEN';
IEN = IEN + 1;
file_name = 'inputdata.dat';
L_n = textread(file_name, '%f', 1, 'headerlines', 1);
fileID = fopen('ResultTempNode.dat');
Coord = fscanf(fileID,formatSpec,[4, numnode]);
fclose(fileID);
fileID = fopen('ResultTempCell.dat');
TS = fscanf(fileID,formatSpec,[4, numcell]);
fclose(fileID);
TT = TS(4, :);

% Post-processing          
X = [Coord(2,IEN(1,:)); Coord(2,IEN(2,:)); Coord(2,IEN(3,:))];
Y = [Coord(3,IEN(1,:)); Coord(3,IEN(2,:)); Coord(3,IEN(3,:))];   
figure(1)
hold on
DX=patch(X,Y,TT);
axis equal;
axis off;
h=colorbar;
set(h,'FontSize',16);
title('Contour Plot using Cells', 'Fontsize', 16);

% Post-processing
T = Coord(4,:);
DX_elem = [T(IEN(1,:));T(IEN(2,:)); T(IEN(3,:))];
figure(2)
hold on
DX=patch(X,Y,DX_elem);
axis equal;
axis off;
h=colorbar;
set(h,'FontSize',16);
title('Contour Plot using Nodes', 'Fontsize', 16);

xi = TS(2, :)./L_n;
yi = TS(3, :)./L_n;
xi = xi'; yi = yi'; TT = TT';
[X,Y,Z] = griddata(xi, yi, TT, linspace(min(xi), max(xi))', linspace(min(yi), max(yi)), 'linear');
figure(3)
[C,H]=contourf(X, Y, Z);
clabel(C, 'Fontsize', 12, 'Color', 'red');
title('Isotherm Line Plot', 'Fontsize', 16);
