clear all;
clc;
format_spec = '%g';
numcell = 0.0; L_n = 0.0;
identifier = {'# number of elements'};
% get numcell
file_name = 'inputmesh.dat';
file_ID = fopen(file_name);
while ~feof(file_ID)
    readLine = fgetl(file_ID);
    if contains(readLine, char(identifier(1)))
        k = strfind(cellstr(readLine), ' ');
        numcell = str2num(readLine(1:k{1}(1)-1));
    end
end
fclose(file_ID);
% get L_n
file_name = 'inputdata.dat';
L_n = textread(file_name, '%f', 1, 'headerlines', 1);
% get temperature data
file_name = 'ResultTempCell.dat';
file_ID = fopen(file_name);
content = fscanf(file_ID, format_spec, [4, numcell]);
x = content(2, :); y = content(3, :); z = content(4, :);
fclose(file_ID);
% Interpolate and get midline temperature
radius = 0.0;
if numcell < 500
    radius = 0.07;
else if numcell < 1500
        radius = 0.04;
    else if numcell < 2000
            radius = 0.03;
        else if numcell < 4000
                radius = 0.025;
            else if numcell < 8000
                    radius = 0.02;
                end
            end
        end
    end
end
step = 0.01; num_point = 1 / step + 1;
xi(num_point) = 1.0; zi(num_point) = 0.0;
x = x./L_n; y = y./L_n;
output_file = 'Output.dat';
fileOP = fopen(output_file, 'w');
for i = 1:num_point
    temp = 0.0; count = 0; sum = 0.0; xi(i) = (i - 1) * step;
    if xi(i) < 0.09 || xi(i) > 0.9
        dist = 3 * radius;
    else if xi(i) <= 0.2 || xi(i) >= 0.8
            dist = 2.5 * radius;
        else if xi(i) < 0.3 || xi(i) > 0.7
                dist = 2 * radius;
            else if xi(i) <= 0.38 || xi(i) >= 0.62
                    dist = 1.6 * radius;
                else
                    dist = radius;
                end
            end
        end
    end
    for k = 1:numcell
        d = sqrt((x(k)-xi(i))*(x(k)-xi(i))+(y(k)-0.5)*(y(k)-0.5));
        if d <= dist
            temp = temp + z(k) / d;
            sum = sum + 1 / d;
            count = count + 1;
        end
    end
    if count > 0
        temp = temp / sum;
        zi(i) = temp;
        fprintf(file_ID, '%g %g %g %g\n', xi(i), count, temp, dist);
    end
end
fclose(fileOP);
plot(xi, zi);
xlabel('x coordinates', 'Fontsize', 16);
ylabel('Midline Temperature(K)', 'Fontsize', 16);
title('Midline Temperature of the Mesh', 'Fontsize', 16);