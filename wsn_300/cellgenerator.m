clc
clear all

cell = zeros(200,3);
for i = 1 : 120
    cell(i,1) = i-1;
    cell(i,2) = ceil((i-1)/2) +120;
    cell(i,3) = 1 ;
end