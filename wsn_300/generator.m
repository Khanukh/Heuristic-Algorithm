data = readtable('converter.csv');
data = table2array(data);
for i = 1 : size(data,1)
    data(i,1) = data(i,1) + 120;
    if data(i,2) >59
        data(i,2) = data(i,2) + 200;
    else
        data(i,2) = data(i,2)+120;
    end
end