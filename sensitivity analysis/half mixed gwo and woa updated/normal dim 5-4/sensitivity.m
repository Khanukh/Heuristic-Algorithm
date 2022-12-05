clc
clear 


xxxx = load ('dim-5-4.mat');

yyyy = cell2mat(struct2cell(xxxx));
energyy = yyyy(:,:,1);

timee = yyyy(:,:,2);
sinkss = yyyy(:,:,3);



