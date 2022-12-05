clc
clear all
clear

xxxx = load ('dim.mat');

yyyy = cell2mat(struct2cell(xxxx));
energyy = yyyy(:,:,1);

timee = yyyy(:,:,2);
sinkss = yyyy(:,:,3);

mean(energyy,'ALL')


