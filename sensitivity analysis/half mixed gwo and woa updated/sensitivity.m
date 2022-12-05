clc
clear all

xxxx = load ('dim.mat');

y = cell2mat(struct2cell(xxxx));
energy = y(:,:,1);
energy = energy ;
time = y(:,:,2);
sinks = y(:,:,3);

mean(energy,'ALL')


