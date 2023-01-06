clc
clear all 
clear
global ssensors;
global ssink_types;
global  ssinks;
global wwsn;
ssensors = xlsread('sensors.csv');
ssink_types = xlsread('sink_types.csv');
ssinks= xlsread('sinks.csv');
wwsn = xlsread('graph.csv');