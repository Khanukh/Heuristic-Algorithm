clc
clear all 
clear
global ssensors;
global ssink_types;
global  ssinks;
global wwsn;
global YY;
global Y_tmp_sensor_total_path_to_sinks;
global Y_Total_communication_energy;
ssensors = xlsread('sensors.csv');
ssink_types = xlsread('sink_types.csv');
ssinks= xlsread('sinks.csv');
wwsn = xlsread('graph.csv');
 
 
% sensor data   vertex name kc maxl cpu ram bw
Y_sensors = ssensors;
% sink_type data   cpu ram bw cost
Y_sink_types = ssink_types;
 
% sinks   vertex x y 
Y_sinks= ssinks;
Y_sink_vertex = min(Y_sinks)+1;
YY = zeros(Y_sink_vertex-1, length(ssinks) );
 
 
 
 
 
 
 
 
% wsn
Y_wsn =wwsn;
Y_wsn(:,1:2)= Y_wsn(:,1:2) +1;
Y_nodes = max( max(Y_wsn(:,1)) ,max(Y_wsn(:,2)));
% creating names
 
 
%% creating Y
Y_G = graph(Y_wsn(:,1),Y_wsn(:,2),Y_wsn(:,3)) ; 
 
 
for i = 1 :  Y_sink_vertex -1
    for j = 1 : Y_nodes - Y_sink_vertex + 1  
            if (   (length( shortestpath(Y_G,i,j + Y_sink_vertex-1) )  - 1 )  <= Y_sensors(i,4) )
                    YY(i,j) = 1 ;
                   
            end
    end
end
 
%% Creating H
tmp_sensor_total_path_to_sinks = zeros(length(Y_sensors),length(Y_sinks));
for i = 1 : length(Y_sensors)
  
   
    for j = 1 : length(Y_sinks)
        if YY(i,j) == 1 
            tmp_sensor_total_path_to_sinks(i,j) =   length (shortestpath(Y_G, i, j + Y_sink_vertex-1)) -2  ;                     
        end       
    end
 
end
Y_tmp_sensor_total_path_to_sinks = tmp_sensor_total_path_to_sinks;
 
HH = sum(Y_tmp_sensor_total_path_to_sinks,2)';
%% creating communication energy
ps_communication = zeros(1,length(Y_sinks)) ;
 
E_elec = 10; % nj/Bit
K_E_elec = 2;% in Kbit
K_total = length(Y_sensors) * 2;
 
Y_Total_communication_energy= zeros(length(Y_sensors),length(Y_sinks));
for i = 1 : length(Y_sensors)
%       display ("next sensor "  +  i )
  
    for j = 1 : length(Y_sinks)
        if YY(i,j) == 1 
            if (HH(i) == 0 )
                continue;
            end
            Y_Total_communication_energy(i,j) =  (K_E_elec*  E_elec ) * exp( (length (shortestpath(Y_G, i, j + Y_sink_vertex-1)) -2 ) /HH(i));
         
           
        end       
    end
    
 
%     display ("end sensor "  +  i + " total energy of " + Total_communication_energy )
end
 
 


