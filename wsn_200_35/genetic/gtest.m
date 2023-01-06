global ssensors;
global ssink_types;
global  ssinks;
global wwsn;
% sensor data   vertex name kc maxl cpu ram bw
sensors = ssensors;

% sink_type data   cpu ram bw cost
sink_types = ssink_types;


% sinks   vertex x y 
sinks= ssinks;
sink_vertex = min(sinks)+1;


x = [ 1 1 4 4]
sink_initial_type =x;


% wsn
wsn = wwsn;
wsn(:,1:2)= wsn(:,1:2) +1;
nodes = max( max(wsn(:,1)) ,max(wsn(:,2)));
% creating names
str_names = string; 
    for i = 1 : nodes
        if(i <sink_vertex ) 
            str_names(i) = "sensor" + i ;
        end
        if(i >=sink_vertex ) 
            str_names(i) = "sink" + (i - sink_vertex +1 ) ;
        end
    end
G = graph(wsn(:,1),wsn(:,2),wsn(:,3),str_names) ; 



%% Y matrix 


% y matrix shortest path for each sink 
Y = zeros(sink_vertex-1, length(sinks) );
for i = 1 :  sink_vertex -1
    for j = 1 : nodes - sink_vertex + 1  
            if (   (length( shortestpath(G,i,j + sink_vertex-1) )  - 1 )  <= sensors(i,4) )
                    Y(i,j) = 1 ;
                    if sink_initial_type(j) == 1
                        Y(i,j) = 0;
                    end
            end
    end
end
%% p1 reliability constraint
p1 = 0 ; 
for i = 1 : sink_vertex -1
    p1 = p1 +  max(0 , sensors(i,3) - sum(Y(i,:) ) );
end
%% p2 max cpu workload  for each sink sum of sensors connected to sink's workload / number of sinks sensor is connected 2
p2 = 0 ; 
sink_sum_load = 0;
for i = 1 : nodes - sink_vertex + 1 
    for j =1 : sink_vertex -1
        if (Y(j,i) == 1)
            sink_sum_load = sink_sum_load + sensors(j,5) / sum(Y(j,:));
        end
        
    end
    
    p2 = p2 +  max(0 , sink_sum_load - sink_types(sink_initial_type(i) ,1) ) ;
    p2_data(i,1) = sink_sum_load;
    p2_data(i,2) = sink_types(sink_initial_type(i) ,1);
    sink_sum_load = 0;
end
%% p3 max bandwith workload
p3 = 0 ; 
sink_sum_load = 0;
for i = 1 : nodes - sink_vertex + 1 
    for j =1 : sink_vertex -1
        if (Y(j,i) == 1 )
            sink_sum_load = sink_sum_load + sensors(j,6) / sum(Y(j,:));
        end
        
    end
    
    p3 = p3 +  max(0 , sink_sum_load - sink_types(sink_initial_type(i) ,2) ) ;
    p3_data(i,1) = sink_sum_load;
    p3_data(i,2) = sink_types(sink_initial_type(i) ,2);
    sink_sum_load = 0;
end
%% p4 max memory workload

p4 = 0 ; 
sink_sum_load = 0;
for i = 1 : nodes - sink_vertex + 1 
    for j =1 : sink_vertex -1
        if (Y(j,i) == 1 )
            sink_sum_load = sink_sum_load + sensors(j,7) / sum(Y(j,:));
        end
        
    end
    
    p4 = p4 +  max(0 , sink_sum_load - sink_types(sink_initial_type(i) ,3) ) ;
    p4_data(i,1) = sink_sum_load;
    p4_data(i,2) = sink_types(sink_initial_type(i) ,3);
    sink_sum_load = 0;
end
%% ps power consumption
ps_communication = zeros(1,length(sinks)) ;

E_elec = 10; % nj/Bit
K_E_elec = 2;% in Kbit
K_total = length(sensors) * 2;
Total_communication_energy = 0;
%%%%% clalculating each sensor max path
for i = 1 : length(sensors)
  
    sensor_total_path_to_sinks = 0;
    for j = 1 : length(sinks)
        if Y(i,j) == 1 
            sensor_total_path_to_sinks = sensor_total_path_to_sinks +   length (shortestpath(G, i, j + sink_vertex-1)) -2  ;                     
        end       
    end
    H(i) = sensor_total_path_to_sinks;
end

%%%%% clalculating each sensor energy
for i = 1 : length(sensors)
%       display ("next sensor "  +  i )
   Total_communication_energy = 0;
    for j = 1 : length(sinks)
        if Y(i,j) == 1 
            if (H(i) == 0 )
                continue;
            end
            Total_communication_energy = Total_communication_energy +    (K_E_elec*  E_elec ) * exp( (length (shortestpath(G, i, j + sink_vertex-1)) -2 ) /H(i));
%              display (" sensor "  +  i  + " to node "  + (j + sink_vertex-1)  +  " = sink "  + j + " energy " +  ( (K_E_elec*  E_elec ) * exp( (length (shortestpath(G, i, j + sink_vertex-1)) -2 ) /H(i))) )
%              display( " with min distance of "  + ( length (shortestpath(G, i, j + sink_vertex-1)) -2 )  +   " considering it's total H is " + H(i) ) 
           
           
        end       
    end
    ps_communication(i) = Total_communication_energy;
%     display ("end sensor "  +  i + " total energy of " + Total_communication_energy )
end

%%%%%%%%%%%%%%%%%%%%%%%    process      %%%%%%%%%%%%%%%%%%%%%%%%%%%%
ps_process = zeros(1,length(sinks)) ;

sink_sum_load = 0;
for i = 1 : nodes - sink_vertex + 1 
    p_idle = sink_types(sink_initial_type(i) ,6) * sink_types( sink_initial_type(i) ,5 );
    p_residual = sink_types(sink_initial_type(i) ,6) - p_idle ;
    for j =1 : sink_vertex -1
        if (Y(j,i) == 1)
            sink_sum_load = sink_sum_load + sensors(j,5) / sum(Y(j,:));
     
        end
        
    end
        if(sink_types(sink_initial_type(i) ,1) ~= 0)
            
        ps_process(i) = sink_sum_load/sink_types(sink_initial_type(i) ,1)*p_residual + p_idle;
        else
            ps_process(i) = 0  ;
        end
       sink_sum_load = 0;     
 
end


%%%%%%%%%%%%%%%%%%%%%%%    end      %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    PS = sum(ps_process) + sum(ps_communication)/(E_elec + K_E_elec); % dimension lessed of energy

%% cost 
cost = 0 ;    
cost = sum(sink_types(sink_initial_type,4));


%% Fitness Function 
alpha = 1 ;
beta = 1000 ;
fobj = alpha*(PS   )  + beta*(p1 + p2 + p3 + p4 )  