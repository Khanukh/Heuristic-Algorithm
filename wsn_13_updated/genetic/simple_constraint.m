function [c, ceq] = simple_constraint(t)


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


x = t;
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

 c1 = p1;%(abs(sqrt(KKeq/M)/(2*pi)-WWn)/WWn*100)-100 ;
     % Compute sd
 c2 = p2;%(abs(M*9.81*1000/KKeq - SSd)/SSd*100)-100;
     % Compute BNEA
  c3 = p3;%(abs(sqrt(4*1.38e-23*300/(M * 1e3 *Q*(Wn*2* pi)^3)) / (Sd / 1e3) - BBNEA)/BBNEA)-700;
  % area

  c4 = p4  ;
   C = [c1;c2;c3;c4];
   
   c = sum(C);
  
   ceq = [];
   