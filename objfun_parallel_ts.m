%*************************************************************************%
%****** Knowledge-guided Fuzzy logic modeling signaling pathway **********%
%**********************        @author: Liu Hui            ***************%
%*************************************************************************%

% Objective function for paralleled Genetic Algorithm 
function [f] = objfun_parallel_ts(x, glb)    

    %set variable to global setting
    node_num=glb.node_num;
    %sample_num=glb.sample_num;
    obv = glb.obv;    
    cyt = glb.cyt;
    gamma = glb.gamma;
    beta = glb.gamma;
    ts = glb.ts;
    w = glb.w;  
    print_flag = glb.print_flag;
    dmax = glb.maxindg;    
    structure_space = glb.space;      
    
    % define temporary variable
    b=zeros(node_num,node_num);
    b_rt= zeros(node_num,node_num);
    b_lg= zeros(node_num,node_num);    
    links = find(structure_space>0);                % space= w + ci_edges
    link_num = length(links);
    
    % parse variable being optimized to network topology and logic gates
    b(links) = x(1:link_num);                
    b_rt(links) = x(link_num+1:2*link_num);   
    b_lg(links) = x(link_num*2+1:3*link_num);     
    %p(links) = x(2*link_num+1:3*link_num);    
   
    
    
    % initialze constant parameters
    p = 0.5;
    n = 2;
    
    % simulation based on current network structure
    fitness = 0;  
    for t=1:length(ts)
        v = (abs(ts{t}))';
        ts_len =size(v,1);
        y = zeros(ts_len,node_num);  
        y(1,:) = v(1,:);             % setting initalize values
        
        tmpy=zeros(node_num, node_num);


        for k=2:ts_len  % simulate the network to maximum step
                y(k,cyt) = y(k-1,cyt);
                for j = 1:node_num
                    if ~ismember(j,cyt) % && ~(ismember(j,ihb) && perturb(k,j)==1)                    
                        for i = 1:node_num % simulate the values of each nodes
                            if b(i,j)>0
                                if b_rt(i,j)>0
                                    tmpy(j,i) = (p^n+1)*y(k-1,i)^n/(y(k-1,i)^n+p^n);
                                elseif b_rt(i,j)<=0
                                    tmpy(j,i) = 1-(p^n+1)*y(k-1,i)^n/(y(k-1,i)^n+p^n);
                                end
                            end                     
                        end

                        %=== parse the b_lg to logical relationship ===% 
                        andgate = 0;
                        orgate = 0;                      
                        inhibit_idx = intersect(find(b_rt(:,j)<=0),find(b(:,j)>0));
                        activat_idx = intersect(find(b_rt(:,j)>0),find(b(:,j)>0));                  
                        for lg = 1:max(b_lg(:,j)) 
                            % parse inhibitors coupled with activitor 
                            atv_effect = intersect(activat_idx,find(b_lg(:,j)==lg));   
                            ihb_effect = intersect(inhibit_idx,find(b_lg(:,j)==lg)); 
                            if ~isempty(atv_effect)             
                                andgate = min(tmpy(j,atv_effect));
                                if ~isempty(ihb_effect)
                                    andgate = min(min(andgate, tmpy(j,ihb_effect)));
                                    inhibit_idx = setdiff(inhibit_idx, ihb_effect);
                                end
                            end
                            orgate = max(orgate,andgate);
                        end
                        % parse inhibitors standing for strong inhibition 
                        if ~isempty(inhibit_idx)                                 
                            y(k,j) = min(min(tmpy(j,inhibit_idx)),orgate);
                            %fprintf('%d\t%d\t%f\t%f\n',k,j,y(k,j),orgate);
                        else
                            %fprintf('%d\t%d\t%f\t%f\n',k,j,y(k,j),orgate);
                            y(k,j) = orgate;
                        end                    
                        %=== end parse logical relationship ===% 
                    end
                end                      
        end               


        % compute fitness between simulation and samples using time series alignment               
            ts1 = y(2:ts_len,obv);
            ts2 = v(2:ts_len,obv);
            dist = cdtw(ts1,ts2,0);
            if isnan(dist)
                fitness = inf;
                disp('fitness is infinity.');
            else
                fitness = fitness+dist;
            end
        
    end
	
    % compute the penalty term 
    penalty = gamma*sum(sum((1-2*w).*b))+beta*sum(max(0,sum(b)-dmax));
    f = fitness+penalty;    

    
    if print_flag==1
        b
        b_rt
        b_lg
        %sw_index
        y(:,obv)
        sum(sum(b))
            
        for i = 1:node_num
            for j = 1:node_num
                if b(i,j)==1
                    if b_rt(i,j)==1
                        fprintf('%s\t%s\t%s\n',nodes{i},'+',nodes{j});
                    else
                        fprintf('%s\t%s\t%s\n',nodes{i},'-',nodes{j});
                    end
                end
            end
        end
    end
end





