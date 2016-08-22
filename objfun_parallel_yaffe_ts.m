%*************************************************************************%
%****** Knowledge-guided Fuzzy logic modeling signaling pathway **********%
%**********************        @author: Liu Hui            ***************%
%*************************************************************************%

% Objective function for paralleled Genetic Algorithm 
function [f] = objfun_parallel_yaffe_ts(x, global_const)    
    %set variable to global setting
    node_num=global_const.node_num;
    sample_num=global_const.sample_num;
    obv = global_const.obv;
    %ihb = global_const.ihb;
    cyt = global_const.cyt;
    gamma = global_const.gamma;
    beta = global_const.gamma;
    v = global_const.v;
    w = global_const.w;
    dmax = global_const.max_indegree;
    initial = global_const.initial;
    nodes = global_const.nodes;
    %w_rt = global_const.w_rt;
    structure_space = global_const.space;
    %perturb = global_const.perturb;   
    print_flag = global_const.print_flag;
    
    % define temporary variable
    b=zeros(node_num,node_num);
    b_rt= zeros(node_num,node_num);
    b_lg= zeros(node_num,node_num);    
    links = find(structure_space>0);                % space= w + ci_edges
    link_num = length(links);
    
    % parse variable being optimized to network     topology and logic gates
    b(links) = x(1:link_num);                
    b_rt(links) = x(link_num+1:2*link_num);   
    b_lg(links) = x(link_num*2+1:3*link_num);     
    %p(links) = x(2*link_num+1:3*link_num);    
   
    
    
    % initialze constant parameters
    burnin = node_num;
    interval = sample_num;
    
    maxPassStep=3*node_num;
    threshold = 1e-3;
    p = 0.5;
    n = 2;
    
    % simulation based on current network structure
      
        
	y = zeros(maxPassStep,node_num);  
    %y(1,:) = ones(1,node_num)*0.5;
    y(1,:)=initial;         % setting initalize values
    tmpy=zeros(node_num, node_num);
        
        
	for k=2:maxPassStep  % simulate the network to maximum step
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

            
          %=== check whether converge to steady state ===%
            if k>=burnin
                latest_ts = y(k-interval:k,:);
                diff = latest_ts-repmat(mean(latest_ts),size(latest_ts,1),1);
                residual = sum(sum(diff.^2));
                if residual<threshold
                    %disp(loop);
                    break;
                end
            end           
	end               
       

    % compute fitness between simulation and samples using time series alignment
    
    %u = y(k-interval:k,:);
    sw_size = 2*sample_num;  % sliding window size    
    sstep = floor(sample_num/2);
    fitness = inf;
    sw_index = 0;
	for j=1:sstep:k-sw_size        
       err = 0;
        for i=1:length(obv)        
            ts1 = y(j:sw_size+j,obv(i));
            ts2 = v(:,obv(i));
            dist = dtw(ts1,ts2,0);
            %disp(dist);
            if isnan(dist)
                err = inf;
            else
                err = err + dist;
            end
        end

%         if isnan(dist)
%         	err = err+inf;
%         else
%             err = err+dist;
%         end
        if err<fitness
            fitness = err;
            sw_index = j;
        end
	end
    %disp(fitness);
    %y(:,obv)
    
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





