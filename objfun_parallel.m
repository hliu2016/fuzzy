%*******************************************************%
%***** Fuzzy logic modeling signaling pathway **********%
%*****           @author: liuhui              **********%
%*******************************************************%

% Objective function for paralleled Genetic Algorithm 
function [f] = objfun_parallel(x, global_const)    
    %fetch global setting from input argument
    ptn_num=global_const.ptn_num;
    expt_num=global_const.expt_num;
    obv = global_const.obv;
    ihb = global_const.ihb;
    cyt = global_const.cyt;
    lambda = global_const.lambda;
    gamma = global_const.gamma;
    v = global_const.v;
    w = global_const.w;
    dmax = global_const.maxindg;
    %w_rt = global_const.w_rt;
    structure_space = global_const.space;
    perturb = global_const.perturb;
               
    b=zeros(ptn_num,ptn_num);
    b_rt= zeros(ptn_num,ptn_num);
    b_lg= zeros(ptn_num,ptn_num);    
    links = find(structure_space>0);                % space= w + ci_edges
    link_num = length(links);
    
    % parse variable being optimized to network topology and logic gates
    b(links) = x(1:link_num);                
    b_rt(links) = x(link_num+1:2*link_num);   
    b_lg(links) = x(link_num*2+1:3*link_num);     
    %p(links) = x(2*link_num+1:3*link_num);
    p = 0.5;
    n = 2;
        
    
    % initialze constant parameters
    burnin = ptn_num;
    interval = ptn_num;
    
    maxPassStep=ptn_num;
    threshold = 1e-6;
    y = ones(expt_num, ptn_num)*0.5;   % initialize node states to 0.5
    % simulating each sample based on current network structure
    for k = 1:expt_num
        tmpy=zeros(ptn_num, ptn_num);
        simu_buffer = zeros(maxPassStep,ptn_num);        
        %initalize input node with respect to each sample
        y(k,cyt)=0;     
        perturb_cyt = intersect(find(perturb(k,:)==1),cyt);
        y(k, perturb_cyt) = 1;
        obv_cyt = intersect(perturb_cyt,obv);
        if ~isempty(obv_cyt)
            y(k, obv_cyt) = v(k,obv_cyt);
        end
        % initialize inhibited node with respect to each sample
        perturb_ihb = intersect(find(perturb(k,:)==1),ihb);
        y(k, perturb_ihb) = 0;

        for loop=1:maxPassStep  % simulate the network to maximum step
            for j = 1:ptn_num
                if ~ismember(j,cyt) && ~(ismember(j,ihb) && perturb(k,j)==1)
                    for i = 1:ptn_num
                        if b(i,j)>0
                            if b_rt(i,j)>0
                                tmpy(j,i) = (p^n+1)*y(k,i)^n/(y(k,i)^n+p^n);
                            elseif b_rt(i,j)<=0
                                tmpy(j,i) = 1-(p^n+1)*y(k,i)^n/(y(k,i)^n+p^n);
                            end
                        end                     
                    end
                    andgate = 0;
                    orgate = 0;  
                    
                    inhibit_idx = intersect(find(b_rt(:,j)<=0),find(b(:,j)>0));
                    activat_idx = intersect(find(b_rt(:,j)>0),find(b(:,j)>0));
                    % parse the b_lg to logical relationship
                    for lg = 1:max(b_lg(:,j)) 
                        % coupled inhibitor with activitor 
                        atv_effect = intersect(activat_idx,find(b_lg(:,j)==lg));   
                        ihb_effect = intersect(inhibit_idx,find(b_lg(:,j)==lg)); 
                        if ~isempty(atv_effect)                            
                            andgate = min(tmpy(j,atv_effect));          
                            if ~isempty(ihb_effect)
                                andgate = min(andgate, tmpy(j,ihb_effect));
                                inhibit_idx = setdiff(inhibit_idx, ihb_effect);
                            end
                        end
                        orgate = max(orgate,andgate);
                    end
                    % outer layer inhibitors stand for strong inhibition 
                    if ~isempty(inhibit_idx)                              
                        y(k,j) = min(min(tmpy(j,inhibit_idx)),orgate);
                    else
                        y(k,j) = orgate;
                    end
                end
            end            
            
            simu_buffer(loop,:) = y(k,:);
            
          % check whether converge to steady state            
            if loop>=burnin
                latest_ts = simu_buffer(loop-interval+1:loop,:);
                diff = latest_ts-repmat(mean(latest_ts),size(latest_ts,1),1);
                err = sum(sum(diff.^2));
                if err<threshold
                    %disp(loop);
                    break;
                end
            end    
            % call unit root test whether oscillation steady state is reached   
%             if loop>=burnin+interval && mod(loop-burnin,interval)==0
%             	flag = true;
%                 for p=1:ptn_num
%                 	ts = simu_buffer(loop-burnin:loop,p);
%                     ur = adf(ts,0,2);
%                     if ur.adf>ur.crit(3)
%                     	flag = false;
%                     end
%                 end
%                 if flag==true   
%                     disp(loop);
%                 	break;
%                 end
%             end            
        end               
    end
    
    % compute the objective value as square error sum plus penalty term
	lse = 0;
    penalty = 0;

    err = y(:,obv)-v(:,obv);
    lse = sum(sum(err.^2));
    penalty = lambda*sum(sum((1-2*w).*b))+gamma*sum(max(0,sum(b)-dmax));
    f = lse+penalty;
    disp(lse);
end


