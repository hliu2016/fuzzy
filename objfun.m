%*******************************************************%
%***** Fuzzy logic modeling signaling pathway **********%
%*****           @author: liuhui              **********%
%*******************************************************%

% define the objective function for Genetic Algorithm 
function [f,lse,b] = objfun(x)
    global  ptn_num expt_num obv ihb cyt gamma v n w space perturb lambda maxPassStep debug;           

    b=zeros(ptn_num,ptn_num);
    b_rt=ones(ptn_num,ptn_num);
    b_lg=zeros(ptn_num,ptn_num);
    p=zeros(ptn_num,ptn_num);
    links = find(space>0);                
    link_num = length(links);
    
    % decode variable x to matrices 
    b(links) = x(1:link_num);                
    b_rt(links) = x(link_num+1:2*link_num);   
    b_lg(links) = x(link_num+1:2*link_num);    
    p(links) = x(3*link_num+1:4*link_num)/10;
    
    
    y = zeros(expt_num, ptn_num);          % predicted value
    tmpy=zeros(expt_num, ptn_num, ptn_num);

    
    for k = 1:expt_num
        tmp = zeros(1,ptn_num);
        for j = 1:length(cyt)
            y(k,cyt(j)) = perturb(k,cyt(j)); 
            tmp(cyt(j)) = perturb(k,cyt(j));
        end
        for j = 1:length(ihb)
            if perturb(k,ihb(j))==1
                y(k,ihb(j)) = 1-perturb(k,ihb(j));
                tmp(ihb(j)) = 1-perturb(k,ihb(j));
            end  
        end
        for loop=1:maxPassStep 
            for j = 1:ptn_num
                if ~ismember(j,cyt) && ~(ismember(j,ihb) && perturb(k,j)==1)
                    for i = 1:ptn_num
                        if b(i,j)>0
                            if b_rt(i,j)>0
                                tmpy(k,j,i) = (p(i,j)^n+1)*y(k,i)^n/(y(k,i)^n+p(i,j)^n);
                            else
                                tmpy(k,j,i) = 1-(p(i,j)^n+1)*y(k,i)^n/(y(k,i)^n+p(i,j)^n);
                            end
                        end                     
                    end
                    andgate = 0;
                    orgate = 0;
                    for lg = 1:max(b_lg(b(:,j)>0,j)) 
                        coreg = intersect(find(b(:,j)>0),find(b_lg(:,j)==lg));                        
                        if ~isempty(coreg)                            
                            andgate = min(tmpy(k,j,coreg));                            
                        end
                        orgate = max(orgate,andgate);
                    end
                    tmp(j) = orgate;
                end
            end            
            y(k,:)=tmp;
        end
    end
    
	lse = 0;
    penalty = 0;
    
    for k = 1:expt_num    
        for i = 1:length(obv)
            lse = lse + (y(k,obv(i))-v(k,obv(i)))^2;            
        end
    end
    %lse
    for j = 1:ptn_num
        for i = 1:ptn_num
            penalty = penalty + (lambda-w(i,j))*b(i,j);
        end
    end    
    f = lse+penalty*gamma;
    
    if debug==1
        addedge = 0;
        deletedge = 0;
        lse
        f
        b     
        for j=1:ptn_num
            for i=1:ptn_num
                addedge = addedge+(1-ceil(w(i,j)))*b(i,j);
                deletedge = deletedge + (1-b(i,j))*ceil(w(i,j));
            end
        end        
        addedge
        deletedge
    end
end


