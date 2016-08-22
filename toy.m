clear

load('ts_synthetic.mat');

ppi = zeros(node_num,node_num);
  

glb.truew = truew;    
glb.node_num = node_num;

glb.obv = obv;    
glb.cyt = cyt;

glb.ts = ts(1:10);

  
opt = gaoptimset('UseParallel','always');
s = RandStream('mt19937ar','Seed','shuffle');
RandStream.setGlobalStream(s);

for g = 1:10    
    gamma = 0.01
    beta = gamma;   
        prior_percent=1; % 50 percent prior edge relative to total true edges
        w=zeros(node_num,node_num);
        edges=find(truew==1);
        for i=1:prior_percent*sum(sum(truew>0))
            index = randi(length(edges));
            while w(edges(index))==1
                index = randi(length(edges));
            end
            w(edges(index))=rand(1,1);
        end

        % introduce noisy edges into the prior information    
        noisy_percent = 0.5; 
        noisy_edge_num = round(noisy_percent*sum(sum(truew>0)));
        for i=1:noisy_edge_num 
            row = randi(node_num);
            col = randi(node_num);
            while ismember(col,cyt) || col==row || w(row,col)>0
                col = randi(node_num);
            end        
            w(row,col) = rand(1,1);
        end

        space = max(w, truew);
        link_num = length(find(space>0));
        nvar = link_num*3;
        vec_intvar = 1:link_num*3;

        lb=zeros(nvar,1);
        ub=ones(nvar,1);

        lg_range = ceil(space)*diag(sum(ceil(space)));
        lb(link_num*2+1:link_num*3)=1;
        ub(link_num*2+1:link_num*3)=lg_range(space>0);
        
        glb.gamma = gamma;
        glb.beta = beta;
        glb.w = w;
        glb.maxindg = 2;        
        glb.space = space;
        glb.print_flag = 0;
        
        result = {};
        for c=1:50               
            fitness = @(x)objfun_parallel_ts(x,glb);
            [objx,objf,exitflag,output] = ga(fitness,nvar,[],[],[],[],lb,ub,[],vec_intvar);
            result = [result; {objf,objx, exitflag,output}];
            
            b=zeros(node_num,node_num);
            links = find(space>0);                % space= w + ci_edges
            link_num = length(links);
            b(links) = objx(1:link_num);   
            penalty = gamma*sum(sum((1-2*w).*b))+beta*sum(max(0,sum(b)- glb.maxindg));
            lse = objf - penalty;
            sd = sum(sum(abs(b-truew)));
            fprintf('%f\t%d\t%d',lse, sd, sum(sum(b)));
        end
        [objf,index] = min(cell2mat(result(:,1)));
        %x = result{index,2};
               
             


end                


[val, idx] = min(cell2mat(result(:,1)))
glb.print_flag = 1;
objfun_parallel_ts(result{idx,2},glb);





