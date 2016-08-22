clear
load('DREAM.mat');





ppi = zeros(ptn_num,ptn_num);
w = testw;
space = max(w, ppi);
link_num = length(find(space>0));
nvar = link_num*3;
vec_intvar = 1:nvar;

lb=zeros(nvar,1);
ub=ones(nvar,1);

lg_range = ceil(space)*diag(sum(ceil(space)));
lb(link_num*2+1:link_num*3)=1;
ub(link_num*2+1:link_num*3)=lg_range(space>0);

%lb(1:link_num)=2;
%ub(1:link_num)=8;

  
glb.ptn_num = ptn_num;
glb.expt_num = expt_num;
glb.obv = obv;    
glb.cyt = cyt;
glb.ihb = ihb; 
glb.perturb = perturb;
glb.v = v;   
glb.w = w;
glb.maxindg = 2;        
glb.space = space;
%glb.print_flag = 0;

% x0=zeros(nvar,1);
% x0(1:ptn*ptn)=reshape(testw,ptn*ptn,1);
% %x0(ptn*ptn+1:ptn*ptn*2)=1;%reshape(test_wrt,ptn*ptn,1);
% tmp = zeros(ptn,ptn);
% for i=1:ptn    
%     tmp(testw(:,i)==1,i)=1:sum(testw(:,i));
% end    
% 
% x0(ptn*ptn+1:nvar)=reshape(tmp,ptn*ptn,1);

glb.lambda = 1.2;
glb.gamma = 1.2;


COUNTER = 50;

options = gaoptimset('UseParallel',true);
s = RandStream('mt19937ar','Seed','shuffle');
RandStream.setGlobalStream(s);

result={};
for c=1:COUNTER
	debug=0;
    fitness = @(x)objfun_parallel(x,glb);	    
    [objx,objf,exitflag,output] = ga(fitness,nvar,[],[],[],[],lb,ub,[],vec_intvar,options);
    result = [result; {objf,objx,exitflag,output}];
    %objfun_testdream(objx);
end  



