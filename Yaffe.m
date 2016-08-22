clear;
%global  node_num sample_num obv ihb cyt gamma v n w w_rt leaf lambda perturb space maxPassStep debug;


load('cell_fate.mat');
%Yaffe=importdata('Yaffe.csv');
%[prior.source,prior.interact,prior.target,prior.score]=textread('cell_fate_prior.txt','%s %s %s','headerlines',1');
[measure_node,norm_data] = normalize(Yaffe,'BT20','T-D',1.5,2);
nodes = sort(unique([prior.source prior.target]));
stimulus = {'TNFR','EGFR','DNA_Damage'};
inhibitors = {'EGFR'};
leaf_nodes = {'4EBP1', 'S6', 'CYCLIN', 'Casp3', 'MYC', 'CDK'};
node_num= length(nodes);
sample_num = 7;
cyt = [];
ihb = [];


perturb = zeros(sample_num,length(nodes));

% for i=1:length(leaf_nodes)
%     leaf(i) = strmatch(leaf_nodes(i),nodes,'exact');
% end

for i=1:length(stimulus)
    cyt(i) = strmatch(stimulus(i),nodes,'exact');
    perturb(:,cyt(i)) = ones(sample_num,1);    
end

for i=1:length(inhibitors)
    ihb(i) = strmatch(inhibitors(i),nodes,'exact');
    perturb(:,ihb(i)) = ones(sample_num,1);    
end


prior_net = zeros(length(nodes),length(nodes));
for i=1:length(prior.interact)
    si = strmatch(prior.source(i),nodes,'exact');
    ti = strmatch(prior.target(i),nodes,'exact');

    if strcmp(prior.interact(i),'+')

        prior_net(si,ti) = 1; %prior.score(i);
    else        
        prior_net(si,ti) = -1; %prior.score(i);
    end
end

obv=[];
v=zeros(sample_num,length(nodes));
initial = ones(1,node_num)*0.5;
for i=1:size(map_n2mp,1)
    node_name = map_n2mp{i,1};
    measured_protein = map_n2mp{i,2};
    %disp(node_name);
    obv(i) = strmatch(node_name, upper(nodes),'exact');
    vi = strmatch(measured_protein, measure_node,'exact');
    v(:,obv(i)) = norm_data(2:end,vi);
    initial(obv(i)) = 1-norm_data(end,vi);
end

initial(strmatch('DNA_Damage',nodes,'exact'))=0.9;

n=3;
w = abs(prior_net);%zeros(length(nodes));
w_rt = prior_net; %ceil(prior_net)+floor(prior_net);
ci = zeros(node_num,node_num);
space = max(abs(prior_net),ci);
link_num = length(find(space>0));
%cytscalar_num = sample_num*length(cyt);
nvar = link_num*3;%+cytscalar_num;
vec_intvar = 1:nvar;

lb=zeros(nvar,1);
ub=ones(nvar,1);
lg_range = ceil(space)*diag(sum(ceil(space)));
lb(link_num*2+1:link_num*3) = 1;
ub(link_num*2+1:link_num*3) = lg_range(space>0);

%lb(link_num*2+1:nvar) = 0.2;
%ub(link_num*2+1:nvar) = 0.8;



% build connectivity constraints
% middle_nodes = setdiff(1:node_num,[cyt leaf]);
% A=zeros(2*length(middle_nodes),nvar);
% B=zeros(2*length(middle_nodes),1);
% con_counter = 1;  
% tmpb = zeros(node_num,node_num);
% tmpb(find(space>0))=1:link_num;
% for i=1:length(middle_nodes)
% 
%     A(con_counter,tmpb(space(:,middle_nodes(i))>0,middle_nodes(i))) = -1;
%     B(con_counter) = -1;
%     con_counter = con_counter + 1;
% 
%     A(con_counter,tmpb(middle_nodes(i),space(middle_nodes(i),:)>0)) = -1;
%     B(con_counter) = -1;
%     con_counter = con_counter + 1;
% end
    
gamma = 0.1
beta = 0.1
global_const.node_num = node_num;
global_const.sample_num = sample_num;
global_const.obv = obv;
global_const.ihb = ihb;
global_const.cyt = cyt;
global_const.gamma = gamma;
global_const.beta = beta;
global_const.v = v;
global_const.w = w;
global_const.max_indegree = 3;
global_const.nodes =nodes;
global_const.w_rt = w_rt;
global_const.space = space;
global_const.perturb = perturb;   
global_const.initial = initial;
global_const.print_flag = 0;
s = RandStream('mt19937ar','Seed','shuffle');
RandStream.setGlobalStream(s);
result={};

    opt = gaoptimset('UseParallel',true);
    fitness = @(x)objfun_parallel_yaffe_ts(x,global_const);
for i=1:50
    i
    debug=0;
    %options = gaoptimset('PopulationSize',200,'Generations',200,'StallGenLimit',50,'UseParallel','always');
    %opt = gaoptimset('plotfcns',{@gaplotbestf},'plotinterval',1);

    [objx,objf,exitflag,output] = ga(fitness,nvar,[],[],[],[],lb,ub,[],vec_intvar,opt);
    %debug=1;
    %objfun_v1(objx);    
    result = [result; {objf,objx, exitflag,output}];
end
[val, idx] = min(cell2mat(result(:,1)))
global_const.print_flag = 1;
objfun_parallel_yaffe_ts(result{idx,2},global_const);