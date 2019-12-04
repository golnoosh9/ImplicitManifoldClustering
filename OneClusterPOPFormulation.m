%%%%Program Variables:
%%%%input: given image of size row,col,bands and parameter k for number of
%%%%smooth curves to cover a cluster
%%%%b(i,k): boolean variables for each curve : dim: row*col*k,
%%%const: b(i,k)^2-b(i,k)=0
%%%%____________

%%%________________
%%%%
%%%%
%%%%eps: curve fit error: in obj: minimize eps
%%%%%%__________
%%%%%Ck: coefficient of kth curve (deg2), dim: bands+bands+choose(bands,2)
%%%%%const: for each pixel i: -eps<b(i,k).Ck.x_raisei<eps, for all curves
%%%%%Ck
%%%%x_raisei: degree 1 and 2 monomials of pixel xi evaluation
%%%%%
%%%%_______________
%%%%implicit var X_raise_k: vars for deg2 monomials of x in cluster
%%%%dim: bands+bands+choose(bands,2)
%%%______________
%%%%b_imp_k: indicator membership for implicit variable in curve k
%%%%const: b_imp_k^2-b_imp_k=0, sum_k(b_imp_k)=1
%%%%const: -eps<b_imp_k.Ck.x_imp<eps, for all curves Ck
%%%%
%%%Add constraint for connecting average of k subclusters
 %%%to implicit variable X_raise conditioned on the corresponding 
 %%b_imp being 1
%%%%%%_____________
%%%%Lk: distance var of implicit X_raise from Ck
%%%%const: -lk<Ck.X_raise<lk
%%%%%%______________
%%%%%for all n variables of X implicit
%%%%-(M(li)+M(lj))<Grad_n(Ci.X_raise)^2-Grad_n(Cj.X_raise)^2<M(li)+M(lj) for large M
%%%%second form: omit the gradient constraints and replace with
%%%%-(l(i)+l(j))^2<(Ci(X+vecEps))^2-(Cj(X+vecEps))^2<(l(i)+l(j))^2


%%%%
%%%%overall obj: minimize: -Sumi(bi)+eps
%%% program dimensions: variable number:
%%% row*col*(k)+1+(2*bands+Choose(bands,2))*(k+1)+k*2
%%%%var order: b(1,1),b(1,2),...,b(1,k),b(2,1)...
%%%%(dataPointNum*k)||C1,...,Ck(mon_num*k)||eps||X_raise_1,...,X_raise_k (mon_num*k)
%%%% ||Lk (k)||vecEps (mon_num)
clear all;
knn=3;
k=2;
numPoints=50;
WeightIntegrityTerms=[];
WeightIntegrityCoeffs=[];
load('SalinasA.mat');
% image=salinasA_corrected;
% [rows cols bands]=size(image);
ReadDataClusterFile;
[dataPointNum bands]=size(data);
% for i=1:bands
%     meanImage=min(min(image(:,:,i)));
%     image(:,:,i)=image(:,:,i)-meanImage*ones(size(image(:,:,i)));
%     maxImage=max(max(image(:,:,i)));
%     image(:,:,i)=image(:,:,i)./maxImage*255;
% end

% image=image(1:5,1:5,1:2);
% [rows cols bands]=size(image);


for i=1:bands
    meanImage=mean(data(:,i));
    data(:,i)=data(:,i)-meanImage*ones(size(data(:,i)));
    maxImage=max(data(:,i));
    data(:,i)=data(:,i)./maxImage;
end

data=data(1:100,:);
sdata=size(data);
while(sdata(1)>numPoints)
data(randi([1,sdata(1)],5,1),:)=[];
sdata=size(data);
end
neighbors = knnsearch(data,data,'K',knn);
[~,idx] = unique(sort(neighbors,2),'rows','stable');
neighbors = neighbors(idx,:);

plot(data(:,1),data(:,2));

[dataPointNum bands]=size(data);

%data=reshape(image,[rows*cols bands]);
%dataPointNum=rows*cols;

mon_num=2*bands+nchoosek(bands,2)+1;
dimVar=dataPointNum*(k)+mon_num*k+1+mon_num*k+k+mon_num;

% syms x1 x2 x3 x4 x5 x6 x7 x8 x9 x10
% VEC=monomials([x1;x2;x3;x4;x5;x6;x7;x8;x9;x10],[0 1 2]);

syms x1 x2 
VEC=monomials([x1;x2],[0 1 2]);

pcount=1;

vars_alldatapoint_weights=[];
AddKNNConstraints;
for j=1:dataPointNum
    
%     eval=(subs(VEC,[x1,x2,x3,x4,x5,x6,x7,x8,x9,x10],[data(j,1),data(j,2),data(j,3),...
%         data(j,4),data(j,5),data(j,6),data(j,7),data(j,8),data(j,9),data(j,10)]));
    
     eval=(subs(VEC,[x1,x2],[data(j,1),data(j,2)]));
    evald=double(eval);
     AddDataWeightsConstraint;
     AddDataPointCurveEvaluationConstraint
end
 AddCoeffNormalizationConstraint;
   AddImplicitEvaluationInCurveConstraint;
   AddImplicitEvaluationBoundConstraint;
  %to be modified%AddSmoothNessConstraint;
  
  coeffs_datapoint_and_data_weights=[];
vars_datapoint_and_data_weights=[];
 vars_implicit_and_data_weights=[];
 AddImplicitVarDataAverageConstraint;
 
s_knn_objective_term=size(knn_cluster_support);
% gen_weight;
objPoly.typeCone = 1;
objPoly.dimVar   = dimVar;
objPoly.degree   = 2;
objPoly.noTerms  = 1+dataPointNum*k+(1+mon_num)*k*dataPointNum;%sum of weights, inv epsilon, curve fi epsilon, 2 curve insertion
vars_eps=zeros(1,dimVar);
vars_eps(dataPointNum*k+mon_num*k+1)=1;%%%for e

objPoly.supports =[vars_eps;...
    vars_alldatapoint_weights;...
    vars_datapoint_and_data_weights;vars_implicit_and_data_weights];

%  total_coeffs=zeros(size(total_coeffs));



objPoly.coef     = [1;-zeros(dataPointNum*k,1)/dataPointNum*k;...
    -coeffs_datapoint_and_data_weights;ones(mon_num*k*dataPointNum,1)];
%%%%var order: b(1,1),b(1,2),...,b(1,k),b(2,1)...
%%%%(dataPointNum*k)||C1,...,Ck(mon_num*k)||eps||X_raise (mon_num)||b_imp_k
%%%%(k) ||Lk (k)||vecEps (mon_num)

lbd=-100*zeros(1,dimVar); 
lbd(dataPointNum*k+1:dataPointNum*k+mon_num*k)=-10;%%coeffs
lbd(1:dataPointNum*k)=-500.9;
lbd(dataPointNum*k+mon_num*k+1+mon_num*k+k+1:dataPointNum*k+mon_num*k+1+mon_num*k+k+mon_num)=0.01;
lbd(dataPointNum*k+mon_num*k+1+1:dataPointNum*k+mon_num*k+1+1+mon_num*k)=-1;%%X_impl (normaliz)


ubd=ones(1,dimVar); 
ubd(1:dataPointNum*k)=500.9;
ubd(dataPointNum*k+1:dataPointNum*k+mon_num*k)=10;%%coeffs
ubd(dataPointNum*k+mon_num*k+1+mon_num*k+k+1:dataPointNum*k+mon_num*k+1+mon_num*k+k+mon_num)=0.1;

dataNum=dataPointNum;
AssignToZero;

dataNum=1;
AssignToOne;

param.scalingSW=0;
%param.SDPsolverEpsilon=0.000000000001;
param.relaxOrder =2;

%param.solver='sdpa'

[param,SDPobjValue,POP,cpuTime,SDPsolverInfo,SDPinfo]=sparsePOP(objPoly,ineqPolySys,lbd,ubd,param);

   solution=SDPinfo.y;
   coeffs_var=solution(cnum+mon_num+2+1:cnum+mon_num+2+mon_num+mon_num*(mon_num-1)/2);
   coeffs_covariance=zeros(mon_num);
   
   count=1;
   for i=1:mon_num
       for j=i:mon_num
           coeffs_covariance(i,j)=coeffs_var(count);
           coeffs_covariance(j,i)=coeffs_var(count);
           count=count+1;
       end
   end
   [V,D] = eig(coeffs_covariance);
   sample_coef=mgd(1,mon_num,coeffs,coeffs_covariance);
   plot_implicit_six(sample_coef);




