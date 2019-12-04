%%%%Program Variables:
 %%%%input: given image of size row,col,bands and parameter k for number of
 %%%%smooth curves to cover a cluster
 %%%%b(i,k): boolean variables for each curve : dim: row*col*k, 
 %%%const: b(i,k)^2-b(i,k)=0
 %%%%____________
 %%%%wi: boolean variable for membership of each point in cluster,
 %%%%dim:row*col
 %%%%const: w(i)^2-w(i)=0
 %%%% const: sumk b(i,k)=wi
 %%%in obj: max sum w(i) 
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
 %%%%implicit var X_raise: vars for deg2 monomials of x in cluster
 %%%%dim: bands+bands+choose(bands,2)
 %%%______________
 %%%%b_imp_k: indicator membership for implicit variable in curve k
 %%%%const: b_imp_k^2-b_imp_k=0, sum_k(b_imp_k)=1
 %%%%const: -eps<b_imp_k.Ck.x_imp<eps, for all curves Ck
 %%%%
 %%%%%%_____________
 %%%%Lk: distance var of implicit X_raise from Ck
 %%%%const: -lk<Ck.X_raise<lk
 %%%%%%______________
 %%%%%for all n variables of X implicit
 %%%%-(M(li)+M(lj))<Grad_n(Ci.X_raise)^2-Grad_n(Cj.X_raise)^2<M(li)+M(lj) for large M
 %%%%second form: omit the gradient constraints and replace with 
 %%%%-(l(i)+l(j))<(Ci(X+vecEps))^2-(Cj(X+vecEps))^2<(l(i)+l(j))
 
 
 %%%%
 %%%%overall obj: minimize: -Sumi(wi)+eps
 %%% program dimensions: variable number:
 %%% row*col*(k)+1+(2*bands+Choose(bands,2))*(k+1)+k*2+mon_num
 %%%%var order: b(1,1),b(1,2),...,b(1,k),b(2,1)...
 %%%%(dataPointNum*k)||C1,...,Ck(mon_num*k)||eps||X_raise (mon_num)||b_imp_k
 %%%%(k) ||Lk (k)||vecEps (mon_num)
 %%%%b_imp_k binary constraint
varsBiSum=[];
for i=1:k
     ineqPolySys{pcount}.typeCone = -1;
        ineqPolySys{pcount}.dimVar =dimVar;
        ineqPolySys{pcount}.degree =1;
        ineqPolySys{pcount}.noTerms = mon_num+1;
        vars16=zeros(1,dimVar);
sec_degree_Coeffs;

ineqPolySys{pcount}.supports = [varst7;vars16];

coefs=ones(1+mon_num,1);
%coefs(mon_num+1:(mon_num)*(mon_num-1)/2+mon_num)=2;
coefs(1+mon_num)=-1;
        ineqPolySys{pcount}.coef= coefs;
        
        pcount=pcount+1

 
end
