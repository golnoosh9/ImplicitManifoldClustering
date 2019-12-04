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
 %%%%overall obj: minimize: -Sumi(wi)+eps+var_knn_evaluationcurvek
 %%% program dimensions: variable number:
 %%% row*col*(k)+1+(2*bands+Choose(bands,2))*(k+1)+k*2+mon_num
 %%%%var order: b(1,1),b(1,2),...,b(1,k),b(2,1)...
 %%%%(dataPointNum*k)||C1,...,Ck(mon_num*k)||eps||X_raise (mon_num)||b_imp_k
 %%%%(k) ||Lk (k)||vecEps (mon_num)
 
 %%%% ad the lk constraints
 knn_cluster_support=[];%%order:
 rowsToDelete=[];
 %%sum(b(neigh(nn),i)^2)*(KNN-1)/KNN-2*sum(bi.bj)
 sneighbors=size(neighbors);
 for i=1:sneighbors(1)-1
     if(ismember(i,sneighbors))
         continue;
     end
     for j=1:knn
         for ii=i+1:sneighbors(1)
            if(ismember(neighbors(i,j),neighbors(ii,:)))
             rowsToDelete=[rowsToDelete;ii];
            end
         end
     end
 end
 neighbors(rowsToDelete,:)=[];
 
 sneighbors=size(neighbors);
 knn_cluster_coeffs=[];
 for i=1:k%%applies for all Ck
    for j=1:sneighbors(1)%%for every point's neighborhood
        for nn=1:knn
            vars_square=zeros(1,dimVar);
            vars_square((neighbors(j,nn)-1)*k+i)=2;
            knn_cluster_support=[knn_cluster_support;vars_square];
            knn_cluster_coeffs=[knn_cluster_coeffs;(knn-1)/knn];
            for nnn=nn+1:knn
                vars_cross=zeros(1,dimVar);
                vars_cross((neighbors(j,nn)-1)*k+i)=1;
                vars_cross((neighbors(j,nnn)-1)*k+i)=1;
                knn_cluster_support=[knn_cluster_support;vars_cross];
                knn_cluster_coeffs=[knn_cluster_coeffs;-2/knn];
            end
            
        end
    end
 end
 