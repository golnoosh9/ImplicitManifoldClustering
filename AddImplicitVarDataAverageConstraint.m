
%%%%
%%%%overall obj: minimize: -Sumi(bi)+eps
%%% program dimensions: variable number:
%%% row*col*(k)+1+(2*bands+Choose(bands,2))*(k+k)+k
%%%%var order: b(1,1),b(1,2),...,b(1,k),b(2,1)...
%%%%(dataPointNum*k)||C1,...,Ck(mon_num*k)||eps||X_raise_1,...,X_raise_k (mon_num*k)
%%%% ||Lk (k)||vecEps (mon_num)
%%% terms=sum(b_i_k)*X_raise_k-(sum(b_i_k*X_data_i))
data_samples=dataPointNum;
monomial_samples=mon_num;
beginloop=1;
 coeffs_datapoint_and_data_weights=[];
vars_datapoint_and_data_weights=[];
 vars_implicit_and_data_weights=[];
 coeffs_datapoint_and_data_weights=zeros(k*data_samples,1)
for i=1:monomial_samples
    coeffs_monomials=[];
    for j=1:k
        

        for d=1:data_samples
            eval=(subs(VEC,[x1,x2],[data(d,1),data(d,2)]));
            evald=double(eval);
            vars1=zeros(1,dimVar);
            vars1((d-1)*k+j)=1;
            vars1(dataPointNum*k+(mon_num)*k+1+(j-1)*mon_num+i)=1;
            %vars1(dataPointNum*k+(mon_num)*k+1+mon_num+j)=1;
            vars_implicit_and_data_weights=...
                [vars_implicit_and_data_weights;vars1];
            coeffs_monomials=[coeffs_monomials;evald(i)];
            if(beginloop==1)
            vars1=zeros(1,dimVar);
            vars1((d-1)*k+j)=1;
          %  vars1(dataPointNum*k+(mon_num)*k+1+mon_num+j)=1;
            vars_datapoint_and_data_weights=...
                [vars_datapoint_and_data_weights;vars1];
        
            
            end
        
        end 
        
%         ineqPolySys{pcount}.typeCone = -1;
%         ineqPolySys{pcount}.dimVar =dimVar;
%         ineqPolySys{pcount}.degree =2;
%         ineqPolySys{pcount}.noTerms = data_samples*2;
%       
%         ineqPolySys{pcount}.supports =...
%         [vars_implicit_data_and_implicit_weights;vars_datapoint_data_and_implicit_weights];
%         
%         ineqPolySys{pcount}.coef= [ones(data_samples,1);...
%            - coeffs_datapoint_data_and_implicit_weights];
%         pcount=pcount+1;
        
       
        
    end
    beginloop=0;
    coeffs_datapoint_and_data_weights=...
                coeffs_datapoint_and_data_weights+coeffs_monomials;
end


