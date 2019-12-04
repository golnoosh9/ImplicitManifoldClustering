%%%%var order: b(1,1),b(1,2),...,b(1,k),b(2,1)...
%%%%(dataPointNum*k)||C1,...,Ck(mon_num*k)||eps||X_raise_1,...,X_raise_k (mon_num*k)
%%%% ||Lk (k)||vecEps (mon_num)
%%%%const: -eps<Ck.X_raise_k<eps
for i=1:k
    
    vars_implicit_curvei=[];
    for j=1:mon_num
        vars1=zeros(1,dimVar);
        vars1(dataPointNum*k+mon_num*k+1+(mon_num)*(i-1)+j)=1;%implicit_monomial
        vars1(dataPointNum*k+(mon_num)*(i-1)+j)=1;%%i-th curve monomial
        vars_implicit_curvei=[vars_implicit_curvei;vars1];
    end
    ineqPolySys{pcount}.typeCone = 1;
    ineqPolySys{pcount}.dimVar =dimVar;
    ineqPolySys{pcount}.degree =2;
    ineqPolySys{pcount}.noTerms = mon_num+1;
    vars_eps=zeros(1,dimVar);
    vars_eps(dataPointNum*k+mon_num*k+1)=1;
    ineqPolySys{pcount}.supports = [vars_implicit_curvei;vars_eps];
    
    ineqPolySys{pcount}.coef     = [-ones(mon_num,1);1];
    pcount=pcount+1;
    
    ineqPolySys{pcount}.typeCone = 1;
    ineqPolySys{pcount}.dimVar =dimVar;
    ineqPolySys{pcount}.degree =2;
    ineqPolySys{pcount}.noTerms = mon_num+1;
    vars_eps=zeros(1,dimVar);
    vars_eps(dataPointNum*k+mon_num*k+1)=1;
    ineqPolySys{pcount}.supports = [vars_implicit_curvei;vars_eps];
    
    ineqPolySys{pcount}.coef     = [ones(mon_num,1);1];
    pcount=pcount+1;
    
    
end


