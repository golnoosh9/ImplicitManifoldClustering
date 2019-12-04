%%%%var order: b(1,1),b(1,2),...,b(1,k),b(2,1)...
%%%%(dataPointNum*k)||C1,...,Ck(mon_num*k)||eps||X_raise_1,...,X_raise_k (mon_num*k)
%%%% ||Lk (k)||vecEps (mon_num)
%%%%const: -eps<Ck.X_raise_k<eps
for i=1:k
    
    vars_implicit_curvei=[];
    g_implicit_curvei=[];
    for j=1:mon_num
        gi=zeros(1,varNum+1);
        gi(dataPointNum*k+mon_num*k+1+(mon_num)*(i-1)+j)=1;%implicit_monomial
        gi(dataPointNum*k+(mon_num)*(i-1)+j)=1;%%i-th curve monomial
        gi(varNum+1)=-100%%for now fixed instead of epsilon
       g_implicit_curvei=[g_implicit_curvei;gi];
    end
    g0=zeros(1,varNum);
    g0(varNum+1)=1;
    pop.G{ineq_count}=[g_implicit_curvei;g0];
    ineq_count=ineq_count+1;
    
    
end


