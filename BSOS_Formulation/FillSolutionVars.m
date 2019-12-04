%%%%var order: b(1,1),b(1,2),...,b(1,k),b(2,1)...
%%%%(dataPointNum*k)||C1,...,Ck(mon_num*k)||eps||X_raise (mon_num)||b_imp_k
%%%%(k) ||Lk (k)||vecEps (mon_num)
sols=psol.YY{1};
sols(1)=[];
B=zeros(dataPointNum,k);
currentInd=1;
C=zeros(mon_num,k);

for i=1:dataPointNum
   B(i,:)= sols(currentInd:currentInd+k-1);
   currentInd=currentInd+k;
end

for i=1:k
    C(:,i)=sols(currentInd:currentInd+mon_num-1);
    currentInd=currentInd+mon_num;
    
end
%%33
eps=sols(currentInd);
currentInd=currentInd+1;
X_implicit=sols(currentInd:currentInd+mon_num*k-1);
currentInd=currentInd+mon_num*k;

lowerInserts=sols(currentInd:currentInd+k-1);
currentInd=currentInd+k;

vectorEpsilons=sols(currentInd:currentInd+mon_num-1);
currentInd=currentInd+mon_num;


