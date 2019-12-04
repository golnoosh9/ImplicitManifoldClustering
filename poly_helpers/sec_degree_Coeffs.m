varst7=[];

for xx=dataPointNum*k+(i-1)*mon_num+1:dataPointNum*k+(i-1)*mon_num+mon_num
    varst2=zeros(1,dimVar);
    varst2(xx)=1;
    %varst2(mon_num+indj)=1; for weight
    varst7=[varst7;varst2];    
    
end
