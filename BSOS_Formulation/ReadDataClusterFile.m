fid = fopen('spiral.txt','r');

tline = fgetl(fid);
ss=size(tline);
dCount=1;
while (  tline~=-1)
X=str2num(tline); 
data(dCount,1)=(X(1));
data(dCount,2)=(X(2));
tline = fgetl(fid);
ss=size(tline);
dCount=dCount+1;
end

