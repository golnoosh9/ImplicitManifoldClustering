K2Inds=find(B(:,2)>0.5);
K1Inds=find(B(:,1)>0.5);
figure;
plot(data(K1Inds,1),data(K1Inds,2));
figure;
plot(data(K2Inds,1),data(K2Inds,2));
figure;
plot(data(:,1),data(:,2));
figure;
plot_implicit_2(C(:,1));
figure;
plot_implicit_2(C(:,2));;


