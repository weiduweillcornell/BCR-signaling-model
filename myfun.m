function F=myfun(x)

data.file=importdata('time_course_data.csv');
data=data.file.data;

T=30;
time_point=[0 0.5 1 3 5 10 20 30];

parameter=10.^x;

options=odeset('events',@events_time);
tic;
[t,y] = ode15s(@BCR_model,[0 T],[0 0 0 0 0 0 parameter(72) 0 0 0 0 0 0 0 0],options,parameter(1:71));
   
   T_pkc=2.5; k_dapk=parameter(60);
   DAG=y(:,9); RAS=y(:,12);pERK=y(:,15);pAKT=y(:,10);
   PKC=(1+k_dapk*(DAG+T_pkc)-sqrt((1+k_dapk*(DAG+T_pkc)).^2-4*k_dapk^2*DAG*T_pkc))./(2*k_dapk); 
   T_raf=1.1; k_raf=parameter(65);K_erra=parameter(69); K_akra=parameter(70);
   RAF=(1+k_raf*(RAS+T_raf)-sqrt((1+k_raf*(RAS+T_raf)).^2-4*k_raf^2*RAS*T_raf))./(2*k_raf)./(1+K_erra*pERK+K_akra*pAKT);
   y_new=[y PKC RAF];
   y_norm=bsxfun(@rdivide, y_new, max(y_new));
   y_norm(isnan(y_norm))=0;
   inter_norm=interp1(t,y_norm,time_point,'spline'); 
   index=[1 6 16 10 17 13 15];
 
F=sum(sum((inter_norm(:,index)-data(:,2+(1:7))).^2));


