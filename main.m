n=1000;
session=1;
warning('off','MATLAB:ode15s:IntegrationTolNotMet');

matlabpool open;

LB=zeros(1,72);
UB=zeros(1,72);


for i=1:2 
    
    LB(i)=-4; UB(i)=0;

end


for i=3:27 
    
    LB(i)=1; UB(i)=5;

end


for i=28:52
    
    LB(i)=-4; UB(i)=3;

end

for i=53:65
    
    LB(i)=-3; UB(i)=4;

end

for i=66:68
    
    LB(i)=-3; UB(i)=3;
    
end

for i=69:71
    
    LB(i)=-1; UB(i)=4;

end

LB(72)=-3;UB(72)=1;

x=zeros(n,72);
f=zeros(1,n);

spmd
    sc = RandStream('CombRecursive','Seed',session*labindex);
    RandStream.setGlobalStream(sc);
end

parfor i=1:n
  

  disp(i);
  options=gaoptimset('TolFun',1E-6,'Generations',Inf);
  [x(i,:),f(i)] = ga(@myfun,72,[],[],[],[],LB,UB,[],options);

end

csvwrite('ga_x.dat',x);
csvwrite('ga_f.dat',f);

matlabpool close;
exit;

