function Eg=poissonSolver3(rho)
Eg=zeros(length(rho),1);
Eg(1)=(rho(1)-rho(end-1))/2;
Eg(end)=Eg(1);
for i=2:length(Eg)
   Eg(i)=Eg(i-1)+(rho(i)-rho(i-1))/2;
end
Eg=Eg-mean(Eg);