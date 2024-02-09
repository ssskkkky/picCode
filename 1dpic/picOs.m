for icase=1:1:1
dt=0.0001;
stepTotal=62800*2;
xnpoints=64+1;
L=2*pi*10/icase;
kn=1;
parPerCell=80;
k = kn*2*pi/L;
tTotal = stepTotal*dt;
parNum = parPerCell*(xnpoints-1);
T = 1.;
n = 1.;
v0 = 0.;
Fm=@(v) n/sqrt(2*pi*T)*exp(-v^2/(2*T));
Fmp=@(v) n/sqrt(2*pi*T)*(-v/T)*exp(-v^2/(2*T));
pert = 0.5;
parMatrix=zeros(parNum,5);
parMatrix(:,1)=rand(parNum,1)*L+sin(k*parMatrix(:,1))*pert;
parMatrix(:,1)=linspace(0,L-L/parNum,parNum);
parMatrix(:,1)=parMatrix(:,1)+(sin(k*parMatrix(:,1))+cos(k*parMatrix(:,1)))*pert;
parMatrix(:,1)=mod(parMatrix(:,1),L);
parMatrix(:,2)=randn(parNum,1)*sqrt(2*T);
% parMatrix(:,1)=L/2;
% parMatrix(:,2)=1;
parMatrix(:,3)=n*L/parNum;
xList=linspace(0,L,xnpoints);
phi=xList;
dx=L/(xnpoints-1);
icount=0;
Ek=zeros(stepTotal,1);
Ef=zeros(stepTotal,1);
% parAll=zeros(stepTotal,parNum,5);

filter=zeros(1,xnpoints-1);
filter(kn+1)=1;
filter(end+1-kn)=1;

for it=0:dt:tTotal-dt
    icount=icount+1;
    parMatrix(:,4)=floor(parMatrix(:,1)/dx)+1;
    parMatrix(:,5)=parMatrix(:,1)/dx-floor(parMatrix(:,1)/dx);
    phi(:)=0;
    for ii=1:1:parNum
        % phi(parMatrix(ii,4))=(1-parMatrix(ii,5))/(parNum*dx)*(n*L);
        % phi(parMatrix(ii,4)+1)=parMatrix(ii,5)/(parNum*dx)*(n*L);
        phi(parMatrix(ii,4))=phi(parMatrix(ii,4))+(1-parMatrix(ii,5))/(dx)*parMatrix(ii,3);
        phi(parMatrix(ii,4)+1)=phi(parMatrix(ii,4)+1)+parMatrix(ii,5)/(dx)*parMatrix(ii,3);
    end
    phi(1)=phi(1)+phi(end);
    rho=phi;
    phi(1:end-1)=poissonSolver(-phi(1:end-1),filter)/(2*pi/L)^2;
    phi(end)=phi(1);
    % phi=real(phi);
    Ex=gradient(real(phi))/(dx);
    Ex=-[real(phi(2)-phi(end-1))/(2*dx),Ex(2:end-1),real(phi(2)-phi(end-1))/(2*dx)];

    % Ex=sin(xList/L*2*pi);
    % Exn=poissonSolver3(rho)*dx;


    parMatrix(:,2)=parMatrix(:,2)+interp1(xList,Ex,parMatrix(:,1))*dt;
    parMatrix(:,1)=mod(parMatrix(:,1)+parMatrix(:,2)*dt,L);
    Ek(icount)=sum(parMatrix(:,2).^2/2.*parMatrix(:,3));
    Ef(icount)=sum(Ex(1:end-1).^2/2)*dx;
    % parAll(icount,:,:)=parMatrix;
    if (mod(icount,100)==0)
    % figure(1)    
    plot(Ex)
    ylim([-1,1])
    % figure(2)
    % plot(xList,real(phi),xList,imag(phi))
    % ylim([-3,3])
    drawnow
    end
end
figure()
tList=(1:stepTotal);
plot(tList,Ek,tList,Ef,tList,Ek+Ef)
figure()
logEf=log(Ef);
[lMaxI,P]=islocalmax(logEf,'MinSeparation',10000,'MinProminence',0.2,'ProminenceWindow',[10000,10000]);
plot(tList,logEf,tList(lMaxI),logEf(lMaxI),'r*')
maxtList=tList(lMaxI);
maxDiffList=logEf(lMaxI);
T=(maxtList(2)-maxtList(1))*dt*2;
omega=2*pi/T;
gamma=(maxDiffList(2)-maxDiffList(1))/(T/2);                                                                                                                                                                                                                                   
omega+gamma*1j
drawnow
omeAll(icase)=omega+gamma*1j;
end
omeAll
% figure(2)
% plot((1:stepTotal+1),Ef)
% figure()
% plot(Ex)
% figure()
% plot((1:xnpoints),rho-parPerCell/(dx),'b*',(1:xnpoints),gradient(Ex)/dx,'c*')
% figure()
% plot(phi)
% for i=1:1:stepTotal
% plot(parAll(i,:,1),parAll(i,:,2),'b.')
% drawnow
% end
