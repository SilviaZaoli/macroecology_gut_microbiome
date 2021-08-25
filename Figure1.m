load('abd_all_clean.mat')
clear counts 

tresh_occup=0;  
tresh_abd=10^(-4);

for i=1:16
   counts{i,1}=sum(abd{i,1}(:,2:end),2);   %total reads
   relabd{i,1}=abd{i,1}(:,2:end-1)./counts{i,1};   %relative abundance
   meanrelabd{i,1}=mean(relabd{i,1});    %mean relative abundance
   occup{i,1}=sum(relabd{i,1}>0)./size(relabd{i,1},1);   %occupancy
   id{i,1}=find(meanrelabd{i,1}>tresh_abd & occup{i,1}>tresh_occup);   %OTUs above threshold
end

%exaample dynamics
i=11;
figure
y=relabd{i,1}(:,id{i,1});
x=abd{i,1}(:,1);
plot(x,y(:,24),'LineWidth',2)
hold on
plot(x,y(:,1),'LineWidth',2)

plot(x,y(:,11),'LineWidth',2)
plot(x,y(:,25),'LineWidth',2)

xlim([0 400])
set(gca,'yscale','log')
xlabel('time (days)')
ylabel('Relative abundance')
pbaspect([1.6,1,1])

%comparison of two moments
i=11;
figure
%plot(y(100,:),y(300,:),'ko')
plot(relabd{i,1}(100,:),relabd{i,1}(300,:),'ko')
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('Relative abundance at time t_1')
ylabel('Relative abundance at time t_2')
set(gca,'DataAspectRatio', [1 1 1]);

%comparison of different individuals
i1=11;
i2=12;
figure
%plot(y(100,:),y(300,:),'ko')
plot(relabd{i1,1}(100,:),relabd{i2,1}(100,:),'ko')
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('Relative abundance for individual 1')
ylabel('Relative abundance for individual 2')
set(gca,'DataAspectRatio', [1 1 1]);

%possible explanations in terms of trajectories
tau=1;
K=[5*10^(-3),5*10^(-5)];
sigma=[1,1.2];

dt=0.01; %discretization time step
T=400/dt; %Total number of steps to simulate 300days
N=2; %individui da simulare
xLan=zeros(N,T);

dW=sqrt(dt)*randn(N,T);

for j=1:2
    R=sqrt(sigma(j)/tau);
    b=1/tau/K(j);
    
    xLan(j,1)=gamrnd(2./sigma(j) -1, sigma(j).*K(j)/2,1,1);
    for t=1:T
        xLan(j,t+1)=xLan(j,t).*(1+dt/tau+R.*dW(j,t))-b*dt*xLan(j,t).^2;
    end
    
end

x1=xLan(:,1:1/dt:end);

%t=1:40000;
%Kt(1,:)=10^(-4)*exp((t-1)/8000);
%Kt(2,:)=10^(-2)*exp((-t+1)/8000);

Kt(1,:)=[1*10^(-4)*ones(1, 30000) 5*10^(-3)*ones(1, 10000)];
Kt(2,:)=[5*10^(-5)*ones(1, 20000) 5*10^(-7)*ones(1, 20000)];

xLan=zeros(N,T);

for j=1:2
    R=sqrt(sigma(j)/tau);
    
    xLan(j,1)=gamrnd(2./sigma(j) -1, sigma(j).*Kt(j,1)/2,1,1);
    for t=1:T-1
        xLan(j,t+1)=xLan(j,t).*(1+dt/tau+R.*dW(j,t))-(1/tau/Kt(j,t+1))*dt*xLan(j,t).^2;
    end
    
end
x2=xLan(:,1:1/dt:end);

figure
subplot(2,1,1)
plot(x1','LineWidth',2)
xlim([0 400])
set(gca,'yscale','log')
xlabel('time (days)')
ylabel('Relative abundance')
pbaspect([1.6,1,1])

subplot(2,1,2)
plot(x2','LineWidth',2)
xlim([0 400])
set(gca,'yscale','log')
xlabel('time (days)')
ylabel('Relative abundance')
pbaspect([1.6,1,1])

times=1:400;
for j=1:400-100
    Phisim=[];
    for t=1:times(end)-times(j)
        t1= t;
        t2= t+times(j);
       
        Phisim=[Phisim  ((x1(:,t1)-x1(:,t2))./(x1(:,t1)+x1(:,t2))).^2];
        
    end
    meanPhisim{1,1}(:,j)=mean(Phisim,2);
    
end

times=1:400;
for j=1:400-100
    Phisim=[];
    for t=1:times(end)-times(j)
        t1= t;
        t2= t+times(j);
       
        Phisim=[Phisim  ((x2(:,t1)-x2(:,t2))./(x2(:,t1)+x2(:,t2))).^2];
        
    end
    meanPhisim{2,1}(:,j)=mean(Phisim,2);
    
end

figure
subplot(2,1,1)
plot(meanPhisim{1,1}','LineWidth',2)
xlim([0 300])
xlabel('\tau (days)')
ylabel('\Phi(\tau)')
pbaspect([1.6,1,1])

subplot(2,1,2)
plot(meanPhisim{2,1}','LineWidth',2)
xlim([0 300])
xlabel('\tau (days)')
ylabel('\Phi(\tau)')
pbaspect([1.6,1,1])
