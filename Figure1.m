%This script produces Figure 1

load('abd_all_clean.mat') %load empirical abundances
%the variable 'abd' contains one array for each individual, in the
%following order: 10 individuals of BIO-ML, M3, F4, A (pre-travel), B
%(pre-Salmonella), A(post-travel), B(post-Salmonella)
%Array rows correspond to samples. The first column contains the sampling
%day, the following contain the counts for each OTU, the last columnn
%contains the unassigned counts. 


% Panel A: Examples of trajectories of relative abundances from individual
% M3

i=11; %individual M3
tresh_abd=10^(-4); % we select abundant OTUs

counts=sum(abd{i,1}(:,2:end),2);   %total reads
relabd=abd{i,1}(:,2:end-1)./counts;   %relative abundance
meanrelabd=mean(relabd);    %mean relative abundance
id=find(meanrelabd>tresh_abd);   %OTUs above threshold

figure
x=abd{i,1}(:,1); % days
y=relabd(:,id);
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


%Panels B and C: Examples of trajectories simulated according to the SLM
%with either constant parameters or a carrying capacity that changes in
%time

tau=1; %time-scale of SLM
K=[5*10^(-3),5*10^(-5)]; %Carrying capacity of the two individuals
sigma=[1,1.2]; % Noise amplitudes of the two individuals

dt=0.01; %discretization time step
T=400/dt; %Total number of steps to simulate 400days
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

x1=xLan(:,1:1/dt:end); %keep only one abundance value per day

% Carrying capacities changing in time
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

%Paneld D and E: Dissimilarity curves for the trajectories above
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
