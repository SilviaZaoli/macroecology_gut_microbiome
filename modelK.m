%This script simulates abundance time-series according to an SLM with carrying capacity K jumping, and performs anlayses on the simulated data

load('abd_all_clean.mat') %load empirical abundances
%the variable 'abd' contains one array for each individual, in the
%following order: 10 individuals of BIO-ML, M3, F4, A (pre-travel), B
%(pre-Salmonella), A(post-travel), B(post-Salmonella)
%Array rows correspond to samples. The first column contains the sampling
%day, the following contain the counts for each OTU, the last columnn
%contains the unassigned counts. 

%Estimation of K and sigma 
for i=1:16
   counts{i,1}=sum(abd{i,1}(:,2:end),2);   %total reads
   relabd{i,1}=abd{i,1}(:,2:end-1)./counts{i,1};   %relative abundance
   meanrelabd{i,1}=mean(relabd{i,1});    %mean relative abundance
   meansq=nanmean(abd{i,1}(:,2:end-1).*(abd{i,1}(:,2:end-1)-1)./(counts{i,1}.*(counts{i,1}-1)));
   sigma{i,1}=2./(1+meanrelabd{i,1}.^2./(meansq-meanrelabd{i,1}.^2));
   K{i,1}=2*meanrelabd{i,1}./(2-sigma{i,1});
   K{i,1}(isnan(sigma{i,1}) & meanrelabd{i,1}==0)=0;
   occup{i,1}=sum(relabd{i,1}>0)./size(relabd{i,1},1);   %occupancy
   id{i,1}=find(sigma{i,1}>0 & sigma{i,1}<Inf & occup{i,1}>0.2);   %OTUs to analyse
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Analysis for Fig. S?: var(xi) does not depend on K
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find species common to all BIOML individuals
id_comm=intersect(id{1,1},id{2,1});
for i=3:10
   id_comm=intersect(id_comm, id{i,1});
end

%Compute \bar{K} for each OTU and the value of xi for that OTU in each
%individual
Kbar=zeros(length(id_comm),1);
xi=zeros(length(id_comm),10);
for i=1:length(id_comm)
    ks=[]; 
   for j=1:10
       ks=[ks, K{j,1}(id_comm(i))];
   end
    Kbar(i)=mean(ks); %mean K for each OTU across all individuals
    xi(i,:)=ks/Kbar(i);
end

%Plot of values of xi for each OTU against \bar{K}
figure
for i=1:length(id_comm)
plot(Kbar(i), xi(i,:),'ko')
hold on 
end
set(gca, 'xscale','log')
ylabel('\xi_i')
xlabel('$\bar{K}$','interpreter','latex')

%Fig. S?A: Plot of var(xi) for each OTU against \bar{K} (with its average in logarithmic bins)
varxi=var(xi');
binedges=10.^[-5,-4,-3,-2,-1,0];
for i=1:5
   binmean(i)=mean(varxi(Kbar>binedges(i)& Kbar<binedges(i+1)));
end

figure
plot(Kbar, varxi,'ok')
hold on
for i=1:5
   plot(10^(log10(binedges(i))+(log10(binedges(i+1))-log10(binedges(i)))/2),binmean(i),'r.','Markersize',20) 
end
set(gca, 'xscale','log')
ylabel('var(\xi_i)')
xlabel('$\bar{K}_i$','interpreter','latex')
pbaspect([1.6 1 1])


%To show that the concave pattern is due to observation bias, we generate
%the same plot for OTU that have all the same var(chi) and that are
%observed in all 10 individuals
N=10; %individuals
M=10000; %OTUs (not all will be observed in all 10 individuals)
barks=lognrnd(-12,3,M,1); %\bar{K} for each OTU

m = 1; % mean of xi
v = 5; % variance of xi
mu = log((m^2)/sqrt(v+m^2));  %mean of underlying gaussian
s = sqrt(log(v/(m^2)+1));     %variance of underlying gaussinan
ks=barks.*(lognrnd(mu,s,M,10)); %values of K for each OTU in each individual

id_obs=[]; %IDs of OTUs that are observed in all 10 individuals (they have K>10^(-5) in all individuals)
for i=1:M
   if sum(ks(i,:)>10^(-5))==N & sum(ks(i,:)<1)==N %we also exclude cases in which an extracted K is larger than 1
       id_obs=[id_obs;i];
   end
end

%for the selected OTUs, we perform the same analysis performed on the data
%(compute K, xi, var(xi)
Kbar=zeros(length(id_obs),1); 
xi=zeros(length(id_obs),10);
for i=1:length(id_obs)
    Kbar(i)=mean(ks(id_obs(i),:));
    xi(i,:)=ks(id_obs(i),:)/Kbar(i);
end

%Plot of values of xi for each OTU against \bar{K}
figure
for i=1:length(id_comm)
plot(Kbar(i), xi(i,:),'ko')
hold on 
end
set(gca, 'xscale','log')
ylabel('\xi_i')
xlabel('$\bar{K}$','interpreter','latex')

varxi=var(xi');
binedges=10.^[-6,-5,-4,-3,-2,-1,0];
for i=1:6
   binmean(i)=mean(varxi(Kbar>binedges(i)&Kbar<binedges(i+1)));
end

%Fig. S?B: Plot of var(xi) for each OTU against \bar{K} (with its average in logarithmic bins)

figure
plot(Kbar, varxi,'ok')
hold on
for i=1:6
   plot(10^(log10(binedges(i))+(log10(binedges(i+1))-log10(binedges(i)))/2),binmean(i),'r.','Markersize',20) 
end
set(gca, 'xscale','log')
ylabel('var(\xi_i)')
xlabel('$\bar{K}_i$','interpreter','latex')
pbaspect([1.6 1 1])
axis([10^-5 1 0 10])

clear barK

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fig. 5A: Repeat the dissimilarity analysis of Figs 2C-D on the simulated data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Parameters for model of carrying capacity:

m = 1; % mean of xi
v = 2; % variance of xi
mu = log((m^2)/sqrt(v+m^2)); %mean of underlying gaussian 
s = sqrt(log(v/(m^2)+1));  %variance of underlying gaussian

%For each OTU \bar{K} is obtained averaging the estimated K among individuals of the
%same dataset
%1) BIO-ML
Ks=[];
for i=1:10
   Ks= [Ks; K{i,1}];
end
for i=1:10
   barK{i,1}= nanmean(Ks);
end
%2) MP
barK{11,1}=nanmean([K{11,1};K{12,1}]);
barK{12,1}=barK{11,1};
%3) DA (B post-Salmonella treated separately because it has different composition from pre-Salmonella)
barK{13,1}=nanmean([K{13,1};K{14,1};K{15,1}]);
barK{14,1}=barK{13,1};
barK{15,1}=barK{13,1};
barK{16,1}=K{16,1};

%fraction of OTU for which I observe a jump in the observation time-window (from analysis in Fig. 3C)
fracjump=[0.1, 0.23, 0.1, 0.15, 0.16,0.18, 0.3, 0.12, 0.25,0.22, 0.68, 0.5, 0.25, 0.52, 0.41,0.67]; 



rng(1) %seed, for reproducibility

% We simulate time-series
for i=1:16 
    Ndays=abd{i,1}(end,1)-abd{i,1}(1,1); %time-window to simulate
    N=length(id{i,1}); %number of OTUs to simulate 
    ss=sigma{i,1}(id{i,1}); % the sigma of each OTU is taken as the empirical one
    %kk=K{i,1}(id{i,1}); % barK for each OTU is taken as the empirical one
    kk=barK{i,1}(id{i,1}); % barK for each OTU 
    
    
    parexp=-Ndays/log(1-fracjump(i)); %parameter of exponential distribution of jumping times such that the correct fraction of jumping OTUs is observed
    jumptimes1=exprnd(parexp, N,1); %time of first jump
    jumptimes2=jumptimes1+exprnd(parexp, N,1); %time of second jump (used only if it falls in the observation window, rarely)
    jumpingOTUs=find(jumptimes1<Ndays); %number of OTU with at least one jump
    
    Ks1=kk.*(lognrnd(mu,s, 1,N)); %values of K for each OTU before first jump
    Ks2=kk.*(lognrnd(mu,s, 1,N)); %values of K for each OTU between first and second jump
    Ks3=kk.*(lognrnd(mu,s, 1,N)); %values of K for each OTU between second and third jump
    
    tau=1; %parameter of SLM
    dt=0.01; %discretization time step
    T=Ndays/dt; %simulation steps
    
    dW=sqrt(dt)*randn(N,T); %realizations of the Gaussian white noise
    
    x=zeros(N,T); %relative abundances to simulate according to SLM
    for j=1:N
        R=sqrt(ss(j)/tau);
        x(j,1)=gamrnd(2./ss(j) -1, ss(j).*Ks1(j)/2,1,1);
        for t=1:min(round(jumptimes1(j))/dt, T)
            b=1/tau/Ks1(j);
            x(j,t+1)=x(j,t).*(1+dt/tau+R.*dW(j,t))-b*dt*x(j,t).^2;
        end
        for t=t+1:min(round(jumptimes2(j))/dt, T)
            b=1/tau/Ks2(j);
            x(j,t+1)=x(j,t).*(1+dt/tau+R.*dW(j,t))-b*dt*x(j,t).^2;
        end
        for t=t+1:T
            b=1/tau/Ks3(j);
            x(j,t+1)=x(j,t).*(1+dt/tau+R.*dW(j,t))-b*dt*x(j,t).^2;
        end
    end
    
    x=x(:,1:1/dt:end); %one value per day
    Nsamp=5.5*10^4; %sampling depth
    xsamp{i,1}=poissrnd(Nsamp*x); %sampled counts
    xsamp{i,1}(:,setdiff(1:Ndays,abd{i,1}(:,1)-abd{i,1}(1,1)+1))=nan; %exclude values of day that were not sampled in the data
end

%We compute normalized Phi(T)
for i=1:16
    meanrelabd=nanmean(xsamp{i,1}/Nsamp,2);
    idXphi=find(meanrelabd>10^(-4)); %we analyze only OTUs with average relative abundance above threshold
    xsamp2{i,1}=xsamp{i,1}(idXphi,:);
    Ndays=abd{i,1}(end,1)-abd{i,1}(1,1); %time-window
    times=1:Ndays;
    for j=1:Ndays-10
        Phisim=[];
        for t=1:times(end)-times(j)
            t1= t;
            t2= t+times(j);
            n1=xsamp2{i,1}(:,t1)' ;
            n2=xsamp2{i,1}(:,t2)' ;
            try
                d = n1-n2;
                sm = n1+n2;
                phi=(d.^2-sm)./(sm.^2-sm);
                Phisim=[Phisim  phi'];
            end
        end
        if ~isempty(Phisim)
            meanPhisim(:,j)=nanmean(Phisim,2);
        else
            meanPhisim(:,j)=nan(size(xsamp2{i,1},1),1);
        end
        
    end
    
    ss=sigma{i,1}(id{i,1}); % sigma of OTUs
    normPhisim{i,1}=meanPhisim'./(ss(idXphi)./(4-ss(idXphi))); %normalized dissimilarity according to expected value
    clear meanPhisim
end

%Find OTUs with increasing dissimilarity
load('threshold_slope.mat')
for i=1:16
    points=normPhisim{i,1}(11:end-1,:);
    
    for j=1:size(points,2)
        lm=fitlm(1:length(points(:,j)),points(:,j));
        slopesim(j)=lm.Coefficients.Estimate(2);
    end
    
    id_slope{i,1}=slopesim>threshold(i); %OTUs with dissimilarity slope above threshold
    clear slopesim
end

%Fig 5A: Phi(T) for individual 'bh', averaged over stable OTUs and over
%non-stable OTUs
i=7;

figure
mean1= nanmean(normPhisim{i,1}(:,id_slope{i,1}),2);
x1=find(~isnan(mean1(1:end-1)));
if sum(id_slope{i,1})>1
    err=nanstd(normPhisim{i,1}(:,id_slope{i,1})');
    x = [x1', fliplr(x1')];
    inBetween = [mean1(x1)'+err(x1), fliplr(mean1(x1)'-err(x1))];
    fill(x, inBetween,[0.5, 0.5,0.5],'EdgeColor','none');  
end
hold on
mean2= nanmean(normPhisim{i,1}(:,~id_slope{i,1}),2);
x2=find(~isnan(mean2(1:end-1)));
err=nanstd(normPhisim{i,1}(:,~id_slope{i,1})');
x = [x2', fliplr(x2')];
inBetween = [mean2(x2)'+err(x2), fliplr(mean2(x2)'-err(x2))];
fill(x, inBetween,[1 0.4 0.6],'EdgeColor','none');
plot(x1,mean1(x1),'-k', 'Linewidth',2)
plot(x2, mean2(x2),'-r', 'Linewidth',2)
plot(1:size(normPhisim{i,1},1), ones(size(normPhisim{i,1},1),1), '--k')
xlabel('\tau')
ylabel('Average of normalised \Phi(\tau)')
title(individuals{1,i})
legend({strcat('Increasing (',num2str(sum(id_slope{i,1})),')'), strcat('Flat (',num2str(sum(~id_slope{i,1})),')')},'location','northwest')
alpha(0.5)
xlim([0 max(x1)])
pbaspect([1.6 1 1])

%Fig. S? 
perc_stable=zeros(16,1);
for i=1:16
   perc_stable(i)=sum(~id_slope{i,1})/length(id_slope{i,1})*100; 
end
individuals2=individuals;
individuals2{1,13}='A bef';
individuals2{1,14}='A aft';
individuals2{1,15}='B bef';
individuals2{1,16}='B aft';

figure
bar(perc_stable,'Facecolor', [1 0.4 0.6],'EdgeColor','none')
ylabel('Percentage Stable OTU')
xticks(1:16)
xticklabels(individuals2)
pbaspect([1.6,1,1])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Analysis for Fig. 5B: reproduction of emipirical correlations in time and
%across individuals with simulated time-series
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Correlations in time
for i=1:16 %We divide time series in two and estimate K and sigma for each half
    abd1=xsamp{i,1}(:,1:round(size(xsamp{i,1},2)/2));
    abd2=xsamp{i,1}(:,round(size(xsamp{i,1},2)/2)+1:end);
    counts1=nansum(abd1);
    counts2=nansum(abd2);
    relabd1=abd1./counts1;
    relabd2=abd2./counts2;
    meanrelabd1=nanmean(relabd1,2);    %mean relative abundance
    meanrelabd2=nanmean(relabd2,2);
    mediasq1=nanmean(abd1.*(abd1-1)./(counts1.*(counts1-1)),2);
    mediasq2=nanmean(abd2.*(abd2-1)./(counts2.*(counts2-1)),2);
    sigma1{i,1}=2./(1+meanrelabd1.^2./(mediasq1-meanrelabd1.^2));
    sigma2{i,1}=2./(1+meanrelabd2.^2./(mediasq2-meanrelabd2.^2));
    K1{i,1}=2*meanrelabd1./(2-sigma1{i,1});
    K2{i,1}=2*meanrelabd2./(2-sigma2{i,1});
end

corKt=zeros(16,1);
corst=zeros(16,1);
for i=1:16
    x=K1{i,1}';
    y=K2{i,1}';
    a=sigma1{i,1}';
    b=sigma2{i,1}';
    ids=a<Inf & b<Inf & a>0 & b>0; %keep only OTUs for which sigma could be estimated
    corst(i)=corr(a(ids)', b(ids)');
    corKt(i)=corr(log(x(ids)'), log(y(ids)'));
end

%correlazione K tra individui
for i=1:15 %estimate K and sigma for each OTU in each individual
   counts=nansum(xsamp{i,1});   %total reads
   relabd=xsamp{i,1}./counts;   %relative abundance
   meanrelabd=nanmean(relabd,2);    %mean relative abundance
  
   meansq=nanmean(xsamp{i,1}.*(xsamp{i,1}-1)./(counts.*(counts-1)),2);
   sigmaest{i,1}=2./(1+meanrelabd.^2./(meansq-meanrelabd.^2));
   Kest{i,1}=2*meanrelabd./(2-sigmaest{i,1});
end


corK=[];
cors=[];
for i=1:10 %between all pairs of BIO-ML individuals
    for j=i+1:10
        ids=intersect(id{i,1},id{j,1}); %OTU common to the two individuals
        [~,idsi]=ismember(ids,id{i,1}); %their position in parameters of individual i
        [~,idsj]=ismember(ids,id{j,1}); %their position in parameters of individual j
        x=Kest{i,1}(idsi)';
        y=Kest{j,1}(idsj)';
        a=sigmaest{i,1}(idsi)';
        b=sigmaest{j,1}(idsj)';
        ids=a<Inf & b<Inf & a>0 & b>0; %keep only OTUs for which sigma could be estimated
        corK=[corK; corr(log(x(ids))', log(y(ids))')];
        cors=[cors; corr(a(ids)', b(ids)')];
    end
end
%between M3 and F4
i=11;
j=12;
ids=intersect(id{i,1},id{j,1});
[~,idsi]=ismember(ids,id{i,1});
[~,idsj]=ismember(ids,id{j,1});
x=Kest{i,1}(idsi)';
y=Kest{j,1}(idsj)';
a=sigmaest{i,1}(idsi)';
b=sigmaest{j,1}(idsj)';
ids=a<Inf & b<Inf & a>0 & b>0;
corK=[corK; corr(log(x(ids))', log(y(ids))')];
cors=[cors; corr(a(ids)', b(ids)')];
%between A pre-travel and B pre-Salmonella
i=13;
j=14;
ids=intersect(id{i,1},id{j,1});
[~,idsi]=ismember(ids,id{i,1});
[~,idsj]=ismember(ids,id{j,1});
x=Kest{i,1}(idsi)';
y=Kest{j,1}(idsj)';
a=sigmaest{i,1}(idsi)';
b=sigmaest{j,1}(idsj)';
ids=a<Inf & b<Inf & a>0 & b>0;
corK=[corK; corr(log(x(ids))', log(y(ids))')];
cors=[cors; corr(a(ids)', b(ids)')];
%Between B pre-Salmonella and A post-travel   
i=14;
j=15;
ids=intersect(id{i,1},id{j,1});
[~,idsi]=ismember(ids,id{i,1});
[~,idsj]=ismember(ids,id{j,1});
x=Kest{i,1}(idsi)';
y=Kest{j,1}(idsj)';
a=sigmaest{i,1}(idsi)';
b=sigmaest{j,1}(idsj)';
ids=a<Inf & b<Inf & a>0 & b>0;
corK=[corK; corr(log(x(ids))', log(y(ids))')];
cors=[cors; corr(a(ids)', b(ids)')];
    

%figure 5B
J(:,:,2)=[ [corKt; nan(length(corK)-length(corKt),1)] corK [corst; nan(length(cors)-length(corst),1)] cors];
load('corrData.mat')
J(:,:,1)=[ [corKt; nan(length(corK)-length(corKt),1)] corK [corst; nan(length(cors)-length(corst),1)] cors];
H=iosr.statistics.boxPlot(J)
ylabel('Correlation')
H.theme='colorall'
H.boxwidth=0.2;
pbaspect([1.6 1 1])
ax = gca;
ax.TickLabelInterpreter='tex';
xticklabels({'K (time)', 'K (cross)','\sigma (time)','\sigma (cross)'})
ylim([0 1])

%Figura S?

figure
i=7;
a=sigma1{i,1};
b=sigma2{i,1};
pos=a>0 & b>0;
x=K1{i,1};
y=K2{i,1};
subplot(2,2,1)
hold on
s=scatter(x(pos),y( pos),'ko')
s.MarkerEdgeAlpha = .5;
plot(10^(-6):10^(-5):1,10^(-6):10^(-5):1,'r--' )
xlabel('K first half')
ylabel('K second half')
set(gca, 'xscale','log','yscale','log')
pbaspect([1,1,1])
axis([1*10^(-6) 1 1*10^(-6) 1])

subplot(2,2,2)
hold on
s=scatter(a,b,'ko')
s.MarkerEdgeAlpha = .5;
plot(0:0.1:2,0:0.1:2,'r--' )
xlabel('\sigma first half')
ylabel('\sigma second half')
pbaspect([1,1,1])
axis([0 2 0 2])

i=2;
j=7;
ids=intersect(id{i,1},id{j,1});
[~,idsi]=ismember(ids,id{i,1});
[~,idsj]=ismember(ids,id{j,1});
x=Kest{i,1}(idsi)';
y=Kest{j,1}(idsj)';
a=sigmaest{i,1}(idsi);
b=sigmaest{j,1}(idsj);
pos=a>0 & b>0;


subplot(2,2,3)
hold on
s=scatter(x(pos),y( pos),'ko')
s.MarkerEdgeAlpha = .5;
plot(10^(-6):10^(-5):1,10^(-6):10^(-5):1,'r--' )
xlabel('K "am" ')
ylabel('K "bh"')
set(gca, 'xscale','log','yscale','log')
pbaspect([1,1,1])
axis([1*10^(-6) 1 1*10^(-6) 1])

subplot(2,2,4)
hold on
s=scatter(a,b,'ko')
s.MarkerEdgeAlpha = .5;
plot(0:0.1:2,0:0.1:2,'r--' )
xlabel('\sigma "am" ')
ylabel('\sigma "bh"')
pbaspect([1,1,1])
axis([0 2 0 2])










