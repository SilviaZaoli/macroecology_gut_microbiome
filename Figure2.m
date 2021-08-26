%This script produces Fig. 2

load('abd_all_clean.mat') %load empirical abundances
%the variable 'abd' contains one array for each individual, in the
%following order: 10 individuals of BIO-ML, M3, F4, A (pre-travel), B
%(pre-Salmonella), A(post-travel), B(post-Salmonella)
%Array rows correspond to samples. The first column contains the sampling
%day, the following contain the counts for each OTU, the last columnn
%contains the unassigned counts. 

%Estimation of K and sigma for each OTU of each individual
tresh_occup=0;  % no threshold on occupancy
tresh_abd=10^(-4); % threshold on average abundance

for i=1:16
   counts{i,1}=sum(abd{i,1}(:,2:end),2);   %total reads
   relabd{i,1}=abd{i,1}(:,2:end-1)./counts{i,1};   %relative abundance
   meanrelabd{i,1}=mean(relabd{i,1});    %mean relative abundance
   occup{i,1}=sum(relabd{i,1}>0)./size(relabd{i,1},1);   %occupancy
   
   mediasq=nanmean(abd{i,1}(:,2:end-1).*(abd{i,1}(:,2:end-1)-1)./(counts{i,1}.*(counts{i,1}-1)));
   sigma{i,1}=2./(1+meanrelabd{i,1}.^2./(mediasq-meanrelabd{i,1}.^2));
   K{i,1}=2*meanrelabd{i,1}./(2-sigma{i,1});
   id{i,1}=find(meanrelabd{i,1}>tresh_abd & occup{i,1}>tresh_occup);   %OTUs above threshold
end


load('meanPhi_alldataset.mat') %load Phi(T), the intra-individual dissimilarity, for all individulas

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure 2A: Examples of Phi(T) curves
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
x=meanPhi{7,1}(id{7,1},:)';
for i=1:10
plot(smooth(x(:,i),10),'LineWidth',2)
hold on
end
plot(smooth(x(:,21),10),'LineWidth',2)
xlim([0 300])
xlabel('\tau (days)')
ylabel('\Phi(T)')
pbaspect([1.6,1,1])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure 2B: Phi(T) curves converge to theoretically expected value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
i=7
plot(sigma{i,1}(id{i,1})./(4-sigma{i,1}(id{i,1})), nanmean(meanPhi{i,1}(id{i,1},10:end)'),'ko')
hold on
plot(0:1,0:1,'--r')
xlabel('$\frac{\sigma_i}{(4-\sigma_i)}$','Interpreter','latex')
ylabel('\Phi_i^{\infty}')
pbaspect([1.6,1,1])  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure 2C: Normalized Phi(T) curves averaged over `flat' and `increasing'
%OTUs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Normalize Phi by its theoretical asymptotic value
for i=1:16
normPhi{i,1}=meanPhi{i,1}(id{i,1},:)'./(sigma{i,1}(id{i,1})./(4-sigma{i,1}(id{i,1})));
end

%Linear fit to find slopes of Phi(T) for each OTU
slope=cell(16,1);
for i=1:16
   points=normPhi{i,1}(11:end-1,:); %discard initial transient
   for j=1:length(id{i,1})   
   lm=fitlm(1:length(points(:,j)),points(:,j));
    slope{i,1}(j)=lm.Coefficients.Estimate(2);
   end
end

%To define a threshold slope for flat Phi(T), we simulate abundance
%time-series with a SLM model with constant parameters, compute the slope
%of Phi(T) for simulated trajectories and define the threshold as the slope
%such that 95% of the simulatad trajectories have a smaller slope

%***WARNING***: long computation, the result can be loaded un-commenting the following line 
%load('threshold_slope.mat')

x=cell(16,1); %simulated relative abundances
meanPhisim=cell(16,1); %simulated Phi(T)
slopesim=cell(16,1); % slopes of simulated Phi(T)
threshold=zeros(16,1);

for i=1:16
    %parameters
    tau=1; %days
    s=sigma{i,1}(id{i,1});  %sigma parameters for analysed OTU of individual i
    Kk=K{i,1}(id{i,1});     %K parameters for analysed OTU of individual i
    
    %simulate time-series
    dt=0.01; %discretization time step
    T=size(meanPhi{i,1},2)/dt; %Total number of steps to simulate  time-series of the same length of the empirical ones
    N=length(id{i,1}); %#OTUs to simulate
    x{i,1}=zeros(N,T);
    
    dW=sqrt(dt)*randn(N,T); %realization of gaussian white noise
    
    for j=1:N
        R=sqrt(s(j)/tau);
        b=1/tau/Kk(j);
        x{i,1}(j,1)=gamrnd(2./s(j) -1, s(j).*Kk(j)/2,1,1); %initial value from stationary distribution
        for t=1:T
            x{i,1}(j,t+1)=x{i,1}(j,t).*(1+dt/tau+R.*dW(j,t))-b*dt*x{i,1}(j,t).^2;
        end
    end
    
    xdaily=x{i,1}(:,1:1/dt:end); %keep only one value per day
    days=abd{i,1}(:,1)-abd{i,1}(1,1)+1; %days when the individual was sampled
    xdaily(:,setdiff(1:size(xdaily,2),days))=nan; %keep only sampled days
    countssamp=nan(size(xdaily,2),1);
    countssamp(days)=counts{i,1};  %sampling depth (this is for all OUTs, but the analyzed OTU that are simulated represent ~100% of counts)
    xsamp=poissrnd(xdaily.*countssamp');  %poisson sampling according to simulated relative abundance
    
    %compute Phi(T)
    times=1:size(meanPhi{i,1},2);
    for j=1:size(meanPhi{i,1},2)-10
        Phisim=[];
        for t=1:times(end)-times(j)
            t1= t;
            t2= t+times(j);
            n1=xsamp(:,t1)' ;
            n2=xsamp(:,t2)' ;
            try 
                %downsample the sample with larger sampling depth
                n1=equalize(n1,n2); 
                n2=equalize(n2,n1);
                d = n1-n2;
                sm = n1+n2;
                phi=(d.^2-sm)./(sm.^2-sm); %Phi(t, T)
                Phisim=[Phisim  phi'];
            end
        end
        if ~isempty(Phisim)
            meanPhisim{i,1}(:,j)=nanmean(Phisim,2);  %average over t
        else
            meanPhisim{i,1}(:,j)=nan(size(xsamp,1),1);
        end
    end
    normPhisim=meanPhisim{i,1}'./(s./(4-s)); %normalization according to theoretical value
    
    %compute slope of Phi(T)
    points=normPhisim(11:end-1,:); %discard initial transient
    for j=1:length(id{i,1})
        lm=fitlm(1:length(points(:,j)),points(:,j));
        slopesim{i,1}(j)=lm.Coefficients.Estimate(2);
    end
    
    %find threshold slope
    threshold(i)=quantile(slopesim{i,1},0.95);
end
save('threshold_slope.mat')

%***end long computation***




%individuate OTUs with increasing Phi
for i=1:16
   id_slope{i,1}=slope{i,1}>threshold(i);
end

% Figure 2C
i=7; %individual bh
figure %2C
mean1= nanmean(normPhi{i,1}(:,id_slope{i,1}),2); % <Phi(T)> averaged over OTUs with increasing Phi(T)
x1=find(~isnan(mean1(1:end-1)));
if sum(id_slope{i,1})>1 %if there are more than one OTU in this group
    err=nanstd(normPhi{i,1}(:,id_slope{i,1})'); % compute standard deviation 
    x = [x1', fliplr(x1')];
    inBetween = [mean1(x1)'+err(x1), fliplr(mean1(x1)'-err(x1))];
    fill(x, inBetween,[0.5, 0.5,0.5],'EdgeColor','none'); % fill area of std
end
hold on
mean2= nanmean(normPhi{i,1}(:,~id_slope{i,1}),2); % <Phi(T)> averaged over OTUs with flat Phi(T)
x2=find(~isnan(mean2(1:end-1)));
err=nanstd(normPhi{i,1}(:,~id_slope{i,1})'); % compute standard deviation 
x = [x2', fliplr(x2')];
inBetween = [mean2(x2)'+err(x2), fliplr(mean2(x2)'-err(x2))];
fill(x, inBetween,[1 0.4 0.6],'EdgeColor','none'); % fill area of std
plot(x1,mean1(x1),'-k', 'Linewidth',2) %Plot <Phi(T)> increasing
plot(x2, mean2(x2),'-r', 'Linewidth',2) %Plot <Phi(T)> flat
plot(1:size(normPhi{i,1},1), ones(size(normPhi{i,1},1),1), '--k')
xlabel('\tau')
ylabel('Average of normalised \Phi(\tau)')
legend({strcat('Increasing (',num2str(sum(id_slope{i,1})),' OTU)'), strcat('Flat (',num2str(sum(~id_slope{i,1})),' OTU)')},'location','northwest')
alpha(0.5)
xlim([0 max(x1)])
pbaspect([1.6,1,1])

%same figure as for 2C, but for all individuals
figure
for i=1:16
    subplot(4,4,i)
    
    mean1= nanmean(normPhi{i,1}(:,id_slope{i,1}),2);
    x1=find(~isnan(mean1(1:end-1)));
    if sum(id_slope{i,1})>1
        err=nanstd(normPhi{i,1}(:,id_slope{i,1})');
        x = [x1', fliplr(x1')];
        inBetween = [mean1(x1)'+err(x1), fliplr(mean1(x1)'-err(x1))];
        fill(x, inBetween,[0.5, 0.5,0.5],'EdgeColor','none');
       
    end
    hold on
    mean2= nanmean(normPhi{i,1}(:,~id_slope{i,1}),2);
    x2=find(~isnan(mean2(1:end-1)));
    err=nanstd(normPhi{i,1}(:,~id_slope{i,1})');
    x = [x2', fliplr(x2')];
    inBetween = [mean2(x2)'+err(x2), fliplr(mean2(x2)'-err(x2))];
    fill(x, inBetween,[1 0.4 0.6],'EdgeColor','none');
    plot(x1,mean1(x1),'-k', 'Linewidth',2)
    plot(x2, mean2(x2),'-r', 'Linewidth',2)
    plot(1:size(normPhi{i,1},1), ones(size(normPhi{i,1},1),1), '--k')
    xlabel('\tau')
    ylabel('Average of normalised \Phi(\tau)')
    title(individuals{1,i})
    legend({strcat('Increasing (',num2str(sum(id_slope{i,1})),')'), strcat('Flat (',num2str(sum(~id_slope{i,1})),')')},'location','northwest')
    alpha(0.5)  
    xlim([0 max(x1)])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure 2D: percentage of OTUs whose Phi(T) is flat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
perc_stable=zeros(16,1);
for i=1:16
   perc_stable(i)=sum(~id_slope{i,1})/length(id_slope{i,1})*100; 
end
individuals2=individuals;
individuals2{1,13}='A bef';
individuals2{1,14}='A aft';
individuals2{1,15}='B bef';
individuals2{1,16}='B aft';

figure %2D
bar(perc_stable,'Facecolor', [1 0.4 0.6],'EdgeColor','none')
ylabel('Percentage Stable OTU')
xticks(1:16)
xticklabels(individuals2)
pbaspect([1.6,1,1])

