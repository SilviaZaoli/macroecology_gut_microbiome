%This scripts produced Figure 3 and related supplementary figures

load('abd_all_clean.mat') %load empirical abundances
%the variable 'abd' contains one array for each individual, in the
%following order: 10 individuals of BIO-ML, M3, F4, A (pre-travel), B
%(pre-Salmonella), A(post-travel), B(post-Salmonella)
%Array rows correspond to samples. The first column contains the sampling
%day, the following contain the counts for each OTU, the last columnn
%contains the unassigned counts. 

thresh_occup=0.2; %threshold on occupancy
for i=1:16
   counts{i,1}=sum(abd{i,1}(:,2:end),2);   %total reads
   relabd{i,1}=abd{i,1}(:,2:end-1)./counts{i,1};   %relative abundance
   meanrelabd{i,1}=mean(relabd{i,1});    %mean relative abundance
   occup{i,1}=sum(relabd{i,1}>0)./size(relabd{i,1},1);   %occupancy
   mediasq=nanmean(abd{i,1}(:,2:end-1).*(abd{i,1}(:,2:end-1)-1)./(counts{i,1}.*(counts{i,1}-1)));
   sigma{i,1}=2./(1+meanrelabd{i,1}.^2./(mediasq-meanrelabd{i,1}.^2));
   K{i,1}=2*meanrelabd{i,1}./(2-sigma{i,1});
   id{i,1}=find(sigma{i,1}>0 & sigma{i,1}<Inf & occup{i,1}>thresh_occup);   %OTUs above threshold
end

%Load Phi_i^{a,b}(T): for each pair of individuals of the same dataset, it
%contains the values of Phi_i{a,b}(T) for each value of T
load('meanPhicross_alldataset2.mat')
%I average it over T, as it does not depend on T, to obtain one value per OTU
meanPhicross2=cell(47,1);
for i=1:47
   meanPhicross2{i,1}=nanmean(meanPhicross{i,1},2)';
end

% Computation of the theoretical expected value of Phi_i{a,b} based on the
% estimated K and sigma parameters in the two individuals
%***WARNING***: long computation, to load the result uncomment the
%following line
%load('tval_all_threshoccup.mat')
tvalcross=cell(47,1);
tvalcrossK=cell(47,1);
tvalcrossS=cell(47,1);
tvalcrossKS=cell(47,1);
for j=1:47
    for i=1:length(ids{j,1})
        s1=sigma{inds(j,1),1}(ids{j,1}(i));
        s2=sigma{inds(j,2),1}(ids{j,1}(i));
        K1=K{inds(j,1),1}(ids{j,1}(i));
        K2=K{inds(j,2),1}(ids{j,1}(i));
        tvalcross{j,1}(i)=theoreticalvaluecross(s1, s2, K1, K2);
        tvalcrossS{j,1}(i)=theoreticalvaluecross((s1+s2)/2, (s1+s2)/2, K1, K2 );
        tvalcrossK{j,1}(i)=theoreticalvaluecross(s1, s2, (K1+K2)/2, (K1+K2)/2);
        tvalcrossKS{j,1}(i)= (s1+s2)/2 /(4-(s1+s2)/2); 
    end
    
end

%For each pair, compute ids of OTUs common to the two individuals
ids=cell(47,1);
ids{1,1}=intersect(id{11,1},id{12,1});
ids{2,1}=intersect(id{13,1},id{15,1});
inds=zeros(47,2);
inds(1,:)=[11 12];
inds(2,:)=[13 15];
count=3;
for i=1:10
   for j=i+1:10
       ids{count,1}=intersect(id{i,1},id{j,1});
       inds(count,:)=[i j];
       count=count+1;
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 3 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Panel A
figure
i=16; %Pair 'am'-'bh'
plot(meanPhicross2{i,1}(ids{i,1}),tvalcross{i,1},'ok')
hold on
plot(0:1,0:1,'--r')
xlabel('Observed dissimilarity \Phi_i^{a,b}')
ylabel('Predicted Dissimilarity')
pbaspect([1,1,1])

% Panel B
figure
for i=1:47
  x=meanPhicross2{i,1}(ids{i,1});
  pos=x>0;
  
  plot(nanmean(x(pos & isfinite(tvalcross{i,1}))),nanmean(tvalcross{i,1}(isfinite(tvalcross{i,1}) & pos)),'ok')
  hold on
  plot(nanmean(x(pos & isfinite(tvalcrossK{i,1}))),nanmean(tvalcrossK{i,1}(isfinite(tvalcrossK{i,1})& pos)),'^r') 
  plot(nanmean(x(pos & isfinite(tvalcrossS{i,1}))),nanmean(tvalcrossS{i,1}(isfinite(tvalcrossS{i,1})& pos)),'dg') 
  plot(nanmean(x(pos & isfinite(tvalcrossKS{i,1}))),nanmean(tvalcrossKS{i,1}(isfinite(tvalcrossKS{i,1})& pos)),'sb') 
end
plot(0:0.1:1,0:0.1:1,'--k' )
axis([0.33 0.5 0 0.5])
xlabel('Observed Dissimilairity <\Phi_i^{a,b}>')
ylabel('Predicted Dissimilarity')
pbaspect([1,1,1])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fig S11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
for i=1:47
  x=meanPhicross2{i,1}(ids{i,1});
  pos= x>0;
  y=tvalcross{i,1}./x;
  plot(nanmean(x(pos)),nanmedian(y(isfinite(y) & pos)),'ok')
  hold on
  y=tvalcrossK{i,1}./x;
  plot(nanmean(x(pos)),nanmedian(y(isfinite(y) & pos)),'^r') 
  y=tvalcrossS{i,1}./x;
  plot(nanmean(x(pos)),nanmedian(y(isfinite(y) & pos)),'dg') 
  y=tvalcrossKS{i,1}./x;
  plot(nanmean(x(pos)),nanmedian(y(isfinite(y) & pos)),'sb')  
end
pbaspect([1.6,1,1])
plot(0:0.1:1,ones(11,1),'--k' )
axis([0.3 0.5 0.1 1.5])
xlabel('<\Phi_i^{a,b}>')
ylabel('Median of f(K_1,K_2, \sigma_1,\sigma_2)/<\Phi_i^{a,b}>')