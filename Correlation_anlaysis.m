%This script computes the correlation of estimated K and sigma parameters
%in time and across individuals, and produces figures S1 and S14

load('abd_all_clean.mat') %load empirical abundances
%the variable 'abd' contains one array for each individual, in the
%following order: 10 individuals of BIO-ML, M3, F4, A (pre-travel), B
%(pre-Salmonella), A(post-travel), B(post-Salmonella)
%Array rows correspond to samples. The first column contains the sampling
%day, the following contain the counts for each OTU, the last columnn
%contains the unassigned counts. 

%Estimation of K and sigma for each OTU of each individual
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

% Estimate K and sigma on two halves of each time-series (to compute
% time-correlations)
for i=1:16
   times=abd{i,1}(:,1);
   half=times(1)+round((times(end)-times(1))/2);
   pos=find(times>=half,1);
   abd1=abd{i,1}(1:pos, 2:end);
   abd2=abd{i,1}(pos+1:end, 2:end);
   counts1=sum(abd1,2);
   counts2=sum(abd2,2);
   relabd1=abd1(:,1:end-1)./counts1;
   relabd2=abd2(:,1:end-1)./counts2;
   meanrelabd1=mean(relabd1);    %mean relative abundance
   meanrelabd2=mean(relabd2);
   mediasq1=nanmean(abd1(:,1:end-1).*(abd1(:,1:end-1)-1)./(counts1.*(counts1-1)));
   mediasq2=nanmean(abd2(:,1:end-1).*(abd2(:,1:end-1)-1)./(counts2.*(counts2-1)));
   sigma1{i,1}=2./(1+meanrelabd1.^2./(mediasq1-meanrelabd1.^2)); %estimate on first half
   sigma2{i,1}=2./(1+meanrelabd2.^2./(mediasq2-meanrelabd2.^2));  %estimate on second half
   K1{i,1}=2*meanrelabd1./(2-sigma1{i,1});  %estimate on first half
   K2{i,1}=2*meanrelabd2./(2-sigma2{i,1});  %estimate on second half
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure S1: bias of sigma for low occupancy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
for i=1:16
    s=scatter(sigma{i,1}(sigma{i,1}>=0),occup{i,1}(sigma{i,1}>=0),'ko')
   s.MarkerEdgeAlpha = .5;
hold on
end
xlabel('\sigma')
ylabel('Occupancy')
pbaspect([1.6,1,1])
axis([-0.02 2 0 1.01])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fig S14: Example of consistency of parameters K and sigma between two 
% halves of a time-series and between two different individuals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Panels A-B
i=7
a=sigma1{i,1}(id{i,1});
b=sigma2{i,1}(id{i,1});
pos=a>0 & b>0;
x=K1{i,1}(id{i,1});
y=K2{i,1}(id{i,1});

subplot(1,2,1)
hold on
s=scatter(x(pos),y( pos),'ko')
s.MarkerEdgeAlpha = .5;
plot(min(x):10^(-5):max(x),min(x):10^(-5):max(x),'k--' )
xlabel('K first half')
ylabel('K second half')
set(gca, 'xscale','log','yscale','log')
pbaspect([1,1,1])
axis([1*10^(-6) 1 1*10^(-6) 1])

subplot(1,2,2)
hold on
s=scatter(a,b,'ko')
s.MarkerEdgeAlpha = .5;
plot(0:0.1:2,0:0.1:2,'k--' )
xlabel('\sigma first half')
ylabel('\sigma second half')
pbaspect([1,1,1])
axis([0 2 0 2])

%Panels C-D
i1=2;
i2=7;
ids=intersect(id{i1,1},id{i2,1});
a=sigma{i1,1}(ids);
b=sigma{i2,1}(ids);
x=K{i1,1}(ids);
y=K{i2,1}(ids);

subplot(1,2,1)
hold on
s=scatter(x,y,'ko')
s.MarkerEdgeAlpha = .5;
plot(10^(-6):10^(-5):1,10^(-6):10^(-5):1,'k--' )
xlabel('K "am" ')
ylabel('K "bh"')
set(gca, 'xscale','log','yscale','log')
pbaspect([1,1,1])
axis([1*10^(-6) 1 1*10^(-6) 1])

subplot(1,2,2)
hold on
s=scatter(a,b,'ko')
s.MarkerEdgeAlpha = .5;
plot(0:0.1:2,0:0.1:2,'k--' )
xlabel('\sigma "am"')
ylabel('\sigma "bh"')
pbaspect([1,1,1])
axis([0 2 0 2])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute correlations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute correlations in time, for each individual (parameters estimated
% in the two halves)
corKt=zeros(16,1);
corst=zeros(16,1);
for i=1:16
    x=K1{i,1}(id{i,1})';
    y=K2{i,1}(id{i,1})';
   
    a=sigma1{i,1}(id{i,1})';
    b=sigma2{i,1}(id{i,1})';
    ids=a<Inf & b<Inf & a>0 & b>0;
    corst(i)=corr(a(ids), b(ids));
    corKt(i)=corr(log(x(ids)), log(y(ids)));
end

%Compute correlations across pairs of individuals of the same dataset
corK=zeros(16,16);
cors=zeros(16,16);
%BIOML
for i=1:10
    for j=i+1:10
        ids=intersect(id{i,1},id{j,1});
        x=K{i,1}(ids)';
        y=K{j,1}(ids)';
        ids2=~isnan(x)&~isnan(y);
    corK(i,j)=corr(log(x(ids2)), log(y(ids2)));
    x=sigma{i,1}(ids)';
    y=sigma{j,1}(ids)';
    ids=x<Inf & y<Inf;
    cors(i,j)=corr(x(ids), y(ids));
    end
end
%MP
i=11;
j=12;
ids=intersect(id{i,1},id{j,1});
x=K{i,1}(ids)';
y=K{j,1}(ids)';
ids2=~isnan(x)&~isnan(y);
corK(i,j)=corr(log(x(ids2)), log(y(ids2)));
x=sigma{i,1}(ids)';
y=sigma{j,1}(ids)';
ids=x<Inf & y<Inf;
cors(i,j)=corr(x(ids), y(ids));
%D
i=13;
j=14;
ids=intersect(id{i,1},id{j,1});
x=K{i,1}(ids)';
y=K{j,1}(ids)';
ids2=~isnan(x)&~isnan(y);
corK(i,j)=corr(log(x(ids2)), log(y(ids2)));
x=sigma{i,1}(ids)';
y=sigma{j,1}(ids)';
ids=x<Inf & y<Inf;
cors(i,j)=corr(x(ids), y(ids));

i=14;
j=15;
ids=intersect(id{i,1},id{j,1});
x=K{i,1}(ids)';
y=K{j,1}(ids)';
ids2=~isnan(x)&~isnan(y);
corK(i,j)=corr(log(x(ids2)), log(y(ids2)));
x=sigma{i,1}(ids)';
y=sigma{j,1}(ids)';
ids=x<Inf & y<Inf;
cors(i,j)=corr(x(ids), y(ids));
       
corK=corK(corK>0); %eliminate empty elements of the matrix
cors=cors(cors>0);


save('prova.mat','corK','cors','corKt', 'corst')