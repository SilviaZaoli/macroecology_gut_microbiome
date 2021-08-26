%This script produces Figure S6,7 and 8

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

individuals2=individuals;
individuals2{1,13}='A bef';
individuals2{1,14}='A aft';
individuals2{1,15}='B bef';
individuals2{1,16}='B aft';

%Identify OTUs with transitions in K
load('Infoloss_5days.mat') %Load precomputed D+ and D- 
thresD=40; %threshold value of D to identify a transition
jumpingOTU=cell(16,1);  %OTUs for which either D+ or D- have a peak above threshold
Dpeaks1=cell(16,1);  %time and height of highest peak of D-
Dpeaks2=cell(16,1);  %time and height of highest peak of D+
for i=1:16
    temp=false(size(Dm{i,1},2),1);
    for j=1:size(Dm{i,1},2)
        m=find(Dm{i,1}(:,j)==max(Dm{i,1}(:,j)));
        Dpeaks1{i,1}(j,:)=[m(1),Dm{i,1}(m(1),j)];
        if any(Dm{i,1}(:,j)>thresD)
            temp(j)=true;
        end
    end
    
    for j=1:size(Dp{i,1},2)
        m=find(Dp{i,1}(:,j)==max(Dp{i,1}(:,j)));
        Dpeaks2{i,1}(j,:)=[m(1),Dp{i,1}(m(1),j)];
        if any(Dp{i,1}(:,j)>thresD)
            temp(j)=true;
        end
    end
jumpingOTU{i,1}=find(temp); 
end

%Indentify OTUs with increasing Phi by fitting linearly the normalized Phi(T) and comparing its slope to the threshold for that individual 

load('threshold_slope.mat') %Load threshold slope used to classify Phi as flat or increasing

%Normalize Phi by its theoretical asymptotic value
load('meanPhi_alldataset.mat') %Load <Phi(T)> for all OTUs of all individuals
for i=1:16
normPhi{i,1}=meanPhi{i,1}(id{i,1},:)'./(sigma{i,1}(id{i,1})./(4-sigma{i,1}(id{i,1})));
end

slope=cell(16,1);
OTUslope=cell(16,1); %ID of OTUs whose Phi(T) is increasing
for i=1:16
   points=normPhi{i,1}(11:end-1,:); %discard initial transient
   for j=1:length(id{i,1})
       lm=fitlm(1:length(points(:,j)),points(:,j));
       slope{i,1}(j)=lm.Coefficients.Estimate(2);
   end
   id_slope{i,1}=slope{i,1}>threshold(i); %Phi is identified as increasing if fitted slope is larger than threshold
   OTUslope{i,1}=find(id_slope{i,1});
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figures S6A-B: examples of abundance time-series with K jumps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=11;  %Individual M3
x=abd{i,1}(:,1)-abd{i,1}(1,1); %days, starting from 0
y=relabd{i,1}(:,id{i,1}); %abundances
ids=jumpingOTU{i,1}; 
y=y(:,ids); %abundances of OTUs with transitions in K

%Panel A
figure
plot(x, y(:,1)/max(y(:,1)),'-') %time series for one OTU, normalized by its maximum 
hold on
plot(x, y(:,18)/max(y(:,18)),'-')
plot(x, y(:,21)/max(y(:,21)),'-')
xlabel('Day')
ylabel('Relative abundance')
pbaspect([1.6 1 1])
xlim([0 100])
ylim([-0.1 1.1])

%Panel B
figure
plot(x, y(:,6)/max(y(:,6)),'-')
hold on
plot(x, y(:,16)/max(y(:,16)),'-')
plot(x, y(:,23)/max(y(:,23)),'-')
xlabel('Day')
ylabel('Relative abundance')
pbaspect([1.6 1 1])
xlim([0 400])
ylim([-0.1 1.1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure S6C: bar plot of OTUs where we detect a jump in K and whose ̃Φ is increasing, where we detect a jump in K but ̃Φ is flat
% and where ̃Φ is increasing but we detect no jump in K
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bars=zeros(16,3);
count=1;
for i=1:16
    bars(count,:)=[length(setdiff(OTUslope{i,1}, jumpingOTU{i,1})), length(intersect(OTUslope{i,1}, jumpingOTU{i,1})),length(setdiff( jumpingOTU{i,1},OTUslope{i,1}))];
    bars(count,:)=bars(count,:)/length(id{i,1})*100';
    count=count+1;
end

figure
bar(bars,'stacked')
xticks(1:16)
xticklabels(individuals2)
legend({'\Phi increasing, no K change','\Phi increasing and K change', 'K change, \Phi flat'})
pbaspect([1.6,1,1])
axis([0 17 0 100])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Table S1: Are OTUs with K transition over-represented in OTU with
%increasing Phi(T)?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:16
    M=length(id{i,1}); % #OTU in this individual
    K=length(jumpingOTU{i,1}); % #OTUs with K transistion
    N=length(OTUslope{i,1}); % #OTU with increasing Phi(T)
    x=length(intersect(jumpingOTU{i,1},OTUslope{i,1})); % #OTUs with K transition AND increasing Phi(T)
    p(i)=hygecdf(x,M,K,N,'upper')+hygepdf(x,M,K,N);  % probability to observe x or more OTUs with both characteristics (p-value of hypergeometric test) 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure S9: histograms of times at which jumps in K are detected
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
count=1;
for i=1:16
    jumptimes1=nan(length(Dpeaks1{i,1}),1); %times when D- identifies a jump
    jumptimes2=nan(length(Dpeaks1{i,1}),1); %times when D+ identifies a jump
    for j=1:length(Dpeaks1{i,1})
        if Dpeaks1{i,1}(j,2)>thresD  &&   Dpeaks2{i,1}(j,2)<thresD %if only D- has a peak above threshold
            jumptimes1(j)=Dpeaks1{i,1}(j,1);
        elseif Dpeaks1{i,1}(j,2)<thresD  &&   Dpeaks2{i,1}(j,2)>thresD %if only D+ has a peak above threshold
            jumptimes2(j)=Dpeaks2{i,1}(j,1);
        elseif Dpeaks1{i,1}(j,2)>thresD  &&   Dpeaks2{i,1}(j,2)>thresD %if both do
            if abs(Dpeaks1{i,1}(j,1)-Dpeaks2{i,1}(j,1))>10 % if the two peaks are more that 10 days apart, I record both as jumps
                jumptimes1(j)=Dpeaks1{i,1}(j,1);
                jumptimes2(j)=Dpeaks2{i,1}(j,1);
            elseif Dpeaks1{i,1}(j,2)>Dpeaks2{i,1}(j,2) % if they are closer, I record only the highest one
                jumptimes1(j)= Dpeaks1{i,1}(j,1);
            elseif Dpeaks1{i,1}(j,2)<Dpeaks2{i,1}(j,2)
                jumptimes2(j)=Dpeaks2{i,1}(j,1);
            end
        end
    end
    
    jumptimes=[jumptimes1(~isnan(jumptimes1));jumptimes2(~isnan(jumptimes2))]; %all times at which a jump is identified
    
    subplot(4,4,count)
    histogram(jumptimes,1:10:abd{i,1}(end,1)-abd{i,1}(1,1))
    pbaspect([1.6,1,1])
    xlabel('days')
    ylabel('OTU with K change')
    count=count+1;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fig. S7: Example of our method to identify jumps in the abundance time series.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=11; %individual M3
x=abd{i,1}(:,1)-abd{i,1}(1,1); %days, starting from 0
y=relabd{i,1}(:,id{i,1}); %abundances
ids=jumpingOTU{i,1}; %IDs of OTU with K jumps
y=y(:,ids);

figure
j=23;
plot(x,y(:,j));
ylabel('Relative abundance')
hold on
yyaxis right
plot(Dm{i,1}(:,ids(j)),'-k');
yyaxis right
plot(Dp{i,1}(:,ids(j)),'-r');
[a,b]=max(Dp{i,1}(:,ids(j)));
plot([b,b],[-100 500],'--k')
ylabel('D_+, D_-')
xlabel('days')
pbaspect([1.6 1 1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure S8: examples of cases where Φ(T) is increasing even though no jump in K is detected or viceversa
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i=11; %individual M3
x=abd{i,1}(:,1)-abd{i,1}(1,1); %days, starting from 0
y=relabd{i,1}(:,id{i,1}); %abundances
ids=setdiff(OTUslope{i,1},jumpingOTU{i,1}); %IDs of OTus with increasing Phi(T) but no K jumps
y=y(:,ids);

%panel A
figure
j=1;
plot(x,y(:,j));
ylabel('Relative abundance')
hold on
yyaxis right
plot(Dm{i,1}(:,ids(j)),'-k');
yyaxis right
plot(Dp{i,1}(:,ids(j)),'-r');
ylabel('D_+, D_-')
xlabel('days')
pbaspect([1.6 1 1])

%panel B
figure
j=22;
plot(x,y(:,j));
ylabel('Relative abundance')
hold on
yyaxis right
plot(Dm{i,1}(:,ids(j)),'-k');
yyaxis right
plot(Dp{i,1}(:,ids(j)),'-r');
ylabel('D_+, D_-')
xlabel('days')
pbaspect([1.6 1 1])

%Panel C
i=12; %individual F4
x=abd{i,1}(:,1)-abd{i,1}(1,1);  %days, starting from 0
y=relabd{i,1}(:,id{i,1}); %abundances
ids=setdiff(jumpingOTU{i,1},OTUslope{i,1});  %IDs of OTus with with detected K jump but flat Phi(T)
y=y(:,ids);

figure
j=3;
plot(x,y(:,j));
ylabel('Relative abundance')
hold on
yyaxis right
plot(Dm{i,1}(:,ids(j)),'-k');
yyaxis right
plot(Dp{i,1}(:,ids(j)),'-r');
ylabel('D_+, D_-')
xlabel('days')
pbaspect([1.6 1 1])

%panel D
i=11; %individual M3
x=abd{i,1}(:,1)-abd{i,1}(1,1);  %days, starting from 0
y=relabd{i,1}(:,id{i,1}); %abundances
ids=setdiff(jumpingOTU{i,1},OTUslope{i,1});  %IDs of OTus with with detected K jump but flat Phi(T)
y=y(:,ids);

figure
j=9;
plot(x,y(:,j));
ylabel('Relative abundance')
hold on
yyaxis right
plot(Dm{i,1}(:,ids(j)),'-k');
yyaxis right
plot(Dp{i,1}(:,ids(j)),'-r');
[a,b]=max(Dp{i,1}(:,ids(j)));
ylabel('D_+, D_-')
xlabel('days')
pbaspect([1.6 1 1])
