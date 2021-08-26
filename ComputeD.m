%This script computes the quantities D+ and D- (Eqs. 9 and 10) needed to detect transitions
%of K

load('abd_all_clean.mat') %load empirical abundances
%the variable 'abd' contains one array for each individual, in the
%following order: 10 individuals of BIO-ML, M3, F4, A (pre-travel), B
%(pre-Salmonella), A(post-travel), B(post-Salmonella)
%Array rows correspond to samples. The first column contains the sampling
%day, the following contain the counts for each OTU, the last columnn
%contains the unassigned counts. 


load('meanPhi_alldataset.mat')
load('meanPhicross_alldataset.mat')

tresh_occup=0;  %No threshold on occupancy
tresh_abd=10^(-4); %threshold on abundance
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

%Estimate K in rolling windows of size w
w=50; % window size
mindays=5; %minimum number of sampled days in a window to estimate of K
Kw1=cell(16,1);
Kw2=cell(16,1);
for i=1:16
   days=abd{i,1}(:,1); 
   Kw1{i,1}=nan(days(end)-days(1)+1,length(id{i,1}));
   Kw2{i,1}=nan(days(end)-days(1)+1,length(id{i,1}));
   for t=days(1):days(end)
       window1=t:t+w; 
       ids1=ismember(days, window1); %sampled days falling in the window
       %forward
       if sum(ids1)>=mindays
           a=abd{i,1}(ids1, 2:end);
           c=sum(a,2); %counts
           rela=a(:,id{i,1})./c; %relative abundance
           media=mean(rela);    %mean relative abundance
           mediasq=nanmean(a(:,id{i,1}).*(a(:,id{i,1})-1)./(c.*(c-1)));
           Kw1{i,1}(t-days(1)+1,:)=2*media./(2-sigma{i,1}(id{i,1}));
       else
           Kw1{i,1}(t-days(1)+1,:)=nan(1,size(id{i,1},2)); %if there are less than minday days in the window, we put nan as the K estimate
       end
       
       %backward
       window2=t-w:t;
       ids2=ismember(days, window2); %sampled days falling in the window
       if sum(ids2)>=mindays
           a=abd{i,1}(ids2, 2:end);
           c=sum(a,2); %counts
           rela=a(:,id{i,1})./c; %relative abundance
           media=mean(rela);    %mean relative abundance
           mediasq=nanmean(a(:,id{i,1}).*(a(:,id{i,1})-1)./(c.*(c-1)));
           Kw2{i,1}(t-days(1)+1,:)=2*media./(2-sigma{i,1}(id{i,1}));
       else
           Kw2{i,1}(t-days(1)+1,:)=nan(1,size(id{i,1},2));%if there are less than minday days in the window, we put nan as the K estimate
       end
   end
   
   %For windows with zero count, for which K=0 has been estimated, we
   %correct the estimate of K as the largest such that no individuals are
   %sampled in that window (see methods)
   s=sigma{i,1}(id{i,1});
   for t=1:days(end)-days(1)+1
       otus=find(Kw1{i,1}(t,:)==0);
       window1=t+days(1)-1:t+w+days(1)-1;
       ids1=ismember(days, window1);
       a=abd{i,1}(ids1, 2:end);
       c=mean(sum(a,2));
       
       for l=1:length(otus)
           si=s(otus(l));
           if ~isnan(si)
           syms k
           temp=vpasolve(1-(c-2/k/(si-2))^(1-2/si)*(k-k*si/2)^(1-2/si)==1/w,10^(-7));
           Kw1{i,1}(t,otus(l))=double(temp);
           end
       end
       
       otus=find(Kw2{i,1}(t,:)==0);
       window2=t-w+days(1)-1:t+days(1)-1;
       ids2=ismember(days, window2);
       a=abd{i,1}(ids2, 2:end);
       c=mean(sum(a,2));
       for l=1:length(otus)
           si=s(otus(l));
           if ~isnan(si)
               syms k
               temp=vpasolve(1-(c-2/k/(si-2))^(1-2/si)*(k-k*si/2)^(1-2/si)==1/w,10^(-7));
               Kw2{i,1}(t,otus(l))=double(temp);
           end
       end
   end
end


%compute D- and D+
Dm=cell(16,1);
Dp=cell(16,1);
for i=1:16
    days=abd{i,1}(:,1); 
    r=(2-sigma{i,1}(id{i,1}))./sigma{i,1}(id{i,1});
    for t=days(1):days(end)
        if any(~isnan(Kw1{i,1}(t-days(1)+1,:))) && any(~isnan(Kw2{i,1}(t-days(1)+1,:)))
            %forward
            window1=t:t+w;
            ids1=ismember(days, window1);
            a=abd{i,1}(ids1, 2:end);
            c=sum(a,2);
            a=a(:,id{i,1});
            th1=Kw1{i,1}(t-days(1)+1,:).*(sigma{i,1}(id{i,1}))/2; %theta parameter for right window
            th2=Kw2{i,1}(t-days(1)+1,:).*(sigma{i,1}(id{i,1}))/2; %theta parameter for backward window
            for j=1:sum(ids1)
                L1(j,:)=nbinpdf(a(j,:),r,1./(c(j)*th1+1)); %likelihood with right parametes
                L2(j,:)=nbinpdf(a(j,:),r,1./(c(j)*th2+1)); %likelihood with wrong parameters
            end
            L1(L1==0)=nan;
            L2(L2==0)=nan;
            Dm{i,1}(t-days(1)+1,:)=nansum(log(L1))-nansum(log(L2));
            clear L1 L2
            %backward
            window2=t-w:t;
            ids2=ismember(days, window2);
            a=abd{i,1}(ids2, 2:end);
            c=sum(a,2);
            a=a(:,id{i,1});
            for j=1:sum(ids2)
                L1(j,:)=nbinpdf(a(j,:),r,1./(c(j)*th1+1)); %likelihood with wrong parameters
                L2(j,:)=nbinpdf(a(j,:),r,1./(c(j)*th2+1));  %likelihood with right parametes
            end
            L1(L1==0)=nan;
            L2(L2==0)=nan;
            Dp{i,1}(t-days(1)+1,:)=nansum(log(L2))-nansum(log(L1));
            clear L1 L2
        else
            Dp{i,1}(t-days(1)+1,:)=nan(1,size(id{i,1},2));
            Dm{i,1}(t-days(1)+1,:)=nan(1,size(id{i,1},2));
        end
    
    end
    
end


