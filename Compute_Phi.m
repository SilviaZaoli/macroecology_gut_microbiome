%This script computes \Phi_i(T) for each OTU of each individual of the
%dataset, as well as \Phi_i^{a.b} between different individuals of the same
%dataset

load('abd_all_clean.mat') %load empirical abundances
%the variable 'abd' contains one array for each individual, in the
%following order: 10 individuals of BIO-ML, M3, F4, A (pre-travel), B
%(pre-Salmonella), A(post-travel), B(post-Salmonella)
%Array rows correspond to samples. The first column contains the sampling
%day, the following contain the counts for each OTU, the last columnn
%contains the unassigned counts. 


% Within individuals
Phi_i=cell(16,600);

T=cell(16,1);
for i=1:16
    times=abd{i,1}(:,1);
    span=times(end)-times(1);
    T{i,1}=1:span;
end

%For each individual k and value of T, Phi_i{k,i} contains all the
%Phi_i(t,T) that can be computed for each OTU (for different values of t)
for k=1:16
    times=abd{k,1}(:,1);
    for i=1:length(T{k,1})
        for t=1:times(end)-T{k,1}(i)
            t1= abd{k,1}(:,1)==t;
            t2= abd{k,1}(:,1)==t+T{k,1}(i);
            if any(t1) && any(t2) %both samples exist
                n1=abd{k,1}(t1,2:end);
                n2=abd{k,1}(t2,2:end);
                Phi_i{k,i}=[Phi_i{k,i}; computephi_i(n1,n2)];
            end
        end
    end
end

%To obtain Phi_i(T), we average over t
meanPhi=cell(16,1);
for k=1:16
    N_OTU=size(Phi_i{k,1},2);
    meanPhi{k,1}=zeros(N_OTU, length(T{k,1})+1);
    for i=1:length(T{k,1})+1
        if size(Phi_i{k,i},1)>5 % the value of Phi_i(T) is computed only if there are at least 5 values of Phi(t,T) to average
            meanPhi{k,1}(:,i)=nanmean(Phi_i{k,i});
        else
            meanPhi{k,1}(:,i)=nan(N_OTU,1);
        end
    end
end

save('meanPhi_alldataset.mat', 'meanPhi')

% Across individuals
pairs=[[11,12]; [13,15]; nchoosek(1:10,2)]; % All pairs of individuals of the same dataset
Tcross=[1:100 101:10:150]; %values of T
meanPhicross=cell(length(pairs),1);

for i=1:length(pairs)
    pair=pairs(i,:);
    %Put first the shortest one
    if length(T{pair(1),1})<length(T{pair(2),1})
        id1=pair(1);
        id2=pair(2);
    else
        id1=pair(2);
        id2=pair(1);
    end
    times1=abd{id1,1}(:,1)-abd{id1,1}(1,1)+1;
    times2=abd{id2,1}(:,1)-abd{id2,1}(1,1)+1;
    for tau=1:length(Tcross)
        temp=[];
        for t=1:times1(end)-Tcross(tau)
            t1=times1==t;
            t2=times2==t+Tcross(tau);
            if any(t1) && any(t2)
                n1=abd{id1,1}(t1,2:end);
                n2=abd{id2,1}(t2,2:end);
                n1=equalize(n1,n2);
                n2=equalize(n2,n1);
                temp=[temp; computephi_i(n1,n2)];
            end
        end
        meanPhicross{i,1}(:,tau)=nanmean(temp);
    end
    
end

save('meanPhicross_alldataset2.mat', 'meanPhicross')
    