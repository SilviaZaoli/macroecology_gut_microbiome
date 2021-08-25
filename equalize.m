function n_new=equalize(n1,n2)
  n=sum(n2);
  if sum(n1)<=n
     n_new=n1;
  else
     list=[];
     for i=1:length(n1)
        list=[list;  repmat(i,n1(i),1)]; %make list of single reads 
     end
     downsampled=randsample(list, n);
     n_new=hist(downsampled,1:length(n1));
  end
end