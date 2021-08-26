function phi=computephi_i(n1,n2)
   n1=equalize(n1,n2);
   n2=equalize(n2,n1);
   d = n1-n2; 
   s = n1+n2;
   d=d(1:end-1); %leave out unassigned
   s=s(1:end-1); %leave out unassigned
   phi=(d.^2-s)./(s.^2-s); 
end