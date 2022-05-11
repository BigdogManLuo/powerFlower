function H=diffP_sita(P,Q,V,sita,G,B,n_pq,n_pv,n_bal)
    n=n_pq+n_pv+n_bal;
    H=zeros(n-n_bal,n-n_bal);
    for i=1:n-n_bal
        for j=1:n-n_bal
           if i==j
               H(i,j)=(V(i)^2)*B(i,j)+Q(i);
           else
               H(i,j)=-V(i)*V(j)*(G(i,j)*sin(sita(i)-sita(j))-B(i,j)*cos(sita(i)-sita(j)));
           end
        end
    end
end