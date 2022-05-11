function N=diffP_V(P,Q,V,sita,G,B,n_pq,n_pv,n_bal)
    n=n_pq+n_pv+n_bal;
    N=zeros(n-n_bal,n_pq);
    for i=1:n-n_bal
        for j=1:n_pq
            if i==j
                N(i,j)=-(V(i)^2)*G(i,j)-P(i);
            else
                N(i,j)=-V(i)*V(j)*(G(i,j)*cos(sita(i)-sita(j))+B(i,j)*sin(sita(i)-sita(j)));
            end
        end
    end
end