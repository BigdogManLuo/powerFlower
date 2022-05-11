function K=diffQ_sita(P,Q,V,sita,G,B,n_pq,n_pv,n_bal)
    n=n_pq+n_pv+n_bal;
    K=zeros(n_pq,n-n_bal);
    for i=1:n_pq
        for j=1:n-n_bal
            if i==j
                K(i,j)=(V(i)^2)*G(i,j)-P(i);
            else
                K(i,j)=V(i)*V(j)*(G(i,j)*cos(sita(i)-sita(j))+B(i,j)*sin(sita(i)-sita(j)));
            end
        end
    end
end