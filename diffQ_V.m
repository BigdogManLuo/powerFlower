function L=diffQ_V(P,Q,V,sita,G,B,n_pq,n_pv,n_bal)
    L=zeros(n_pq,n_pq);
    for i=1:n_pq
        for j=1:n_pq
           if i==j
               L(i,j)=(V(i)^2)*B(i,j)-Q(i);
           else
               L(i,j)=-V(i)*V(j)*(G(i,j)*sin(sita(i)-sita(j))-B(i,j)*cos(sita(i)-sita(j)));
           end
        end
    end
end