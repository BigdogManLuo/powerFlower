function [V,sita,V_vis,sita_vis]=powerFlowPQCalcu(n,n_pq,n_pv,n_bal,V0,sita0,G,B,P_s,Q_s,eps)
   %% 迭代初值初始化
    V=V0;
    sita=sita0;
    P=zeros(1,length(P_s));
    Q=zeros(1,length(Q_s));
    temp1=0;temp2=0;
    deltaP=100;
    deltaQ=100;
    eps=1e-3;
    B1=B(1:n-n_bal,1:n-n_bal);
    B2=B(1:n_pq,1:n_pq);
    k=1;%迭代计数器
    
    %% LU分解
    [L1,U1]=lu(B1);
    [L2,U2]=lu(B2);
    %% 迭代
    while max(abs([deltaP;deltaQ]))>eps
        %计算潮流偏移量
        for i=1:n
            for j=1:length(V)
                temp1=temp1+V(j)*(G(i,j)*cos(sita(i)-sita(j))+B(i,j)*sin(sita(i)-sita(j)));
            end
            P(i)=V(i)*temp1;
            temp1=0;
        end
        for i=1:n
            for j=1:length(V)
                temp2=temp2+V(j)*(G(i,j)*sin(sita(i)-sita(j))-B(i,j)*cos(sita(i)-sita(j)));
            end
            Q(i)=V(i)*temp2;
            temp2=0;
        end
        deltaP=P_s'-P(1:n_pq+n_pv)';
        deltaQ=Q_s'-Q(1:n_pq)';
         %迭代求解线性方程组
         c1=Gauss_solve(L1,deltaP./V(1:n-n_bal)');
         deltaSita=-Gauss_solve(U1,c1);
         c2=Gauss_solve(L2,deltaQ./V(1:n_pq)');
         deltaV=-Gauss_solve(U2,c2);
         %deltaSita=-Gauss_solve(B1,deltaP./V(1:n-n_bal)');
         %deltaV=-Gauss_solve(B2,deltaQ./V(1:n_pq)');
          %更新节点电压相角
        for i=1:length(deltaSita)             %避开广播机制
            sita(i)=sita(i)+deltaSita(i);
        end
        for j=1:length(deltaV(:,1))
            V(j)=V(j)+deltaV(j,1);
        end
        V_vis(k,:)=V(1,:);
        sita_vis(k,:)=sita(1,:);
        k=k+1;
    end
    sita=rad2deg(sita);
    sita_vis=rad2deg(sita_vis);
end