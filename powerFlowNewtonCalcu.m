function [V,sita,V_vis,sita_vis]=powerFlowNewtonCalcu(n,n_pq,n_pv,n_bal,V0,sita0,G,B,P_s,Q_s,eps)
    %迭代初值初始化
    V=V0;
    sita=sita0;
    P=zeros(1,length(P_s));
    Q=zeros(1,length(Q_s));
    temp1=0;temp2=0;
    deltaP=100;
    deltaQ=100;
    k=1; %计数器
    %% 牛顿迭代
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
        %计算雅可比矩阵
        H=diffP_sita(P,Q,V,sita,G,B,n_pq,n_pv,n_bal);
        N=diffP_V(P,Q,V,sita,G,B,n_pq,n_pv,n_bal);
        K=diffQ_sita(P,Q,V,sita,G,B,n_pq,n_pv,n_bal);
        L=diffQ_V(P,Q,V,sita,G,B,n_pq,n_pv,n_bal);
        J=[H,N;K,L];

        %求解线性方程组 高斯消元法
        b=[deltaP;deltaQ];
        X=Gauss_solve(-J,b);
        %更新数据
        deltaSita=X(1:n_pq+n_pv);
        deltaV=X(n_pq+n_pv+1:length(X))./V(1:n_pq);
        for i=1:length(deltaSita)             %避开广播机制
            sita(i)=sita(i)+deltaSita(i);
        end
        for j=1:length(deltaV(:,1))
            V(j)=V(j)+deltaV(j,1);
        end
        %保存迭代结果
        V_vis(k,:)=V(1,:);
        sita_vis(k,:)=sita(1,:);
        %J_vis(k)=J;
        k=k+1;
    end
    sita=rad2deg(sita);
    sita_vis=rad2deg(sita_vis);
end