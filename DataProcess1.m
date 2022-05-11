function [Nodes,P_s,Q_s,V0,theta0,n,n_pq,n_pv,n_bal]=DataProcess1(Nodes_old,NodeType_old,Load,LD_P,LD_Q,Gen,GE_P,GE_V,GE_theta)

%节点数
n=length(Nodes_old);%总节点数
%考虑转移节点数量
n_pq=length(Load);%PQ节点数
n_pv=length(Gen)-1;%PV节点数
n_bal=1;%平衡节点数

%预赋值
V0=ones(1,n);
theta0=zeros(1,n);
P_s=[];
Q_s=[];

%节点顺序替换 按PQ PV SLACK 重新整理节点表格
index_PQ=find(NodeType_old=="PQ");
index_PV=find(NodeType_old=="PV");
index_Slack=find(NodeType_old=="Slack");
Nodes=[Nodes_old(index_PQ),Nodes_old(index_PV),Nodes_old(index_Slack)];

%各个节点的注入功率
for i=1:n_pq  %遍历新的节点表格Nodes  赋值所有的PQ节点
        index1=find(Load==Nodes(i));
        P_s=[P_s,LD_P(index1)];
        Q_s=[Q_s,LD_Q(index1)];
end

for i=n_pq+1:n_pq+n_pv %赋值PV节点
    index=find(Gen==Nodes(i));
    P_s=[P_s,GE_P(index)];
    V0(i)=GE_V(index);
end
index_temp=find(Gen==Nodes(length(Nodes))); %经过重排序后，Nodes的最后一个节点为SlACK
theta0(length(Nodes))=GE_theta;
V0(length(Nodes))=GE_V(index_temp);