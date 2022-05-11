function [G,B]=GenerateY(Nodes,AC_I,AC_J,AC_Z,AC_y,Trans_I,Trans_J,Trans_k)
Y=zeros(length(Nodes),length(Nodes));
% 生成节点导纳矩阵
%计算自导纳
for i=1:length(Nodes) %遍历每个节点
    index1=find(AC_I==Nodes(i));
    index2=find(AC_J==Nodes(i));
    for k=1:length(index1)
        index_trans=find(Trans_I==Nodes(i));
        if isempty(index_trans)==1   %不含变压器
            Y(i,i)=Y(i,i)+1/AC_Z(index1(k))+AC_y(index1(k));
        else %如果当前支路含变压器且当前节点是在k侧
            Y(i,i)=Y(i,i)+1/(Trans_k(index_trans)^2*AC_Z(index1(k)))+AC_y(index1(k));
        end
    end
    for k=1:length(index2)
        index_trans=find(Trans_I==Nodes(i));
        if isempty(index_trans)==1   %不含变压器
            Y(i,i)=Y(i,i)+1/AC_Z(index2(k))+AC_y(index2(k));
        else %如果当前支路含变压器且当前节点是在k侧
            Y(i,i)=Y(i,i)+1/(Trans_k(index_trans)^2*AC_Z(index2(k)))+AC_y(index2(k));
        end
    end
end

%计算互导纳
for i=1:length(AC_I) %遍历所有的边
    index3=find(Nodes==AC_I(i)); %找到当前索引i下对应哪两个节点
    index4=find(Nodes==AC_J(i));
    %一个变压器对一个发电机和母线 不会出现一个节点对两条变压器的情况
    index3_transI=find(Trans_I==Nodes(index3));
    index3_transJ=find(Trans_J==Nodes(index3));
    %一个节点两侧都不含变压器支路，一定不含变压器
    if isempty(index3_transI)==1&&isempty(index3_transJ)==1
        Y(index3,index4)=-1/AC_Z(i);
        Y(index4,index3)=-1/AC_Z(i);
    elseif isempty(index3_transI)==0&&isempty(index3_transJ)==1  %如果节点在变压器I侧
        for j=1:length(index3_transI) % 寻找在J侧的节点是否有index4号
            if Trans_J(index3_transI(j))==Nodes(index4) %如果找到了
                K=Trans_k(index3_transI(j));%计算变比
                Y(index3,index4)=-1/(K*AC_Z(i));
                Y(index4,index3)=Y(index3,index4);
                break;
            end
            %如果没找到 说明此次检索的支路不含变压器
            if j==length(index3_transI)
                Y(index3,index4)=-1/(AC_Z(i));
                Y(index4,index3)=Y(index3,index4);
            end
        end
    elseif isempty(index3_transJ)==0&&isempty(index3_transI)==1
        for j=1:length(index3_transJ) % 寻找在J侧的节点是否有index4号
            if Trans_I(index3_transJ(j))==Nodes(index4) %如果找到了
                K=Trans_k(index3_transJ(j));%计算变比
                Y(index3,index4)=-1/(K*AC_Z(i));
                Y(index4,index3)=Y(index3,index4);
                break;
            end
            %如果没找到 说明此次检索的支路不含变压器
            if j==length(index3_transJ)
                Y(index3,index4)=-1/(AC_Z(i));
                Y(index4,index3)=Y(index3,index4);
            end
        end
    end   
end
Y(3,3)=Y(3,3)+0.04j;
G=real(Y);B=imag(Y);
end