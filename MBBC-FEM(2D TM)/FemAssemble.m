function [A,B,P,mesh]=FemAssemble(mesh,physic)

[xt,yt,pt,IntegralOrder]=GetGuassPoints(2);

mesh.NbrElement=length(mesh.Elements);
mesh.Dof=mesh.NbrNodes;
NbrNon=mesh.NbrElement*9;
NumNonA=0;
Ai=zeros(NbrNon,1);
Aj=zeros(NbrNon,1);
Av=complex(zeros(NbrNon,1));
NumNonB=0;
Bi=zeros(NbrNon,1);
Bj=zeros(NbrNon,1);
Bv=complex(zeros(NbrNon,1));

NbrBBC=length(mesh.SrcIndex);
mesh.NbrBBC=NbrBBC;
mesh.mapping=zeros(NbrBBC,1);
for i=1:NbrBBC
    x2=mesh.Nodes(mesh.DstIndex(i),1);
    for j=1:NbrBBC
        y1=mesh.Nodes(mesh.SrcIndex(j),2);
        if abs(y1-x2)<0.005
            mesh.mapping(i)=mesh.SrcIndex(j);
            break;
        end
    end
end

for n=1:mesh.NbrElement
    x=mesh.Nodes(mesh.Elements(n,:),1);
    y=mesh.Nodes(mesh.Elements(n,:),2);

    Jac=zeros(3,3);JacS=zeros(3,3);
    Jac(1,1)=x(2)-x(1);Jac(1,2)=y(2)-y(1);
    Jac(2,1)=x(3)-x(1);Jac(2,2)=y(3)-y(1);Jac(3,3)=1;
    InvJac=inv(Jac);
    JacS(1,1)=InvJac(2,2);JacS(1,2)=-InvJac(2,1);
    JacS(2,1)=-InvJac(1,2);JacS(2,2)=InvJac(1,1);JacS(3,3)=1;
    DetJac=abs(det(Jac));
    
    Ez=zeros(3,3,IntegralOrder);curlEz=zeros(3,3,IntegralOrder);
    temp=zeros(3,1);
    for i=1:IntegralOrder
        for j=1:3
            temp(3)=BF_Ez(j,xt(i),yt(i));temp(1)=0;temp(2)=0;
            Ez(j,:,i)=temp;
            [temp(1),temp(2)]=BF_curlEz(j,xt(i),yt(i));temp(3)=0;
            curlEz(j,:,i)=JacS*temp;
        end
    end
    
    domain=mesh.Domains(n);
    epsilonr=physic.epsilonr(domain);
    
    Sz=zeros(3,3);Tz=zeros(3,3);
    for i=1:3
        for j=1:3
            for k=1:IntegralOrder
                Sz(i,j)=Sz(i,j)+pt(k)*DetJac*dot(curlEz(i,:,k),curlEz(j,:,k));
                Tz(i,j)=Tz(i,j)+pt(k)*DetJac*epsilonr*dot(Ez(i,:,k),Ez(j,:,k));
            end
        end
    end
    
    for i=1:3
        for j=1:3
            MappingIndexSi=mesh.Elements(n,i);
            MappingIndexSj=mesh.Elements(n,j);
            
            NumNonA=NumNonA+1;
            Ai(NumNonA)=MappingIndexSi;
            Aj(NumNonA)=MappingIndexSj;
            Av(NumNonA)=Sz(i,j);
            
            NumNonB=NumNonB+1;
            Bi(NumNonB)=MappingIndexSi;
            Bj(NumNonB)=MappingIndexSj;
            Bv(NumNonB)=Tz(i,j);
        end
    end
end

A = sparse(Ai,Aj,Av);
B = sparse(Bi,Bj,Bv);

physicIndex=-exp(-1i*physic.kx*1);
P=zeros(mesh.NbrNodes,mesh.NbrNodes);
j=1;
for i=1:mesh.NbrNodes
    if mesh.DstIndex(j)==i
        j=j+1;
    else
        P(i,i)=1;
    end
end
for i=1:NbrBBC
    if i==19
        P(mesh.DstIndex(i),1)=physicIndex*physicIndex;
    else
        P(mesh.DstIndex(i),mesh.mapping(i))=physicIndex;
    end
end
newIndex=unique(sort([mesh.PECIndex;mesh.DstIndex]));
origInex=1:mesh.Dof;
for i=1:length(newIndex)
    origInex(newIndex(length(newIndex)-i+1))=[];
end
P=P(:,origInex);

A=P'*A*P;
B=P'*B*P;

end