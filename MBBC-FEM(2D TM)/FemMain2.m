clc
clear all
%k=[0.35*k0 0.35*k0] chi=-1

%read mesh data
load MeshData0.mat
Mesh.NbrNodes=length(Mesh.Nodes);

%treat mesh
%mirror bloch boundary
SrcBoundaryFlag=1;
DstBoundaryFlag=3;
%src
index=find(Mesh.Edgeflag==SrcBoundaryFlag(1));
index=sort(index);
SrcBoundary=Mesh.Edges(index,:);
Mesh.SrcBoundary=sort(SrcBoundary,2);
Mesh.SrcIndex=unique(sort([Mesh.SrcBoundary(:,1);Mesh.SrcBoundary(:,2)]));
%dst
index=find(Mesh.Edgeflag==DstBoundaryFlag(1));
index=sort(index);
DstBoundary=Mesh.Edges(index,:);
Mesh.DstBoundary=sort(DstBoundary,2);
Mesh.DstIndex=unique(sort([Mesh.DstBoundary(:,1);Mesh.DstBoundary(:,2)]));
%PEC
PECBoundaryFlag=[2 4 5 6];
index=find(Mesh.Edgeflag==PECBoundaryFlag(1));
for i=2:length(PECBoundaryFlag)
    index=[index;find(Mesh.Edgeflag==PECBoundaryFlag(i))];
end
index=sort(index);
PECBoundary=Mesh.Edges(index,:);
Mesh.PECBoundary=sort(PECBoundary,2);
Mesh.PECIndex=unique(sort([Mesh.PECBoundary(:,1);Mesh.PECBoundary(:,2)]));
clearvars SrcBoundary DstBoundary PECBoundary index i

%physics
physic.c_const=299792458;
physic.lam0=1;
physic.k0=pi/physic.lam0;
%material
physic.epsilonr=[1 1.5*1.5]';
physic.mur=[1 1]';

numberSolve=4;
physic.kx=physic.k0*0.35;
physic.ky=physic.k0*0.35;
%matrix assembly
[A,B,P,Mesh]=FemAssemble(Mesh,physic);

%solution
r=symrcm(A);
A=A(r,r);
B=B(r,r);
targetk0=1.5*physic.k0;
sigma=targetk0*targetk0;
[x,GAMA_sq]=eigs(A,B,numberSolve,sigma);
solverk0=diag(sqrt(GAMA_sq));
solverlam0=2*pi./solverk0;
disp("eigenfrequency:")
solverff0=physic.c_const./solverlam0
q=zeros(length(r),1);
for i=1:length(r)
    q(r(i))=i;
end
for i=1:numberSolve
    x(:,i)=x(q,i);
end
x=P*x;

%plot
Ez=x;
nn=3;%numble of mode
plotEz(Mesh.Nodes,Mesh.Elements,abs(Ez(:,nn)),0.2);%plot field
plotTri(Mesh.Nodes,Mesh.Elements);%plot mesh