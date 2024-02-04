clc
clear all

%get mesh
model=mphload('Mesh.mph');
[~, m2]=mphmeshstats(model);

%treat mesh data
Mesh.Nodes=m2.vertex';
Elements=m2.elem{2}+1;
Mesh.Elements = sort(Elements)';
Mesh.Domains=m2.elementity{2};
Mesh.Edges=m2.elem{1}'+1;
Mesh.Edgeflag=m2.elementity{1};

save('MeshData0','Mesh');