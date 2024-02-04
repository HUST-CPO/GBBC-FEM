function  plotEz(node,elem,Ez,alpha)
% %绘制Ez
figure(2);
h2=trisurf(elem,node(:,1),node(:,2),real(Ez));
shading interp
set(h2,'edgecolor','black','edgealpha',alpha);
view(2); axis equal; axis tight; axis off;
colorbar;
colormap jet;
caxis([min(real(Ez)) max(real(Ez))]);

end