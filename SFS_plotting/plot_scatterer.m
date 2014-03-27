function plot_scatterer(xq,R)

hold on;
for idx=length(R)
  rectangle('Position',[xq(idx,1)-R(idx),xq(idx,2)-R(idx), 2*R(idx), 2*R(idx)],...
  'Curvature',[1 1], ...
  'Linewidth',2);
end
hold off;
