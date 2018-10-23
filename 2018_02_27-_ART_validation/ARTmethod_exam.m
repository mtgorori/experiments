dim = 2;
w = [-2, 1; 1, 1];
p = [2; 5];
figure;
grid on
axis equal
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
hold on
fimplicit(@(f1,f2) w(1,1)*f1 + w(1,2)*f2 - 2)
fimplicit(@(f1,f2) w(2,1)*f1 + w(2,2)*f2 - 5)
OU = w(1,:)/norm(w(1,:));
quiver(0,0,OU(1,1),OU(1,2),...
       'linewidth',1.5,...
       'MaxHeadSize',1.0)