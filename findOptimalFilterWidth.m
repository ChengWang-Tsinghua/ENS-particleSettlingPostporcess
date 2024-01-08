function [s,m,w] = findOptimalFilterWidth(track,Ilong)

mycolormap = mycolor('#063970','#e28743');%('#063970','#eeeee4','#e28743')
color3 = [mycolormap(1,:);mycolormap((size(mycolormap,1)+1)/2,:);mycolormap(end,:)];

%% Particle: find proper filter width
[s(1), m(1), w]=findFilterWidth_PTV(track(Ilong),'X');
[s(2), m(2), w]=findFilterWidth_PTV(track(Ilong),'Y');
[s(3), m(3), w]=findFilterWidth_PTV(track(Ilong),'Z');

%% Particle: find optimal filter length
figure;
yyaxis left
loglog(w,s(1).vx,'d-',Color=color3(1,:));hold on;
loglog(w,s(2).vx,'d-',Color=color3(2,:));
loglog(w,s(3).vx,'d-',Color=color3(3,:));
hold off

yyaxis right
loglog(w,s(1).ax,'o-',Color=color3(1,:));hold on;
loglog(w,s(2).ax,'o-',Color=color3(2,:));
loglog(w,s(3).ax,'o-',Color=color3(3,:));

legend('$u_x$','$u_y$','$u_z$','$a_x$','$a_y$','$a_z$');
% title('$\sigma(w)$')
xlabel('$w$')
yyaxis left
ylabel('$\sigma_{v}$')
yyaxis right
ylabel('$\sigma_{a}$')
grid on
axis padded