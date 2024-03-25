

%% PLOT FIRST COMPONENT VELOCITY (only vertex DoFs)
figure()
for k=1:size(vertices,1)
    hold on;
    index=(vertices(k,1:3))';
    uu_tmp=[uh1(index);uh1(index(1))];
    vert_temp=[xv(index),yv(index); xv(index(1)),yv(index(1))];
    fill3(vert_temp(:,1),vert_temp(:,2),uu_tmp,uu_tmp); 
end
title('Prima componente velocita con ',name)
view(3)
grid on
colorbar
hold off   
% 
% %% PLOT SECOND COMPONENT VELOCITY (only vertex DoFs)
% figure(2)
% for k=1:size(vertices,1)
%     hold on;
%     index=(vertices(k,1:3))';
%     uu_tmp=[uh2(index);uh2(index(1))];
%     vert_temp=[xv(index),yv(index); xv(index(1)),yv(index(1))];
%     fill3(vert_temp(:,1),vert_temp(:,2),uu_tmp,uu_tmp); 
% end
% view(3)
% grid on
% colorbar
% hold off 

%% PLOT QUIVER VELOCITY (only vertex DoFs)
figure()
title('Velocità quiver con ', name)
quiver(xv,yv,uh1(1:nver),uh2(1:nver));

%% PLOT PRESSURE
figure()
title('Pressioni con ',name)
for k=1:size(vertices,1)
    hold on;
    index=(vertices(k,1:3))';
    pp_tmp=[ph(index);ph(index(1))];
    vert_temp=[xv(index),yv(index); xv(index(1)),yv(index(1))];
    fill3(vert_temp(:,1),vert_temp(:,2),pp_tmp,pp_tmp); 
end
view(12,60)
grid on
colorbar
hold off