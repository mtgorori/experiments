%% 初期化
close all;
clear;
clc;
%% センサ設置
t_size = 100.e-3;%トランスデューサのサイズ
t_num = 256;%トランスデューサ数
t_pos = zeros(2, t_num);%センサ位置ベクトル
t_pos(1,1:t_num/2) = -t_size/2:t_size/(t_num/2-1):t_size/2 ;%素子水平方向距離[m]
t_pos(2,1:t_num/2) = t_size/2;
t_pos(1,t_num/2+1:t_num) = t_pos(1,1:t_num/2);
t_pos(2,t_num/2+1:t_num) = -t_size/2;
%% セル設定
cell_num = 1024;
p_size = 150.e-3;%計算領域の長さ
cell_size = p_size / (cell_num-1);
x_grid = -p_size/2 : cell_size : p_size/2;
y_grid = -p_size/2 : cell_size : p_size/2;
[X, Y] = meshgrid(x_grid, y_grid);%表示用
ds = cell_size/2;%弧長（音線）の微小変化分:分母に2,4,8を入れて検討した．2018/07/23
t_angle = atan2(t_pos(2,:), t_pos(1,:));%センサ位置角
%% 媒質設定
% 音速分布
v_water = 1500;
v_fat = 1480;
w_fat = (-50.e-3<Y) & (Y<0);
w_water = not(w_fat);
v_dist = v_water.*w_water + v_fat.*w_fat;%表示用
% 屈折率分布
n = v_water./v_dist;%表示用
n_cal = n';%計算用
% 3x3の平均値フィルターをかけスム‐シング ←ここをはずす
% h = ones(3,3)*1/9;
% n_cal = filter2(h,n_cal);
% 媒質表示1(音速分布)
figure;
imagesc(x_grid*1e3,y_grid*1e3,v_dist);
hold on
plot(t_pos(1,1:t_num/2)*1000,t_pos(2,1:t_num/2)*1000,'r','LineWidth',3);
plot(t_pos(1,t_num/2+1:end)*1000,t_pos(2,t_num/2+1:end)*1000,'r','LineWidth',3);
hold off
axis equal;
axis tight;
colorbar;
c = colorbar;
c.Label.String = '[m/s]';
set(gca,'YDir','normal');
xlabel('x方向[mm]')
ylabel('y方向[mm]')
% exportfig('H:\result\2018_06_25_raytracing_validation\2. 平均化フィルタなし; 平行層\2018_07_23_sound_velocity_2layers','png',[300,300]);
% 媒質表示2(屈折率分布)
figure;
imagesc(x_grid*1e3,y_grid*1e3,n);
hold on
plot(t_pos(1,1:t_num/2)*1000,t_pos(2,1:t_num/2)*1000,'r','LineWidth',3);
plot(t_pos(1,t_num/2+1:end)*1000,t_pos(2,t_num/2+1:end)*1000,'r','LineWidth',3);
hold off
axis equal;
axis tight;
colorbar;
c = colorbar;
c.Label.String = '相対屈折率[水の屈折率=1.0]';
set(gca,'YDir','normal');
xlabel('x方向[mm]')
ylabel('y方向[mm]')
% exportfig('H:\result\2018_06_25_raytracing_validation\2. 平均化フィルタなし; 平行層\2018_07_23_relative_refractive_index_2layers','png',[300,300]);
%% 単純照射法を用いた最速経路の推定
for ii = 1
    pos_tr = [t_pos(1,ii), t_pos(2,ii)];%送信素子位置ベクトル
    pos_ray = pos_tr; %音線位置ベクトル
    num_ray_head = 1; %音線先頭更新回数
    theta0 = pi*(3/2+0.25); %初期射出角[rad]
    while(1)%音線作成ループ
        x(num_ray_head) = pos_ray(1);
        y(num_ray_head) = pos_ray(2);
        if num_ray_head>1 %この場合分けを行わないと音線方向ベクトルが更新されない．
            %dx,dy : the change of x and y
            dx = x(num_ray_head)-x(num_ray_head-1);
            dy = y(num_ray_head)-y(num_ray_head-1);
        else
            dx = ds*cos(theta0);
            dy = ds*sin(theta0);
        end
        ix = round((x(num_ray_head)+p_size/2)/cell_size+1);%ループごとに変化している．切り上げを行っている．
        jy = round((y(num_ray_head)+p_size/2)/cell_size+1);%音線構築ループの各ステップにおける音線上の点を示すグリッド番号
        nx = (n_cal(ix+1,jy)-n_cal(ix-1,jy))/2/cell_size;%nx,ny : the partial difference of n
        ny = (n_cal(ix,jy+1)-n_cal(ix,jy-1))/2/cell_size;
        detx = (x(num_ray_head)+p_size/2)/cell_size+1-ix;
        dety = (y(num_ray_head)+p_size/2)/cell_size+1-jy;
        if detx>=0
            ix2 = ix+1;
        else
            ix2 = ix-1;
        end
        if dety>=0
            jy2 = jy+1;
        else
            jy2 = jy-1;
        end
        lx1 = abs((x(num_ray_head)+p_size/2)/cell_size+1-ix);lx2 = abs((x(num_ray_head)+p_size/2)/cell_size+1-ix2);
        ly1 = abs((y(num_ray_head)+p_size/2)/cell_size+1-jy);ly2 = abs((y(num_ray_head)+p_size/2)/cell_size+1-jy2);
        n_inter = n(ix,jy)*lx2*ly2+n(ix2,jy)*lx1*ly2+n(ix,jy2)*lx2*ly1+n(ix2,jy2)*lx1*ly1;
        DS = sqrt((dx+1/2/n_inter*(nx-(nx*dx/ds)*dx/ds)*cell_size^2)^2+(dy+1/2/n_inter*(ny-(ny*dy/ds)*dy/ds)*cell_size^2)^2);
        dsx = (dx+1/2/n_inter*(nx-(nx*dx/ds)*dx/ds)*cell_size^2)/DS*ds;
        dsy = (dy+1/2/n_inter*(ny-(ny*dy/ds)*dy/ds)*cell_size^2)/DS*ds;
        pos_ray(1) = x(num_ray_head)+dsx;
        pos_ray(2) = y(num_ray_head)+dsy;
        %境界の設定（計算の終了条件）
        if pos_ray(2)<-t_size/2 || pos_ray(1)==0 || pos_ray(1)==p_size/2
            %             re_distance = sqrt((r(1)-re(1))^2+(r(2)-re(2))^2);
            break
        end
        num_ray_head = num_ray_head+1;
    end
    figure;
    imagesc(x_grid*1e3,y_grid*1e3,n);
    hold on
    plot(t_pos(1,1:t_num/2)*1e3,t_pos(2,1:t_num/2)*1e3,'r','LineWidth',3);
    plot(t_pos(1,t_num/2+1:end)*1e3,t_pos(2,t_num/2+1:end)*1e3,'r','LineWidth',3);
    set(gca,'Ydir','Normal');
    plot(pos_tr(1)*1e3,pos_tr(2)*1e3,'*');
    plot(pos_ray(1)*1e3,pos_ray(2)*1e3,'+');
    plot(x*1e3,y*1e3,'k');
    hold off
    axis equal;
    axis tight;
    colorbar;
    c = colorbar;
    c.Label.String = '相対屈折率[水の屈折率=1.0]';
    xlabel('x方向[mm]')
    ylabel('y方向[mm]')
    exportfig('H:\result\2018_06_25_raytracing_validation\2. 平均化フィルタなし; 平行層\2018_07_23_sound_ray_2layers_ds_2','png',[300,300]);
end