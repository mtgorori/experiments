%% ������
close all;
clear;
clc;
%% �Z���T�ݒu
t_size = 100.e-3;%�g�����X�f���[�T�̃T�C�Y
t_num = 256;%�g�����X�f���[�T��
t_pos = zeros(2, t_num);%�Z���T�ʒu�x�N�g��
t_pos(1,1:t_num/2) = -t_size/2:t_size/(t_num/2-1):t_size/2 ;%�f�q������������[m]
t_pos(2,1:t_num/2) = t_size/2;
t_pos(1,t_num/2+1:t_num) = t_pos(1,1:t_num/2);
t_pos(2,t_num/2+1:t_num) = -t_size/2;
%% �Z���ݒ�
cell_num = 1024;
p_size = 150.e-3;%�v�Z�̈�̒���
cell_size = p_size / (cell_num-1);
x_grid = -p_size/2 : cell_size : p_size/2;
y_grid = -p_size/2 : cell_size : p_size/2;
[X, Y] = meshgrid(x_grid, y_grid);%�\���p
ds = cell_size/4;%�ʒ��i�����j�̔����ω���
t_angle = atan2(t_pos(2,:), t_pos(1,:));%�Z���T�ʒu�p
%% �}���ݒ�
% �������z
v_water = 1500;
v_fat = 1480;
w_fat = (-50.e-3<Y) & (Y<0);
w_water = not(w_fat);
v_dist = v_water.*w_water + v_fat.*w_fat;%�\���p
% ���ܗ����z
n = v_water./v_dist;%�\���p
n_cal = n';%�v�Z�p
% 3x3�̕��ϒl�t�B���^�[�������X���]�V���O
h = ones(3,3)*1/9;
n_cal = filter2(h,n_cal);
% �}���\��1(�������z)
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
xlabel('x����[mm]')
ylabel('y����[mm]')
% exportfig('H:\result\2018_06_25_raytracing_validation\2018_06_25_sound_velocity_2layers','png',[300,300]);
% �}���\��2(���ܗ����z)
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
c.Label.String = '���΋��ܗ�[���̋��ܗ�=1.0]';
set(gca,'YDir','normal');
xlabel('x����[mm]')
ylabel('y����[mm]')
% exportfig('H:\result\2018_06_25_raytracing_validation\2018_06_25_relative_refractive_index_2layers','png',[300,300]);
%% �P���Ǝ˖@��p�����ő��o�H�̐���
for ii = 1
    pos_tr = [t_pos(1,ii), t_pos(2,ii)];%���M�f�q�ʒu�x�N�g��
    pos_ray = pos_tr; %�����ʒu�x�N�g��
    num_ray_head = 1; %�����擪�X�V��
    theta0 = pi*(3/2+0.25); %�����ˏo�p[rad]
    while(1)%�����쐬���[�v
        x(num_ray_head) = pos_ray(1);
        y(num_ray_head) = pos_ray(2);
        if num_ray_head>1 %���̏ꍇ�������s��Ȃ��Ɖ��������x�N�g�����X�V����Ȃ��D
            %dx,dy : the change of x and y
            dx = x(num_ray_head)-x(num_ray_head-1);
            dy = y(num_ray_head)-y(num_ray_head-1);
        else
            dx = ds*cos(theta0);
            dy = ds*sin(theta0);
        end
        ix = round((x(num_ray_head)+p_size/2)/cell_size+1);%���[�v���Ƃɕω����Ă���D�؂�グ���s���Ă���D
        jy = round((y(num_ray_head)+p_size/2)/cell_size+1);%�����\�z���[�v�̊e�X�e�b�v�ɂ����鉹����̓_�������O���b�h�ԍ�
        nx = (n_cal(ix+1,jy)-n_cal(ix-1,jy))/2/cell_size;%�������X�g�b�v����g���K�[�͑��ɂ���͂��D�R�����g�A�E�g
        ny = (n_cal(ix,jy+1)-n_cal(ix,jy-1))/2/cell_size;%nx,ny : the partial difference of n
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
        %���E�̐ݒ�i�v�Z�̏I�������j
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
    plot(x*1e3,y*1e3,'k')
    hold off
    axis equal;
    axis tight;
    colorbar;
    c = colorbar;
    c.Label.String = '���΋��ܗ�[���̋��ܗ�=1.0]';
    xlabel('x����[mm]')
    ylabel('y����[mm]')
    exportfig('H:\result\2018_06_25_raytracing_validation\2018_06_25_sound_ray_2layers_2','png',[300,300]);
end