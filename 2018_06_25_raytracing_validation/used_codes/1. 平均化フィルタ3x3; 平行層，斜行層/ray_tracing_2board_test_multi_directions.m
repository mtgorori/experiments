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
v_water = 1540;
v_fat = 1420;
w_fat = (-50.e-3<Y) & (Y<0);
w_water = not(w_fat);
v_dist = v_water.*w_water + v_fat.*w_fat;%�\���p
% ���ܗ����z
n = v_water./v_dist;%�\���p
n_cal = n';%�v�Z�p
% 3x3�̕��ϒl�t�B���^�[�������X���]�V���O
h = ones(5,5)*1/25;
n_cal = filter2(h,n_cal);
% �}���\��1(�������z)
figure;
imagesc(x_grid*1e3,y_grid*1e3,v_dist);
hold on
plot(t_pos(1,1:t_num/2)*1000,t_pos(2,1:t_num/2)*1000,'r','LineWidth',3);
plot(t_pos(1,t_num/2+1:end)*1000,t_pos(2,t_num/2+1:end)*1000,'r','LineWidth',3);
hold off
colorbar;
c = colorbar;
c.Label.String = '[m/s]';
set(gca,'YDir','normal');
xlabel('x����[mm]')
ylabel('y����[mm]')
% �}���\��2(���ܗ����z)
figure;
imagesc(x_grid*1e3,y_grid*1e3,n);
hold on
plot(t_pos(1,1:t_num/2)*1000,t_pos(2,1:t_num/2)*1000,'r','LineWidth',3);
plot(t_pos(1,t_num/2+1:end)*1000,t_pos(2,t_num/2+1:end)*1000,'r','LineWidth',3);
hold off
colorbar;
c = colorbar;
c.Label.String = '[m/s]';
set(gca,'YDir','normal');
xlabel('x����[mm]')
ylabel('y����[mm]')
%% �P���Ǝ˖@��p�����ő��o�H�̐���
pos_re = [t_pos(1,:); t_pos(2,:)];%��M�f�q�ʒu�x�N�g��
for ii = 1:1
    %��������
    pos_tr = [t_pos(1,ii), t_pos(2,ii)];%���M�f�q�ʒu�x�N�g��
    num_angle = 100;%�ˏo�p�̑���(������)�F���[�v���Ƃɑ�������
    condition_1 = (ii<=t_num/2);%���M�f�q���㕔�Ɉʒu����Ƃ�True
    condition_2 = (ii > t_num/2);%���M�f�q�������Ɉʒu����Ƃ�True
    ind_re = ones(1,t_num/2);%��M�f�q�y�A��������΂����v�f������0�ɒu�������
    %�w�i�`�ʗp
    figure;
    imagesc(x_grid*1e3,y_grid*1e3,n);
    hold on
    plot(t_pos(1,1:t_num/2)*1e3,t_pos(2,1:t_num/2)*1e3,'r','LineWidth',3);
    plot(t_pos(1,t_num/2+1:end)*1e3,t_pos(2,t_num/2+1:end)*1e3,'r','LineWidth',3);
    plot(pos_tr(1)*1e3,pos_tr(2)*1e3,'*');
    caxis([0.9 1.1]);set(gca,'Ydir','Normal');
    colorbar;
    xlabel('x����[mm]')
    ylabel('y����[mm]')
    
    while(1)%�����p�x���}���[�v
        initial_angle = condition_1*linspace(pi,2*pi,num_angle)+...
            condition_2*linspace(0,pi,num_angle);%�����ˏo�p[rad]�F���M�f�q�̑����镽�ɂ���Ċp�x�͈̔͂��قȂ�
        for ind_angle = 1:num_angle
            pos_ray = pos_tr; %�����ʒu�x�N�g���i�������j
            num_ray_head = 1; %�����擪�X�V��(������)
            while(1)%�����쐬���[�v
                x(num_ray_head) = pos_ray(1);
                y(num_ray_head) = pos_ray(2);
                if num_ray_head>1 %���̏ꍇ�������s��Ȃ��Ɖ��������x�N�g�����X�V����Ȃ��D
                    %dx,dy : the change of x and y
                    dx = x(num_ray_head)-x(num_ray_head-1);
                    dy = y(num_ray_head)-y(num_ray_head-1);
                else
                    dx = ds*cos(initial_angle(ind_angle));%���̕����̓��[�v�̍Ō�Ɏ����Ă���Ώ���������ȗ��ł���Ǝv���D2018/06/21
                    dy = ds*sin(initial_angle(ind_angle));
                end
                ix = round((x(num_ray_head)+p_size/2)/cell_size+1);%���[�v���Ƃɕω����Ă���D�؂�グ���s���Ă���D
                jy = round((y(num_ray_head)+p_size/2)/cell_size+1);%�����\�z���[�v�̊e�X�e�b�v�ɂ����鉹����̓_�������O���b�h�ԍ�
                if (ix <=5  || ix >= cell_num-5 || jy <=5 || jy >= cell_num-5)%�������v�Z���E�ɓ��B�����ꍇ�C����j�����ă��[�v�������D
                    % ix == 1�ȂǂłȂ��̂́C�p�f�B���O�����̌����ɑΉ����邽��
                    clear x y
                    break
                end
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
                lx1 = abs((x(num_ray_head)+p_size/2)/cell_size+1-ix);
                lx2 = abs((x(num_ray_head)+p_size/2)/cell_size+1-ix2);
                ly1 = abs((y(num_ray_head)+p_size/2)/cell_size+1-jy);
                ly2 = abs((y(num_ray_head)+p_size/2)/cell_size+1-jy2);
                n_inter = n(ix,jy)*lx2*ly2+n(ix2,jy)*lx1*ly2+n(ix,jy2)*lx2*ly1+n(ix2,jy2)*lx1*ly1;
                DS = sqrt((dx+1/2/n_inter*(nx-(nx*dx/ds)*dx/ds)*cell_size^2)^2+...
                    (dy+1/2/n_inter*(ny-(ny*dy/ds)*dy/ds)*cell_size^2)^2);
                dsx = (dx+1/2/n_inter*(nx-(nx*dx/ds)*dx/ds)*cell_size^2)/DS*ds;
                dsy = (dy+1/2/n_inter*(ny-(ny*dy/ds)*dy/ds)*cell_size^2)/DS*ds;
                pos_ray(1) = x(num_ray_head)+dsx;
                pos_ray(2) = y(num_ray_head)+dsy;
                %���E�̐ݒ�i�v�Z�̏I�������j
                if (condition_1&&(pos_ray(2)<-t_size/2)) || (condition_2&&(pos_ray(2)>t_size/2))
                    distance2re = abs(pos_ray(1)-pos_re(1,1:t_num/2));
                    [min_distance2re, I] = min(distance2re);
                    if (min_distance2re < 3.e-4) && (ind_re(I) == 1)
                        %��������p
                        ind_re(I) = 0;
                        %�`�ʗp
                        plot(pos_ray(1)*1e3,pos_ray(2)*1e3,'+');
                        plot(x*1e3,y*1e3,'k');
                    end
                    break
                end
                num_ray_head = num_ray_head+1;
            end
        end
        if sum(ind_re)==0%�S��M�f�q�ɑΉ����鉹�����������Ă���΂����Ń��[�v�𔲂���D
            break
        end
        num_angle = num_angle + 100;
    end
end