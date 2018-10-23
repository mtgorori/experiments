%% �T�v
%EMCL����16 %, IMCL����2 %, EMCL�N���X�^���R�̃��f���ɑ΂���
%���ω����l��bent-ray�Ɋ�Â��ĎZ�o���邽�߂̃v���O�����@2018/07/25
%% ������
close all;
clear;
clc;
%% �f�[�^�Ăяo��
cd('H:\result\2018_04_27_analyzeLipidModel')
myfilename = sprintf('2018_04_27_TOFdata_Group1ERate16IRate02Num3-2');
load(myfilename)
%% �Z���T�ݒu
param = makeParam( 0.125,0.125,400,400,0 );
param.source.point_map = cast(linspace(1,100,100),'int8');
t_size = 40.e-3;%�g�����X�f���[�T�̃T�C�Y
t_num = 200;%�g�����X�f���[�T��
t_pos = zeros(2, t_num);%�Z���T�ʒu
t_pos(1,1:t_num/2) = -t_size/2:t_size/(t_num/2-1):t_size/2 ;%�f�q������������[m]
t_pos(2,1:t_num/2) = t_size/2;
t_pos(1,t_num/2+1:t_num) = t_pos(1,1:t_num/2);
t_pos(2,t_num/2+1:t_num) = -t_size/2;
diff_t = t_pos(1,2) - t_pos(1,1);
sensor.mask=t_pos;
%% �Z���ݒ�
cell_num = 1024;
p_size = 50.e-3;%�v�Z�̈�̒���
cell_size = p_size / (cell_num-1);
x_grid = -p_size/2 : cell_size : p_size/2;
y_grid = -p_size/2 : cell_size : p_size/2;
[X, Y] = meshgrid(x_grid, y_grid);%�\���p
ds = cell_size/2;%�ʒ��i�����j�̔����ω���
t_angle = atan2(t_pos(2,:), t_pos(1,:));%�Z���T�ʒu�p
%% �}���ݒ�
% �������z
%medium.sound_speed�̔z��T�C�Y�����������߂Ƀ��T�C�Y�������s��
% v_dist= resize(medium.sound_speed',[cell_num,cell_num]);%�\���p
rateIMCL = 2;
v_muscle = 1580;%�����O�̉����l
v_fat = 1450;%�����O�̉����l�D"makeRandomMedium"�Q��
baseSoundSpeed = v_muscle + (v_fat - v_muscle)*rateIMCL/100;
radius = 6.51*1e-3 / cell_size;%[grid points]
disc1 = makeDisc(cell_num, cell_num, cell_num*121/400,cell_num*137/400,radius);
disc2 = makeDisc(cell_num, cell_num, cell_num*105/400,cell_num*261/400,radius);
disc3 = makeDisc(cell_num, cell_num, cell_num*262/400,cell_num*153/400,radius);
discs = disc1 + disc2 + disc3;
v_dist = not(discs) * baseSoundSpeed + discs * v_fat;
% ���ܗ����z
n = v_dist./baseSoundSpeed;%�\���p
% h = ones(2,2)/4;
% n = filter2(h,n);
n_cal = n';%�v�Z�p
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
%% �ő��o�H���Ƃ�f�q�y�A���o
[Min_time, inds_fastest]= min(tof_data(101:end,:));
inds_fastest = inds_fastest+100;
%% �P���Ǝ˖@��p�����ő��o�H�̐���
pos_re = [t_pos(1,:); t_pos(2,:)];%��M�f�q�ʒu�x�N�g��
length_tr2re = zeros(1,t_num/2);
%�w�i�`�ʗp
figure;
imagesc(x_grid*1e3,y_grid*1e3,n);
hold on
plot(t_pos(1,1:t_num/2)*1e3,t_pos(2,1:t_num/2)*1e3,'r','LineWidth',3);
plot(t_pos(1,t_num/2+1:end)*1e3,t_pos(2,t_num/2+1:end)*1e3,'r','LineWidth',3);
caxis([0.9 1.1]);set(gca,'Ydir','Normal');
colorbar;
xlabel('x����[mm]')
ylabel('y����[mm]')
for ii = 1:t_num/2
    %��������
    ind_fastest = inds_fastest(1,ii);
    pos_tr = [t_pos(1,ii), t_pos(2,ii)];%���M�f�q�ʒu�x�N�g��
    num_angle = 2000;%�ˏo�p�̑���(������)�F���[�v���Ƃɑ�������
    condition_1 = (ii<=t_num/2);%���M�f�q���㕔�Ɉʒu����Ƃ�True
    condition_2 = (ii > t_num/2);%���M�f�q�������Ɉʒu����Ƃ�True
    ind_search = 1;%�o�H���݂������Ƃ���0�ƂȂ�D
    %���M�f�q�̈ʒu���v���b�g
    plot(pos_tr(1)*1e3,pos_tr(2)*1e3,'*');
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
                if (ix <=5  || ix >= cell_num-5 || jy <=5 || jy >= cell_num-5)%�������z��O�̗̈�ɐi�o�����ꍇ�C����j�����ă��[�v�������D
                    clear x y
                    break
                end
                nx = (n_cal(ix+1,jy)-n_cal(ix-1,jy))/2/cell_size;
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
                if (condition_1&&(pos_ray(2)<=-t_size/2+ds/2)) || (condition_2&&(pos_ray(2)>=t_size/2-ds/2))
                    distance2re = abs(pos_ray(1)-pos_re(1,ind_fastest));
                    [min_distance2re, ~] = min(distance2re);
                    if (min_distance2re < diff_t/2)
                        %��������p
                        ind_search = 0;
                        %�`�ʗp
                        x(num_ray_head+1) = pos_ray(1);
                        y(num_ray_head+1) = pos_ray(2);
                        plot(pos_ray(1)*1e3,pos_ray(2)*1e3,'+');
                        plot(x*1e3,y*1e3,'k');
                        %�����Œ���������o���D
                        length_tr2re(1,ii) = ds * length(x);
                    end
                    clear x y
                    break
                end
                num_ray_head = num_ray_head+1;
            end
            if ind_search==0%�S��M�f�q�ɑΉ����鉹�����������Ă���΂����Ń��[�v�𔲂���D
                break
            end
        end
        if ind_search==0%�S��M�f�q�ɑΉ����鉹�����������Ă���΂����Ń��[�v�𔲂���D
            break
        end
        num_angle = num_angle + 1000;
    end
end