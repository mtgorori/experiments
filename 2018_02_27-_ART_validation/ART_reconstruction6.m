clear
close all
tic
%%%% �����p�����[�^ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%����
v_water = 1540;%���̉���[m/s]
v_fat = 1420;%���b�̉���[m/s]
%���b�̈�`��
fat_cx = 20.e-3; fat_cy = 0; % ���S�ʒu
fat_radius = 10.e-3;
fd = fat_radius;
%�Z���T�ݒu
t_size = 100.e-3;
t_num = 128;%�g�����X�f���[�T��
t_pos = zeros(2, t_num);%�Z���T�ʒu
t_pos(1,1:t_num/2) = -t_size/2:t_size/(t_num/2-1):t_size/2 ;%�f�q������������[m]
t_pos(2,1:t_num/2) = t_size/2;
t_pos(1,t_num/2+1:t_num) = t_pos(1,1:t_num/2);
t_pos(2,t_num/2+1:t_num) = -t_size/2;
%�Z���ݒ�
cell_num = 256;
p_size = 50.e-3;
cell_size = p_size*2 / cell_num;
x_grid = -p_size: cell_size : p_size;
y_grid = x_grid;
%�����������z
v_dist = v_water*ones(length(x_grid),length(y_grid));
for i = 1:length(x_grid)
    for j = 1:length(y_grid)
        if (x_grid(i)-fat_cx)^2 + (y_grid(j)-fat_cy)^2 <= fd^2
            v_dist(j,i) = v_fat;
        end
    end
end
%%%% �������z %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% �������e�f�[�^�Ăяo��%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('../result/projection_map_03(right_shifted_20mm).mat','projection')
%%%% �����č\�� %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iteration = 30;
pjt_est = zeros(t_num,t_num);
p = zeros(2,length(x_grid)+length(y_grid));
reI = zeros(length(x_grid),length(y_grid));%�č\���摜
reI_store = zeros(length(x_grid),length(y_grid),iteration);
us = zeros(length(x_grid),length(y_grid),iteration);
MSE = zeros(1,iteration);
for it = 1:iteration
    for ii = 1:t_num
        %���M�f�q�̍��W�ݒ�(x_tr,y_tr)
        x_tr = t_pos(1,ii);
        y_tr = t_pos(2,ii);
        pos_tr = [x_tr, y_tr];
        for jj = 1:t_num
            if (ii == jj) || (1<=ii)&&(ii<=t_num/2) && (1<=jj)&&(jj<=t_num/2) || (t_num/2<ii && t_num/2<jj)
                continue
            else
                %��M�f�q�̍��W�ݒ�(x_re,y_re)
                p_length = zeros(length(x_grid),length(y_grid));
                x_re = t_pos(1,jj);
                y_re = t_pos(2,jj);
                pos_re = [x_re, y_re];
                all_length = norm(pos_tr - pos_re);
                zz = zeros(length(x_grid),length(y_grid));
                %�f�q�����钼���̕�����
                x_line = ((x_tr-x_re) / (y_tr-y_re)) * (y_grid-y_re) + x_re;
                y_line = ((y_tr-y_re) / (x_tr-x_re)) * (x_grid-x_re) + y_re;
                %�e�f�q�ʒu���O���b�h�ɓ��Ă͂߂�
                [~,x_tr_index] = min(abs(x_tr - x_grid));
                [~,y_tr_index] = min(abs(y_tr - y_grid));
                [~,x_re_index] = min(abs(x_re - x_grid));
                [~,y_re_index] = min(abs(y_re - y_grid));
                %�ϕ��J�n�ʒu�C�I���ʒu�̌���i�����`�̈�):�o�H�����钷���`
                x_start = min(x_tr_index,x_re_index);
                x_end = max(x_tr_index,x_re_index);
                y_start = min(y_tr_index,y_re_index);
                y_end = max(y_tr_index,y_re_index);
                %�v�Z�̈撆�̊i�q�Ƃ̌�_���i�[(p)
                p(1,1:length(y_start:y_end)) = x_line(y_start:y_end);
                p(2,1:length(y_start:y_end)) = y_grid(y_start:y_end);
                p(1,length(y_start:y_end)+1:length(y_start:y_end)+length(x_start:x_end)) = x_grid(x_start:x_end);
                p(2,length(y_start:y_end)+1:length(y_start:y_end)+length(x_start:x_end)) = y_line(x_start:x_end);
                p(:,length(y_start:y_end)+length(x_start:x_end)+1:length(p)) = [];
                rm_ind = find(p(1,:)>t_size/2 | p(1,:)<-t_size/2 | p(2,:)>t_size/2 | p(2,:)<-t_size/2);
                p(:,rm_ind) = [];
                 p=rmmissing(p,2);%NaN���o�����߁C�Ώ��D
                %�������W�ݒ�(x�̒l���ŏ��ȓ_�j
                if max(p(1,:)) - min(p(1,:)) < cell_size % argx(min(p))������������ł��Ȃ����Ƃ�h�����߁D
                    [~,b] = min(p(2,:));
                    p_init = p(:,b);
                    while ~isempty(p)
                        %������������Z���̌��o
                        x_ind = find(abs(p_init(1,1) - x_grid) == 0, 1);
                        y_ind = find(abs(p_init(2,1) - y_grid) == 0, 1);
                        if isempty(x_ind)
                            [~,x_ind] = min(abs(p_init(1,1) - x_grid - cell_size/2));
                        end
                        if isempty(y_ind)
                            [~,y_ind] = min(abs(p_init(2,1) - y_grid - cell_size/2));
                        end
                        %�����̂�������̒[�_�̌��o
                        p(:,b) = [];%�[�_�̏���
                        [~,b] = min(p(2,:));%��������̒[�_�̃C���f�b�N�X���o
                        p_neighbor = p(:,b);%��������̒[�_���o
                        %��^�����Z�o
                        p_length(y_ind,x_ind) = norm(p_init - p_neighbor);
%                         pjt_est(ii,jj) = pjt_est(ii,jj) + p_length(x_ind,y_ind)*reI(ii,jj);
                        p_init = p_neighbor;
                        zz(y_ind,x_ind) = 1;
                    end
                else
                    [~,b] = min(p(1,:));
                    p_init = p(:,b);
                    while ~isempty(p)
                        %������������Z���̌��o
                        x_ind = find(abs(p_init(1,1) - x_grid) == 0, 1);
                        y_ind = find(abs(p_init(2,1) - y_grid) == 0, 1);
                        if isempty(x_ind)
                            [~,x_ind] = min(abs(p_init(1,1) - x_grid - cell_size/2));
                        end
                        if isempty(y_ind)
                            [~,y_ind] = min(abs(p_init(2,1) - y_grid - cell_size/2));
                        end
                        %�����̂�������̒[�_�̌��o
                        p(:,b) = [];%�[�_�̏���
                        [~,b] = min(p(1,:));%��������̒[�_�̃C���f�b�N�X���o
                        p_neighbor = p(:,b);%��������̒[�_���o
                        %��^�����Z�o
                        p_length(y_ind,x_ind) = norm(p_init - p_neighbor);
%                         pjt_est(ii,jj) = pjt_est(ii,jj) + reI(x_ind,y_ind);
                        p_init = p_neighbor;
                        zz(y_ind,x_ind) = 1;
                    end
                end
            end
            pjt_est(jj,ii) = sum(sum(reI.*zz));%(ii,jj)��ii���s�C���Ȃ킿y�����Cjj����C���Ȃ킿x�����Ȃ̂œ]�u����(x,y)�ɒu��������K�v���������D
            pjt_act = projection(jj,ii);
            det_I = (pjt_act - pjt_est(jj,ii))/all_length;
            reI = reI + det_I.*p_length;
        end
    end
     reI_store(:,:,it) = reI;
    us(:,:,it) = cell_size./(reI+cell_size/v_water);
    MSE(1,it) = immse(us(:,:,it), v_dist);
    if (it>1)&&((MSE(1,it-1)-MSE(1,it))<0.1)
    It = it-1;    
    break
    end
end

%%%% �掿�]�� %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
plot(MSE(1,1:it));
xlabel('�X�V��')
ylabel('���ϓ��덷[m/s]')
axis tight
exportfig('../result/ART_reconstruction6_MSE','png',[250,200]);
%%%% �v���t�@�C���\���@%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
plot(x_grid*1000,us(cell_num/2, :,it))
hold on
plot(x_grid*1000,v_dist(cell_num/2,:),'r');
xlabel('x����[mm]')
ylabel('����[m/s]')
legend('�č\������','��������','Location','southwest')
hold off
exportfig('../result/ART_reconstruction6_profile','png',[300,300]);
%%%% �`��p����1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
subplot(1,2,1);
imagesc(x_grid*1000, y_grid*1000, v_dist);
hold on
scatter(t_pos(1,:)*1000,t_pos(2,:)*1000,15,'r','filled');
hold off
colorbar;
set(gca,'YDir','normal');
xlabel('x����[mm]')
ylabel('y����[mm]')
caxis([1350 1600])
axis equal;
axis tight
subplot(1,2,2)
imagesc(x_grid*1000,y_grid*1000,us(:,:,it));
colorbar;
axis equal;
axis tight;
set(gca,'YDir','normal');
xlabel('x����[mm]')
ylabel('y����[mm]')
caxis([1350 1600])
exportfig('../result/ART_reconstruction6_1','png',[800,400]);

%%%%%%%%%%%%%
% figure
% imagesc(x_grid*1000,y_grid*1000,v_dist-us(:,:,It));
% colorbar;
% axis equal;
% axis tight;
% set(gca,'YDir','normal');
% xlabel('x����[mm]')
% ylabel('y����[mm]')
% caxis([1350 1600])
%%%% �`��p����2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%��������������C���̌㔭�U����l�q�������D
figure
for i = 1:4
subplot(2,3,i);
imagesc(x_grid*1000,y_grid*1000,us(:,:,i));
colorbar;
axis equal;
axis tight;
set(gca,'YDir','normal');
title(['iteration count:',num2str(i),''])
xlabel('x����[mm]')
ylabel('y����[mm]')
caxis([1350 1600])
end
i = 15;
subplot(2,3,5);
imagesc(x_grid*1000,y_grid*1000,us(:,:,i));
colorbar;
axis equal;
axis tight;
set(gca,'YDir','normal');
title(['iteration count:',num2str(i),''])
xlabel('x����[mm]')
ylabel('y����[mm]')
caxis([1350 1600])
i = it;
subplot(2,3,6);
imagesc(x_grid*1000,y_grid*1000,us(:,:,i));
colorbar;
axis equal;
axis tight;
set(gca,'YDir','normal');
title(['iteration count:',num2str(i),''])
xlabel('x����[mm]')
ylabel('y����[mm]')
caxis([1350 1600])
exportfig('../result/ART_reconstruction6_iter','png',[1000,450]);