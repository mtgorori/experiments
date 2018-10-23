function [us,RMSE]=SART_reconstruction(projection)
tic
%%%% �����p�����[�^ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CT�`��
ring_cx = 0; ring_cy = 0;% ���S�ʒu
ring_radius = 50.e-3;% ���a
rr =  ring_radius;
%����
v_water = 1540;%���̉���[m/s]
v_fat = 1420;%���b�̉���[m/s]
%���b�̈�`��
fat_cx = 20.e-3; fat_cy = 0; % ���S�ʒu
fat_radius = 10.e-3;
fd = fat_radius;
%�Z���T�ݒu
t_num = 256;%�g�����X�f���[�T��
t = linspace(0, ((t_num-1)/t_num)*2*pi, t_num);%�Z���T�ʒu�p
t_pos = zeros(2, t_num);%�Z���T�ʒu
t_pos(1,:) = rr*cos(t)+ring_cx;
t_pos(2,:) = rr*sin(t)+ring_cy;
%�Z���ݒ�
cell_num = 256;
cell_size = ring_radius*2 / cell_num;
x_grid = -ring_radius: cell_size : ring_radius;
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
%%%% �����č\�� %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pjt_est = zeros(t_num,t_num);
iteration = 30;
p = zeros(2,length(x_grid)+length(y_grid));
reI = zeros(length(x_grid),length(y_grid));
reI_store = zeros(length(x_grid),length(y_grid),iteration);
us = zeros(length(x_grid),length(y_grid),iteration);
det_I = zeros(length(x_grid),length(y_grid));
RMSE = zeros(1,iteration);

for it = 1:iteration
    count = ones(length(x_grid));
    for ii = 1:t_num
        %���M�f�q�̍��W�ݒ�(x_tr,y_tr)
        x_tr = t_pos(1,ii);
        y_tr = t_pos(2,ii);
        pos_tr = [x_tr, y_tr];
        for jj = 1:t_num
            if ii == jj
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
                rm_ind = sqrt(p(1,:).^2 + p(2,:).^2)>rr;
                p(:,rm_ind) = [];
                %�������W�ݒ�(x�̒l���ŏ��ȓ_�j
                if max(p(1,:)) - min(p(1,:)) < 3*cell_size % argx(min(p))������������ł��Ȃ����Ƃ�h�����߁D
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
                        p_init = p_neighbor;
                        zz(y_ind,x_ind) = 1;
                    end
                end
            end
            pjt_est(jj,ii) = sum(sum(reI.*zz));%(ii,jj)��ii���s�C���Ȃ킿y�����Cjj����C���Ȃ킿x�����Ȃ̂œ]�u����(x,y)�ɒu��������K�v���������D
            pjt_act = projection(jj,ii);
            count = count + zz;
            Det = ((pjt_act - pjt_est(jj,ii))/all_length).*p_length;
            det_I = det_I +Det;
        end
    end
    reI = reI + det_I./count;
    reI_store(:,:,it) = reI;
    us(:,:,it) = cell_size./(reI+cell_size/v_water);
    det_I = zeros(length(x_grid),length(y_grid));
    RMSE(1,it) = sqrt(immse(us(:,:,it), v_dist));
    if (it>1)&&((RMSE(1,it-1)-RMSE(1,it))<0.1)
        break
    end
    
end
toc
end

