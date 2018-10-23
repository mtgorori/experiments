tic
%make_layered_medium
%�T�v�����F�ؓ��C���b�C�ؓ��̎O�w�\��(�g�����X�f���[�T�ɑ΂��Đ��������j�D���b�w�̌��������`�ɑ�������ƁCROI���̕��ω��������`�ɑ�������D
%�ړI�F���̃f�[�^����ɂ��ĉ������z�č\���≹��������s���āC�^�l�Ƃ̍���]������D
%   �ڍא����������ɋL�q
%�@input:Layer.thickness�[���b�w�̌���[mm]
%�@output:Medium.sound_speed�[�}���̉������z[m/s]
%               Medium.density�[�}���̖��x���z[kg/m3]
%               Layer.sound_speed_ave�[ROI���̕��ω���[m/s]
% ����
% v_muscle = 1580;%�ؓ��̉���[m/s]
% v_fat = 1450;%���b�̉���[m/s]
% ���x
% den_muscle = 1040;%�ؓ��̖��x[kg/m3]
% den_fat = 920;%���b�̖��x[kg/m3]
clear 
close all
% �����p�����[�^
Sensor.sizeTotal = 100.e-3;%�g�����X�f���[�T��[m]
Sensor.num = 256;%�g�����X�f���[�T��[-]
thickness_Min = 4;%[mm]
thickness_Max = 10;%[mm]
iteration = 30;
% �O���b�h�ݒ�
[ Grid ] = make_grid( Sensor );
% �Z���T�ݒ�
[ Sensor ] = make_sensor_2board( Sensor );

us = zeros(length(Grid.x),length(Grid.y),iteration,thickness_Max);
RMSE = zeros(1,iteration,thickness_Max);

for thickness = thickness_Min:2:thickness_Max%���b�w�̌���[mm]
    % �}���ݒ�
    [ Medium, Layer ] = make_layered_medium_vertical( thickness, Grid );
    % ���e�f�[�^�쐬
    [ projection ] = getProjectionData( Sensor, Grid, Medium );
    % �����č\��
    pjt_est = zeros(Sensor.num,Sensor.num);
    p = zeros(2,length(Grid.x)+length(Grid.y));%��_���W�Q
    reI = zeros(length(Grid.x),length(Grid.y));
    reI_store = zeros(length(Grid.x),length(Grid.y),iteration);
    det_I = zeros(length(Grid.x),length(Grid.y));
    for it = 1:iteration
        count = ones(length(Grid.x));
        for ii = 1:Sensor.num
            %���M�f�q�̍��W�ݒ�(x_tr,y_tr)
            x_tr = Sensor.pos(1,ii);
            y_tr = Sensor.pos(2,ii);
            pos_tr = [x_tr, y_tr];
            for jj = 1:Sensor.num
                if (ii == jj) || (1<=ii)&&(ii<=Sensor.num/2) && (1<=jj)&&(jj<=Sensor.num/2) || (Sensor.num/2<ii && Sensor.num/2<jj)
                    continue
                else
                    %��M�f�q�̍��W�ݒ�(x_re,y_re)
                    p_length = zeros(length(Grid.x),length(Grid.y));
                    x_re = Sensor.pos(1,jj);
                    y_re = Sensor.pos(2,jj);
                    pos_re = [x_re, y_re];
                    all_length = norm(pos_tr - pos_re);
                    zz = zeros(length(Grid.x),length(Grid.y));
                    %�f�q�����钼���̕�����
                    x_line = ((x_tr-x_re) / (y_tr-y_re)) * (Grid.y-y_re) + x_re;
                    y_line = ((y_tr-y_re) / (x_tr-x_re)) * (Grid.x-x_re) + y_re;
                    %�e�f�q�ʒu���O���b�h�ɓ��Ă͂߂�
                    [~,x_tr_index] = min(abs(x_tr - Grid.x));
                    [~,y_tr_index] = min(abs(y_tr - Grid.y));
                    [~,x_re_index] = min(abs(x_re - Grid.x));
                    [~,y_re_index] = min(abs(y_re - Grid.y));
                    %�ϕ��J�n�ʒu�C�I���ʒu�̌���i�����`�̈�):�o�H�����钷���`
                    x_start = min(x_tr_index,x_re_index);
                    x_end = max(x_tr_index,x_re_index);
                    y_start = min(y_tr_index,y_re_index);
                    y_end = max(y_tr_index,y_re_index);
                    %�v�Z�̈撆�̊i�q�Ƃ̌�_���i�[(p)
                    p(1,1:length(y_start:y_end)) = x_line(y_start:y_end);
                    p(2,1:length(y_start:y_end)) = Grid.y(y_start:y_end);
                    p(1,length(y_start:y_end)+1:length(y_start:y_end)+length(x_start:x_end)) = Grid.x(x_start:x_end);
                    p(2,length(y_start:y_end)+1:length(y_start:y_end)+length(x_start:x_end)) = y_line(x_start:x_end);
                    p(:,length(y_start:y_end)+length(x_start:x_end)+1:length(p)) = [];
                    rm_ind = p(1,:)>Sensor.sizeTotal/2 | p(1,:)<-Sensor.sizeTotal/2 | p(2,:)>Sensor.sizeTotal/2 | p(2,:)<-Sensor.sizeTotal/2;
                    p(:,rm_ind) = [];
                    p=rmmissing(p,2);%NaN���o�����߁C�Ώ��D
                    %�������W�ݒ�(x�̒l���ŏ��ȓ_�j
                    if max(p(1,:)) - min(p(1,:)) < Grid.size % argx(min(p))������������ł��Ȃ����Ƃ�h�����߁D
                        [~,b] = min(p(2,:));
                        p_init = p(:,b);
                        while ~isempty(p)
                            %������������Z���̌��o
                            x_ind = find(abs(p_init(1,1) - Grid.x) == 0, 1);
                            y_ind = find(abs(p_init(2,1) - Grid.y) == 0, 1);
                            if isempty(x_ind)
                                [~,x_ind] = min(abs(p_init(1,1) - Grid.x - Grid.size/2));
                            end
                            if isempty(y_ind)
                                [~,y_ind] = min(abs(p_init(2,1) - Grid.y - Grid.size/2));
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
                            x_ind = find(abs(p_init(1,1) - Grid.x) == 0, 1);
                            y_ind = find(abs(p_init(2,1) - Grid.y) == 0, 1);
                            if isempty(x_ind)
                                [~,x_ind] = min(abs(p_init(1,1) - Grid.x - Grid.size/2));
                            end
                            if isempty(y_ind)
                                [~,y_ind] = min(abs(p_init(2,1) - Grid.y - Grid.size/2));
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
        us(:,:,it,thickness) = Grid.size./(reI+Grid.size/Medium.v0);
        det_I = zeros(length(Grid.x),length(Grid.y));
        RMSE(1,it,thickness) = sqrt(immse(us(:,:,it,thickness), Medium.sound_speed));
    end
end
toc