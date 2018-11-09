tic
%%%% �����p�����[�^ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CT�`��
ring_cx = 0; ring_cy = 0;% ���S�ʒu
ring_radius = 50.e-3;% ���a
rr =  ring_radius;
%�Z���T�ݒu
t_num = 16;%�g�����X�f���[�T��
t = linspace(0, ((t_num-1)/t_num)*2*pi, t_num);%�Z���T�ʒu�p
t_pos = zeros(2, t_num);%�Z���T�ʒu
t_pos(1,:) = rr*cos(t)+ring_cx;
t_pos(2,:) = rr*sin(t)+ring_cy;
%�Z���ݒ�
cell_num = 16;
cell_size = ring_radius*2 / (cell_num-1);
x_grid = -ring_radius: cell_size : ring_radius;
y_grid = x_grid;

p = zeros(2,length(x_grid)+length(y_grid));
C = zeros(t_num^2, cell_num^2);

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
                C(t_num*(ii-1)+jj,cell_num*(x_ind-1)+y_ind) = C(t_num*(ii-1)+jj,cell_num*(x_ind-1)+y_ind) + 1;
                p_init = p_neighbor;
            end
        end
    end
end
A = sum(C,1);
figure;
plot(A);
axis tight
ylim([0 45]);
xlabel('��f')
ylabel('�o�H���x')
exportfig( 'H:\result\matrix_rank_visualization\2018_03_02_matrix_rank_visualization_ring_16_2','png',[300,150]);