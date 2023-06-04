function [null_d, null_nd, null_d_pilots, null_nd_pilots] = bf_null_v(bf_d, bf_r, bf_d_pilots, bf_r_pilots)

% For data sc
null_d = zeros(48, 4);
null_nd = zeros(48, 4);
for i=1:48
 null_d_space = null(repmat((bf_d(i,:)), [4,1]));
 temp_w_1 = conj(bf_r(i,:))*null_d_space(:,1);
 temp_w_2 = conj(bf_r(i,:))*null_d_space(:,2);
 temp_w_3 = conj(bf_r(i,:))*null_d_space(:,3);
 temp_w = temp_w_1*null_d_space(:,1)+temp_w_2*null_d_space(:,2)+temp_w_3*null_d_space(:,3);
 null_d(i,:) = temp_w/vecnorm(temp_w, 2,1);


 null_nd_space = null(repmat((bf_r(i,:)), [4, 1]));  % conjugate nullspace   
 temp_w_1 = conj(bf_d(i,:))*null_nd_space(:,1);
 temp_w_2 = conj(bf_d(i,:))*null_nd_space(:,2);
 temp_w_3 = conj(bf_d(i,:))*null_nd_space(:,3);
 temp_w = temp_w_1*null_nd_space(:,1)+temp_w_2*null_nd_space(:,2)+temp_w_3*null_nd_space(:,3);
 null_nd(i,:) = temp_w/vecnorm(temp_w, 2,1);
end

% For pilots sc
null_d_pilots = zeros(4, 4);
null_nd_pilots = zeros(4, 4);
for i=1:4
 null_d_space_pilots = null(repmat((bf_d_pilots(i,:)), [4,1]));
 temp_w_1 = conj(bf_r_pilots(i,:))*null_d_space_pilots(:,1);
 temp_w_2 = conj(bf_r_pilots(i,:))*null_d_space_pilots(:,2);
 temp_w_3 = conj(bf_r_pilots(i,:))*null_d_space_pilots(:,3);
 temp_w = temp_w_1*null_d_space_pilots(:,1)+temp_w_2*null_d_space_pilots(:,2)+temp_w_3*null_d_space_pilots(:,3);
 null_d_pilots(i,:) = temp_w/vecnorm(temp_w, 2,1);


 null_nd_space_pilots = null(repmat((bf_r_pilots(i,:)), [4, 1]));  % conjugate nullspace   
 temp_w_1 = conj(bf_d_pilots(i,:))*null_nd_space_pilots(:,1);
 temp_w_2 = conj(bf_d_pilots(i,:))*null_nd_space_pilots(:,2);
 temp_w_3 = conj(bf_d_pilots(i,:))*null_nd_space_pilots(:,3);
 temp_w = temp_w_1*null_nd_space_pilots(:,1)+temp_w_2*null_nd_space_pilots(:,2)+temp_w_3*null_nd_space_pilots(:,3);
 null_nd_pilots(i,:) = temp_w/vecnorm(temp_w, 2,1);
end
end
