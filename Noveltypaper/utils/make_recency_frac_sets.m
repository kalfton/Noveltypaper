function frac_sets = make_recency_frac_sets(Fractal_info, frac_sets)
% fractal set for recency index calculation, the input can be 1 or 2, 
%if it is 2 the second one should be a struct, and the program will write the new fractal sets into the struct

shortIFI_ind_2nd = find(ismember(Fractal_info(2,:),[4,5]) & Fractal_info(6,:)<1.5 & Fractal_info(3,:)==2);
shortIFI_ind_3rd = find(ismember(Fractal_info(2,:),[4,5]) & Fractal_info(6,:)<1.5 & Fractal_info(3,:)==3);
num_short_2nd_pos = numel(shortIFI_ind_2nd);
num_short_3rd_pos = numel(shortIFI_ind_3rd);

longIFI_ind_2nd = find(ismember(Fractal_info(2,:),[4,5]) & Fractal_info(6,:)>1.5 & Fractal_info(6,:)<inf & Fractal_info(3,:)==2);
longIFI_ind_3rd = find(ismember(Fractal_info(2,:),[4,5]) & Fractal_info(6,:)>1.5 & Fractal_info(6,:)<inf & Fractal_info(3,:)==3);
num_long_2nd_pos = numel(longIFI_ind_2nd);
num_long_3rd_pos = numel(longIFI_ind_3rd);

num_2nd_pos = min([num_short_2nd_pos, num_long_2nd_pos, num_short_3rd_pos, num_long_3rd_pos]);
num_3rd_pos = num_2nd_pos;

% Sort the indices by the IFI, aka, Fractal_info(6,:)

[~,I]=sort(Fractal_info(6, shortIFI_ind_2nd));
shortIFI_ind_2nd = shortIFI_ind_2nd(I);
[~,I]=sort(Fractal_info(6, shortIFI_ind_3rd));
shortIFI_ind_3rd = shortIFI_ind_3rd(I);
[~,I]=sort(Fractal_info(6, longIFI_ind_2nd));
longIFI_ind_2nd = longIFI_ind_2nd(I);
[~,I]=sort(Fractal_info(6, longIFI_ind_3rd));
longIFI_ind_3rd = longIFI_ind_3rd(I);

shortIFI_ind_2nd = shortIFI_ind_2nd(1:num_2nd_pos);
shortIFI_ind_3rd = shortIFI_ind_3rd(1:num_3rd_pos);
longIFI_ind_2nd = longIFI_ind_2nd(end-num_2nd_pos+1:end);
longIFI_ind_3rd = longIFI_ind_3rd(end-num_3rd_pos+1:end);

frac_sets.shortIFI_ind = [shortIFI_ind_2nd, shortIFI_ind_3rd];
frac_sets.longIFI_ind = [longIFI_ind_2nd, longIFI_ind_3rd];


end