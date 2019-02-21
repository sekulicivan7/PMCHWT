
function set_ed_2_tri=edge2trian(obj)


Ne = size(obj.edges,2);
Nt = size(obj.topol,2);

set_ed_v1 = obj.trian(1,:);
i_neg_v1 = find(set_ed_v1<0);
set_ed_v1(i_neg_v1) = abs(set_ed_v1(i_neg_v1)) + Ne;

set_ed_v2 = obj.trian(2,:);
i_neg_v2 = find(set_ed_v2<0);
set_ed_v2(i_neg_v2) = abs(set_ed_v2(i_neg_v2)) + Ne;

set_ed_v3 = obj.trian(3,:);
i_neg_v3 = find(set_ed_v3<0);
set_ed_v3(i_neg_v3) = abs(set_ed_v3(i_neg_v3)) + Ne;

%%% set_ed_2_tri
set_ed_2_tri = zeros(1,3*Nt);
set_ed_2_tri(1:3:end) = set_ed_v1;
set_ed_2_tri(2:3:end) = set_ed_v2;
set_ed_2_tri(3:3:end) = set_ed_v3; 