function normalizer=calc_normal(A,data,theta_c)
data_size=size(data,1);
normalizer=0;
for i=1:data_size
    normalizer=normalizer+calc_p_all_x(A,data(i,:),theta_c);
end
end