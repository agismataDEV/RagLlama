for i = 15:15
A1 = importdata(['..\David_21deg_x_Lagrangian_Only_3\output_file_pressure',int2string_ICS(i),'.txt']);
A2 = importdata(['..\David_21deg_x_Lagrangian_Only_2\output_file_pressure',int2string_ICS(i),'.txt']);
sum(sum(A1-A2>1E-10))
end

