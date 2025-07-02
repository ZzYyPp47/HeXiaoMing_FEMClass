function [A,b] = fixed_presure(A,b,func_real_p,nodes_before_fixed,fixed_index,Pb_for_p)
A(nodes_before_fixed + fixed_index,:) = 0;
A(nodes_before_fixed + fixed_index, nodes_before_fixed + fixed_index) = 1;
b(nodes_before_fixed + fixed_index) = func_real_p(Pb_for_p(1,fixed_index),Pb_for_p(2,fixed_index),0,0);
end