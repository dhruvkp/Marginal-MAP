function [theta_i,theta_c] = simple_ising_overcomplete_theta(A,J,X)
% input: adjacency matrix, edge weight
% output: singleton parameters, clique parameters
n_vertices=size(A,1);
theta_i=zeros(n_vertices,2);
theta_c=zeros(n_vertices,n_vertices,2,2);
for i=1:n_vertices
    for j=1:n_vertices
        if A(i,j)==0
            continue
        end 
        for xi=1:2
            for xj=1:2
                if xi==xj
                   theta_c(i,j,xi,xj)=exp(J);
                end
            end
        end
    end
end