function r_RightBoundary = r_RightBoundary(r, N_phi, N_theta, r_mer_vec, r_equ)

r_grid = reshape(r, [N_phi, N_theta]);

r_RightBoundary = ones(N_phi+2, 1);
r_RightBoundary(1, 1) = r_equ(N_theta+2);
r_RightBoundary(end, 1) = r_mer_vec(N_theta+2);

for k=1:N_phi
  
  	if k == 1
  		r_RightBoundary(k+1, 1) = r_equ(N_theta+2-1) + r_grid(k+1, end) - r_grid(k, end-1);
  	elseif k == N_phi
  		r_RightBoundary(k+1, 1) = r_grid(k-1, end) + r_mer_vec(N_theta+2-1) - r_grid(k, end-1);
  	else
		r_RightBoundary(k+1, 1) = r_grid(k-1, end) + r_grid(k+1, end) - r_grid(k, end-1);
  	end

end


end