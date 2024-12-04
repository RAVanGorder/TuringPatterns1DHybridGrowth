function W = extract_gen_fourier_coff(u, xvec, cutoff)
% Given a solution vector u, this code computes the magnitude of each of the spatial modes

W = zeros(1,cutoff);
len_int = xvec(end);

for i=1:cutoff
    efun = sqrt(2/len_int) * cos(pi*i*xvec/len_int);
    W(i) = trapz(xvec, u.*efun);
end 
end
