function W = extract_gen_fourier_coff(u, xvec, cutoff)
% Given a solution vector u with, computes the magnitude of each of the 
% 1-cutoff generalised fourier modes

W = zeros(1,cutoff);
len_int = xvec(end);

for i=1:cutoff
    efun = sqrt(2/len_int) * cos(pi*i*xvec/len_int);
    W(i) = trapz(xvec, u.*efun);
end 
end