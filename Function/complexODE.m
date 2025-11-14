function dzc = complexODE(z,coe,SC) 
    z1 = zeros(379,1);
    for k1 = 1:379
        z1(k1) = z(k1) + 1i*z(379+k1);
    end
    dz = construct_data(z1,coe,SC); 
    dzc = [real(dz);imag(dz)];
end
