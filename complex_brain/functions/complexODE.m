function dzc = complexODE(t,z,coe,SC) 
    z1 = zeros(379,1);
    for k1 = 1:379
        z1(k1) = z(k1) + 1i*z(379+k1);
    end
    dz = construct_data2(z1,t,coe,SC); 
    dzc = [real(dz);imag(dz)];
end
