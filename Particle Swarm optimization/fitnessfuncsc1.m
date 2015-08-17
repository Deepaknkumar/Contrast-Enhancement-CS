 function z = fitnessfuncsc1(u,hi,lamda,gamma,dh,sm)
            dh = dh*transpose(u);
            dh = sum(dh);
            s1 = sum((transpose(u)-hi).^2);
            s2 = sum((transpose(u)-sm).^2);
            s3 = sum(dh.^2);
            z = s1 + lamda*s2 + gamma*s3;