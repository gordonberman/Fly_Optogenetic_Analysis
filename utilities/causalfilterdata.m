function out = causalfilterdata(data,alpha)

    L = length(data);	
    if iscolumn(data)
	data = [data;zeros(L,1)];
    else 
	data = [data zeros(1,L)];
    end	

    xx = ((1:2*L) - round(L));
    if iscolumn(data)
        xx = xx';
    end
    
    g = alpha^2.*xx.*exp(-alpha.*xx);
    g(g<0) = 0;
    %g = exp(-.5.*xx.^2./sigma^2)/sqrt(2*pi*sigma^2);
    out = fftshift(ifft(fft(data).*fft(g)));
    
    
    out = out(1:L);