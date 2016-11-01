function s = getcontourlines(c)

    sz = size(c,2);    
    ii = 1;             
    jj = 1;             

    while ii < sz      
        n = c(2,ii);    
        s(jj).v = c(1,ii);        
        s(jj).x = c(1,ii+1:ii+n);
        s(jj).y = c(2,ii+1:ii+n);
        ii = ii + n + 1;         
        jj = jj + 1;             
    end

end