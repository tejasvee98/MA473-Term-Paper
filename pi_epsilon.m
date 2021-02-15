function x=pi_epsilon(y,epsilon)
    c_0 = (35/256)*epsilon; c_1 = 0.5 ; c_2 = 35/(64*epsilon); c_4 = -35/(128*epsilon^3);
    c_6= 7/(64*epsilon^5); c_8 = -5/(256*epsilon^7);
    c_3=0;c_5=0;c_7=0;c_9=0;
    p=[c_9,c_8,c_7,c_6,c_5,c_4,c_3,c_2,c_1,c_0];
    x=y;
    if (y>=epsilon)
        x=y;
    elseif (abs(y)<epsilon)
        x=polyval(p,y);
    elseif (y<=-epsilon) 
        x=zeros(1,length(y));
    end
end