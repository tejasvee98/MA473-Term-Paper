function x=pi_epsilon(y,epsilon)
    c_0 = (35/256)*epsilon; c_1 = 0.5 ; c_2 = 35/(64*epsilon); c_4 = -35/(128*epsilon^3);
    c_6= 7/(64*epsilon^5); c_8 = -5/(256*epsilon^7);
    c_3=0;c_5=0;c_7=0;c_9=0;
    x=0;
    if (y>=epsilon)
        x=y;
    elseif (abs(y)<epsilon)
        x=c_0 + c_1*y + c_2*y^2 + c_3*y^3 + c_4*y^4 + c_5*y^5 + c_6*y^6 + c_7*y^7 + c_8*y^8 + c_9*y^9; 
    elseif (y<=-epsilon) 
        x=0;
    end
end