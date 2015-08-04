function Atc = atc_input(t,dataset)

    switch dataset
        case '16_8'
            t0 = 0;
            p1=120; 
            p2=120; 
        case '14_7'
            t0 = 0;
            p1=120; 
            p2=120; 
        case '14_9' % this dataset is split into 3 segments           
            SecondPartTime = 10 *(60+20); %10 on/off periods
            ThirdPartTime = SecondPartTime + 20*(60+10);
            if(t < SecondPartTime)
                t0 = 0;
                p1=60; 
                p2=20; 
            elseif(t < ThirdPartTime) 
                t0 = SecondPartTime;
                p1=60; 
                p2=10; 
            else
                t0 = ThirdPartTime;
                p1 = 0;
                p2 = 0;
            end           
        case '13_9' % this dataset is split into 2 parts
            SecondPartTime = 7* (60 + 60); % 7 initial periods
            if(t < SecondPartTime)
                t0 = 0;
                p1=60;
                p2=60; 
            else 
                t0 = SecondPartTime;
                p1 = 60;
                p2 = 30;
            end
    end
   
    p=p1+p2;
    tp = mod(t-t0,p); 
    
    Atc(tp <= p1) = 2; %ng/mL
    Atc(tp > p1) = 0; %ng/mL     
end