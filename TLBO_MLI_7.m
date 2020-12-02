%clc;
%clear all;
%close all;


ma=0.01:0.01:1;

for aa=1:length(ma)

    XX=linspace(0,pi/2,4);
    xmin=XX(1,1:3);
    xmax=XX(1,2:4);

    n=length(xmax);
    N=10*n;
    maxit=5000;
    cost_t=0;
    for i=1:N

        Ipop(i,:)=(rand(1,n).*(xmax-xmin)+xmin);
        x=Ipop(i,:);

        if x(1)<0.0175
           cost_t=1e10;
        end
        if x(1)>=x(2)
            cost_t=1e10;
        end
        if x(2)>=x(3)
            cost_t=1e10;
        end

        if x(3)>((pi/2)-0.0175)
            cost_t=1e10;
        end
        if x(2)-x(1)<0.0175
            cost_t=1e10;
        end
        if x(3)-x(2)<0.0175
            cost_t=1e10;
        end

        cost(i,1)=cost_t+(100*((3*ma(aa))-(cos(x(1))+cos(x(2))+cos(x(3))))/(3*ma(aa))).^4 +...
                     (1/5)*((50*(1/5)*(cos(5*x(1))+cos(5*x(2))+cos(5*x(3)))/(cos(x(1))+cos(x(2))+cos(x(3)))).^2)+...
                     (1/7)*((50*(1/7)*(cos(7*x(1))+cos(7*x(2))+cos(7*x(3)))/(cos(x(1))+cos(x(2))+cos(x(3)))).^2); % Fitness Function

        cost_t=0;
    end

     Ipop_mix=[Ipop,cost];
     Ipop_sort=sortrows(Ipop_mix,n+1);
     best=Ipop_sort(1,:);
     worst=Ipop_sort(N,:);
     xmin=[0 0 0];
     xmax=[pi/2 pi/2 pi/2]; 

    for iter=1:maxit

        for i=1:N

            MM=mean(Ipop);
            F=round(1+rand());    
            Dx=rand(1,n).*(best(1,1:n)-F*MM);

            Ipop_TA=Ipop(i,:)+Dx;
            Ipop_TA=min([xmax;Ipop_TA]);
            Ipop_TA=max([xmin;Ipop_TA]);

            x=Ipop_TA;

            if x(1)<0.0175
               cost_t=1e10;
            end
            if x(1)>=x(2)
                cost_t=1e10;
            end
            if x(2)>=x(3)
                cost_t=1e10;
            end

            if x(3)>((pi/2)-0.0175)
                cost_t=1e10;
            end
            if x(2)-x(1)<0.0175
                cost_t=1e10;
            end
            if x(3)-x(2)<0.0175
                cost_t=1e10;
            end

            cost_TA=cost_t+(100*((3*ma(aa))-(cos(x(1))+cos(x(2))+cos(x(3))))/(3*ma(aa))).^4 +...
                     (1/5)*((50*(1/5)*(cos(5*x(1))+cos(5*x(2))+cos(5*x(3)))/(cos(x(1))+cos(x(2))+cos(x(3)))).^2)+...
                     (1/7)*((50*(1/7)*(cos(7*x(1))+cos(7*x(2))+cos(7*x(3)))/(cos(x(1))+cos(x(2))+cos(x(3)))).^2); % Fitness Function

            cost_t=0;

            if cost_TA<cost(i,1)
               Ipop(i,:)=Ipop_TA;
               cost(i,1)= cost_TA;
            end 


            k=randperm(N);
            k1=find(k~=i);
            kk=k(k1);

            if cost(i,1)<cost(kk(1),1)
                Dx_new=rand(1,n).*(Ipop(i,:)-Ipop(kk(1),:));
                Ipop_st=Ipop(i,:)+Dx_new;
            else
                Dx_new=rand(1,n).*(Ipop(kk(1),:)-Ipop(i,:));
                Ipop_st=Ipop(i,:)+Dx_new;
            end

            Ipop_st=min([xmax;Ipop_st]);
            Ipop_st=max([xmin;Ipop_st]);

            x=Ipop_st;

            if x(1)<0.0175
               cost_t=1e10;
            end
            if x(1)>=x(2)
                cost_t=1e10;
            end
            if x(2)>=x(3)
                cost_t=1e10;
            end

            if x(3)>((pi/2)-0.0175)
                cost_t=1e10;
            end
            if x(2)-x(1)<0.0175
                cost_t=1e10;
            end
            if x(3)-x(2)<0.0175
                cost_t=1e10;
            end

            cost_st=cost_t+(100*((3*ma(aa))-(cos(x(1))+cos(x(2))+cos(x(3))))/(3*ma(aa))).^4 +...
                     (1/5)*((50*(1/5)*(cos(5*x(1))+cos(5*x(2))+cos(5*x(3)))/(cos(x(1))+cos(x(2))+cos(x(3)))).^2)+...
                     (1/7)*((50*(1/7)*(cos(7*x(1))+cos(7*x(2))+cos(7*x(3)))/(cos(x(1))+cos(x(2))+cos(x(3)))).^2); % Fitness Function

              cost_t=0;

            if cost_st<cost(i,1)
               Ipop(i,:)=Ipop_st;
               cost(i,1)=cost_st;
            end 

        end

        Ipop_mix=[Ipop,cost];
        Ipop_sort=sortrows(Ipop_mix,n+1);
        best=Ipop_sort(1,:);
        worst=Ipop_sort(N,:);
        final(iter,:)=best;

    end

    %% Results
    V1star=3*ma(aa);
    V1= cos(best(1))+cos(best(2))+cos(best(3));
    V1nesbi=V1/V1star;
    V5=abs((1/5)*(cos(5*best(1))+cos(5*best(2))+cos(5*best(3)))/V1);
    V7=abs((1/7)*(cos(7*best(1))+cos(7*best(2))+cos(7*best(3)))/V1);
    V11=abs((1/11)*(cos(11*best(1))+cos(11*best(2))+cos(11*best(3)))/V1);
    V13=abs((1/13)*(cos(13*best(1))+cos(13*best(2))+cos(13*best(3)))/V1);
    V17=abs((1/17)*(cos(17*best(1))+cos(17*best(2))+cos(17*best(3)))/V1);
    V19=abs((1/19)*(cos(19*best(1))+cos(19*best(2))+cos(19*best(3)))/V1);
    V23=abs((1/23)*(cos(23*best(1))+cos(23*best(2))+cos(23*best(3)))/V1);
    V25=abs((1/25)*(cos(25*best(1))+cos(25*best(2))+cos(25*best(3)))/V1);
    V29=abs((1/29)*(cos(29*best(1))+cos(29*best(2))+cos(29*best(3)))/V1);
    V31=abs((1/31)*(cos(31*best(1))+cos(31*best(2))+cos(31*best(3)))/V1);
    V35=abs((1/35)*(cos(35*best(1))+cos(35*best(2))+cos(35*best(3)))/V1);
    V37=abs((1/37)*(cos(37*best(1))+cos(37*best(2))+cos(37*best(3)))/V1);
    V41=abs((1/41)*(cos(41*best(1))+cos(41*best(2))+cos(41*best(3)))/V1);
    V43=abs((1/43)*(cos(43*best(1))+cos(43*best(2))+cos(43*best(3)))/V1);
    V47=abs((1/47)*(cos(47*best(1))+cos(47*best(2))+cos(47*best(3)))/V1);
    V49=abs((1/49)*(cos(49*best(1))+cos(49*best(2))+cos(49*best(3)))/V1);

    THD=100.*sqrt((V5^2)+(V7^2)+(V11^2)+(V13^2)+(V17^2)+(V19^2)+(V23^2)+(V25^2)+(V29^2)+(V31^2)+(V35^2)+(V37^2)+(V41^2)+(V43^2)+(V47^2)+(V49^2));


    best(1,1:3)=best(1,1:3)/pi*180;
    Switching_Angles=best(1,1:3);
    fitness_best=best(1,end);
    Modulation_Index=ma(aa);


    Objective_Function_Value(aa)=fitness_best;
    
    Fundamental_Harmonic(aa)=V1nesbi;
    
    V5th(aa)=V5;
    V7th(aa)=V7;
    V11th(aa)=V11;
    V13th(aa)=V13;
    V17th(aa)=V17;
    V19th(aa)=V19;
    V23th(aa)=V23;
    V25th(aa)=V25;
    
    teta1(aa)=best(1);
    teta2(aa)=best(2);
    teta3(aa)=best(3);
    
    
    Total_Harmonic_Distortion(aa)=THD;

end



    figure(6);
    semilogy(ma,Objective_Function_Value,'LineWidth',3)
    xlabel('Modulation Index')
    ylabel('ObjectiveFunction Value')
    grid
    
    figure(7);
    plot(ma,Fundamental_Harmonic.*100,'LineWidth',3)
    ylim([0 120])
    xlabel('Modulation Index')
    ylabel('Fundamental Harmonic')
    grid
    
    figure(8);
    plot(ma,V5th,'b','LineWidth',2)
    hold on
    plot(ma,V7th,'r','LineWidth',2)
    hold on
    legend('5th','7th')
    xlabel('Modulation Index')
    ylabel('Harmonics')
    grid
    
    
    figure(9);
    plot(ma,teta1,'b','LineWidth',2)
    hold on
    plot(ma,teta2,'r','LineWidth',2)
    hold on
    plot(ma,teta3,'m','LineWidth',2)
    legend('\theta1','\theta2','\theta3')
    xlabel('Modulation Index')
    ylabel('Optimum Switching Angles (Degree)')
    grid  
    
    figure(10);
    plot(ma,Total_Harmonic_Distortion,'LineWidth',2)
    xlabel('Modulation Index')
    ylabel('Total Harmonic Distortion (%)')
    grid