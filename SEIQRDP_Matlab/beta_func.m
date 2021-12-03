function Beta = beta_func(t,beta1,beta2,beta3,wave)
wave_1_end = days(datetime([2020,10,8])-datetime([2020,4,18]));
wave_2_end = days(datetime([2021,3,30])-datetime([2020,4,18]));
% t_true = t;
if wave == 1
    t_true = t; % Get true time from 0 date
elseif wave == 2
    t_true = t + wave_1_end; % Offset true time from zero date since t will start at 0
elseif wave ==3
    t_true = t + wave_2_end;
end
% data starts at 18 April 2020
% level5 = 0 + days(datetime([2020,5,1]) - datetime([2020,4,18]));
level5 = 13;
% level4 = level5 + days(datetime([2020,6,1]) - datetime([2020,5,1]));
level4 = 44;
% level3 = level4 + days(datetime([2020,8,18]) - datetime([2020,6,1]));
level3 = 122;
% level2 = level3 + days(datetime([2020,9,21]) - datetime([2020,8,18]));
level2 = 156;
% level1 = level2 + days(datetime([2020,12,28]) - datetime([2020,9,21]));
level1 = 254;
% level3_a = level1 + days(datetime([2021,2,28]) - datetime([2020,12,28]));
level3_a = 316;
% level1_a = level3_a + days(datetime([2021,5,30]) - datetime([2021,2,28]));
level1_a = 591;
% level2_a = level1_a + days(datetime([2021,6,15]) - datetime([2021,5,30]));
level2_a = 607;
% level3_a2 = level2_a + days(datetime([2021,6,27]) - datetime([2021,6,15]));
level3_a2 = 619;
% level4_a = level3_a2 + days(datetime([2021,7,25]) - datetime([2021,6,27]));
level4_a = 647;
% level3_a3 = level4_a + days(datetime([2021,9,12]) - datetime([2021,7,25]));
level3_a3 = 696;
% level2_a2 = level3_a3 + days(datetime([2021,9,30]) - datetime([2021,9,12]));
level2_a2 = 714;
% level1_a2 = level2_a2 + days(datetime("today") - datetime([2021,9,30]));
level1_a2 = 775;
if t_true< level5
%     disp("In level 5")
    Beta_m = 0.468;
elseif (t_true>=level5)&(t_true<level4)
%     disp("In level 4")
    Beta_m = 0.681;
elseif(t_true>=level4)&(t_true<level3)
%     disp("In level 3")
%     if wave==2
%         Beta_m = 0.95;
%     else
      Beta_m = 1.0;
%     end
elseif (t_true>=level3)&(t_true<level2)
%     disp("In level 2")
    Beta_m = 1.319;
elseif (t_true>=level2)&(t_true<level1)
%     disp("In level 1")
    Beta_m = 1.638;
elseif (t_true>=level1)&(t_true<level3_a)
    Beta_m = 0.681;
elseif (t_true>=level3_a)&(t_true<level1_a)
    Beta_m = 1.319;
elseif (t_true>=level1_a)&(t_true<level2_a)
    Beta_m = 1.0;
elseif (t_true>=level2_a)&(t_true<level3_a2)
    Beta_m = 1.0;
elseif (t_true>=level3_a2)&(t_true<level4_a)
    Beta_m = 0.468;
elseif (t_true>=level4_a)&(t_true<level3_a3)
    Beta_m = 1.681;
elseif (t_true>=level3_a3)&(t_true<level2_a2)
    Beta_m = 1.0;
elseif (t_true>=level2_a2)&(t_true<level1_a2)
    Beta_m = 1.319;
else
    Beta_m = 1.0;
end
% Beta_m
Beta = Beta_m*(beta1 + beta2*exp(-beta3*t));
end