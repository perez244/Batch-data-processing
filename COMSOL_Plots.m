% %plot figure S11
% figure('Renderer', 'painters', 'Position', [300 300 1200 500])
% plot(F,S11_50,F,S11_70,F,S11_100,F,S11_250,F,S11_365,F,S11_500,F,S11_750,F,S11_1000,'--r','LineWidth',1.5)
% set(gca, 'FontName', 'Consolas')
% xlabel('Frequency(GHz)','FontSize',11,'FontWeight','bold')
% ylabel('Reflection(dB)','FontSize',11,'FontWeight','bold')
% title('S11 for varying Cu_2O thickness ')
% grid on
% lgd = legend('50nm','70nm','100nm','250nm','365nm','500nm','750nm','1000nm');
% lgd.FontSize = 14;
% 
% % % plot S21
% figure('Renderer', 'painters', 'Position', [200 300 1200 500])
% plot(F,S21_50,F,S21_70,F,S21_100,F,S21_250,F,S21_365,F,S21_500,F,S21_750,F,S21_1000,'--r','LineWidth',1.5)
% set(gca, 'FontName', 'Consolas')
% xlabel('Frequency(GHz)','FontSize',11,'FontWeight','bold')
% ylabel('Transmission(dB)','FontSize',11,'FontWeight','bold')
% title('S21 for varying Cu_2O thickness ')
% grid on
% lgd2 = legend('50nm','70nm','100nm','250nm','365nm','500nm','750nm','1000nm');
% lgd2.FontSize = 14;

% figure('Renderer', 'painters', 'Position', [300 300 1200 500])
% plot(F,S11_500,F,newS11_500,'--r','LineWidth',1.5)
% set(gca, 'FontName', 'Consolas')
% xlabel('Frequency(GHz)','FontSize',11,'FontWeight','bold')
% ylabel('Reflection(dB)','FontSize',11,'FontWeight','bold')
% title('S11 for varying Cu_2O thickness ')
% grid on
% lgd = legend('500nm','new 500nm');
% lgd.FontSize = 14;

% % plot S21
figure('Renderer', 'painters', 'Position', [200 300 1200 500])
plot(F,S21_250_A,F,S21_250,'--o','LineWidth',1.5)
set(gca, 'FontName', 'Consolas')
xlabel('Frequency(GHz)','FontSize',11,'FontWeight','bold')
ylabel('Transmission(dB)','FontSize',11,'FontWeight','bold')
title('S21 for Cu_2O and CuO ')
grid on
lgd2 =legend('\epsilon = 1','\epsilon = 18.1');
lgd2.FontSize = 14;



