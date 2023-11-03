clc;
clear;

for i = 1 : 1 : 5
    % 創建一個簡單的圖像
    x = 1:10;
    y = x.^2;
    plot(x, y);
    
    % 指定要儲存的圖像檔案名稱和路徑
    file_name =  sprintf('plot_%d.png', i);  % 檔案名稱可以是您喜歡的名稱，包括路徑
    file_path = 'C:\Users\JOU\Desktop\git\Special-project\Cycloidal Gear\gif';  % 將路徑替換為您想要儲存的路徑
    
    % 創建完整的路徑
    full_file_path = fullfile(file_path, file_name);
    
    % 使用saveas函數將圖像儲存為指定的檔案類型（例如.png、.jpg等）
    saveas(gcf, full_file_path);
    
    % 清除畫面（可選）
    clf;

end
