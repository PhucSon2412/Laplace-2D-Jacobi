clc             % Xóa cửa sổ lệnh
clear           % Xóa tất cả các biến khỏi workspace
close all       % Đóng tất cả các figure đang mở

% Nhập chiều dài chiều rộng của miền
a = input('Nhập chiều dài:');
b = input('Nhập chiều rộng:');

% Nhập số điểm lưới và sai số
nx = input('Nhập số lượng điểm lưới theo trục x:');
ny = input('Nhập số lượng điểm lưới theo trục y:');
tolerance = input('Nhập ngưỡng sai số:');

% Đặt kích thước của miền hình chữ nhật     
a = a*pi;
b = b*pi;

% Khoảng cách giữa các điểm lưới theo chiều x và y
dx = a/nx;
dy = b/ny;

% Tạo các vector chứa tọa độ x và y
x = 0:dx:a;
y = 0:dy:b;

% Thiết lập các điều kiện biên cho ma trận U
Co = 1 / (2 * (dx^2 + dy^2));         % Hệ số Courant cho phương pháp số
U = zeros(ny+1, nx+1);                % Khởi tạo ma trận U
U(1, :) = 0;
U(ny+1, :) = sin(x) / sin(a);
U(:, 1) = 0;
U(:, nx+1) = sinh(y) / sinh(b);

% Thêm biến sai số và giá trị ngưỡng
error = Inf;                            % Khởi tạo sai số ban đầu là vô cùng
iter_number = 0;                        % Khởi tạo số lần lặp
max_iterations = 10000;                 % Số lần lặp tối đa để tránh vòng lặp vô hạn
errors = [];                            % Mảng để lưu giá trị error sau mỗi lần lặp
while error > tolerance && iter_number < max_iterations
    U_old = U;                        % Lưu lại ma trận U của lần lặp trước
    maxU_old = max(U_old, [], 'all');
    nU_old = U_old / maxU_old;
    % Tính giá trị của hàm U tại mỗi điểm lưới dựa trên phương pháp lặp số cho phương trình Laplace
    for i = 2:ny
        for j = 2:nx
            U(i, j) = Co * (dx^2 * (U(i+1, j) + U(i-1, j)) + dy^2 * (U(i, j+1) + U(i, j-1)));
        end
    end
    maxU = max(U, [], 'all');    % Tìm giá trị lớn nhất trong ma trận U
    nU = U / maxU;               % Chuẩn hóa ma trận U để giá trị lớn nhất bằng 1
    E = nU_old - nU;                 % Tính ma trận sai số
    % Tính sai số giữa hai lần lặp liên tiếp
    error = max(max(E));
    errors = [errors, error];    % Lưu giá trị error vào mảng
    iter_number = iter_number + 1;  % Tăng số lần lặp
end

% Vẽ biểu đồ Numerical
subplot(2, 1, 1);
    contourf(nU, 200, 'linecolor', 'none')
    title('Numerical')
    xlabel('Chiều dài(mm)')
    ylabel('Chiều rộng(mm)')
    colormap(jet(256))
    colorbar
    clim([-1, 1])

% Vẽ đồ thị biểu diễn giá trị error sau mỗi lần lặp
subplot(2, 1, 2);
plot(errors, 'LineWidth', 2);
xlabel('Iteration Number');
ylabel('Error');
title('Sai số qua từng lần lặp(V)');
grid on;

% Vẽ đồ thị Numerical 3D
figure();
[X, Y] = meshgrid(1:size(nU, 2), 1:size(nU, 1));
surf(X, Y, nU, 'EdgeColor', 'none');
title('Numerical');
xlabel('Chiều dài(mm)');
ylabel('Chiều rộng(mm)');
zlabel('Điện thế');
colormap(jet(256));
colorbar;
clim([-1, 1]);

% Hiển thị số lần lặp và sai số cuối cùng
disp(['Số lần lặp: ', num2str(iter_number)])
disp(['Sai số cuối cùng: ', num2str(error)])
