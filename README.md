# Chương trình Tách Nguồn Âm Thanh Thiếu Xác Định (Underdetermined Source Separation)

Đây là một chương trình **MATLAB** triển khai thuật toán **tách nguồn âm thanh** (*audio source separation*) trong **điều kiện thiếu xác định**, nghĩa là **số lượng nguồn âm thanh (sources)** lớn hơn **số lượng microphone**.

Trong cấu hình mặc định, chương trình **mô phỏng 3 nguồn âm thanh và sử dụng 2 microphone (stereo mix)**.

Thuật toán hoạt động theo **phương pháp lặp lại (iterative)**, sử dụng **Phân tích Thành phần Độc lập (ICA – Independent Component Analysis)** kết hợp với **mặt nạ thời gian-tần số (Time-Frequency Masks)** và các **tiêu chí dừng dựa trên năng lượng và tính hai tai (binaural criteria)**.

---

## 1. Khởi tạo và Cấu hình (Initialization)

Phần này thiết lập các tham số cơ bản và cấu hình cho STFT.

```matlab
%% main.m
% Chương trình chính: Underdetermined source separation
% Tương thích MATLAB 2016
clear all; close all; clc;
format compact;

dis = 1;
if dis, disp('Initialisation...'); end

% Cấu hình ngẫu nhiên để lặp lại được kết quả
rng('default');
rng(1,'twister');

% Tham số chung
M = 2;            % số microphone
u = 0.5;          % cardioid control
N = 3;            % số nguồn để tạo mixture khi evalu=1
th = 1;
stopthresholdini = 3000;
TC1 = 0.1;
TC2 = 0.03;
numlags = 1;
thepow = 20;
minpow = 30;

evalu = 1; % chạy chế độ tách + đánh giá
```

**Chức năng chính:**  
- **Tham số chung:** Đặt số microphone (`M=2`), số nguồn (`N=3`) và ngưỡng dừng (`stopthresholdini`, `th`).  
- **evalu = 1:** Bật chế độ tạo hỗn hợp ngẫu nhiên và đánh giá hiệu suất.  
- **Tham số STFT:** Xác định kích thước FFT (`NFFT=2048`), cửa sổ (`WINDOW` – mặc định là Hamming) và độ chồng lấn (`NOVERLAP`) để phân tích phổ.

---

## 2. Tải / Tạo Dữ liệu Hỗn hợp (Load / Create Sources and Stereo Mix)

Phần này tạo ra hỗn hợp âm thanh (stereo mix) `X` từ các nguồn âm thanh mẫu, mô phỏng quá trình thu âm thực tế.

```matlab
%% === Load / create sources and stereo mix ===
if evalu
    % Load 12 source files (A-L)
    s = [];
    [s(:,1),fs]  = audioread('ukma.wav'); % A
    [s(:,2),~]   = audioread('frma.wav'); % B
    ... (tải thêm 10 file nữa)
```

**Chức năng chính:**  
- **Tạo hỗn hợp:** Chọn ngẫu nhiên các nguồn và hướng, sau đó sử dụng ma trận hỗn hợp `A` (tính bằng `calcA`) để tạo ra hỗn hợp stereo `X`.  
- **Ngưỡng năng lượng:** Ước tính ngưỡng năng lượng (`thE`, `minpower`) để xác định tín hiệu giọng nói đáng tin cậy.  
- **Mặt nạ lý tưởng:** Tính toán mặt nạ lý tưởng (`imaskL`, `imaskR`) để đánh giá khi `evalu=1`.

---

## 3. Vòng lặp Tách nguồn Chính (Main Separation Loop)

Đây là phần lõi của thuật toán, nơi quá trình tách nguồn lặp đi lặp lại.

- **Làm trắng PCA (whitening):** Chuẩn bị dữ liệu để ICA hoạt động ổn định hơn.  
- **ICA:** Sử dụng `icaML` để tìm các thành phần độc lập.  
- **Tạo / Áp dụng Mặt nạ:** Dùng `applymasks` để tách tín hiệu.  
- **Tiêu chí dừng:**  
  - `enerstop`: Dừng khi năng lượng quá thấp.  
  - `oneortwo_cond`: Dừng khi tín hiệu đạt điều kiện hai tai.

---

## 4. Hậu xử lý và Đánh giá (Post-processing and Evaluation)

Sau khi tách xong, chương trình tính toán hiệu suất:

- **comparemasks:** So sánh mặt nạ đầu ra với mặt nạ lý tưởng.  
- **calcELNR:** Tính các chỉ số SNR, năng lượng và cải thiện chất lượng tách.  
- **Lưu trữ:** Kết quả được lưu thành `.mat`.

---

## Các File Chức năng Bắt buộc (Required .m Files)

| File | Mô tả |
|------|-------|
| `calcA.m` | Tính ma trận hỗn hợp A dựa trên góc và mô hình microphone |
| `icaML.m` | Thực hiện ICA (Independent Component Analysis) |
| `idealmask.m` | Tính mặt nạ lý tưởng (ideal TF mask) |
| `colorimask.m` | Hiển thị mặt nạ lý tưởng có màu |
| `applymasks.m` | Áp dụng mặt nạ để tách tín hiệu |
| `oneortwo_cond.m` | Tính tiêu chí dừng hai tai |
| `enerstop.m` | Tính tiêu chí dừng theo năng lượng |
| `getfinalmask.m` | Tạo mặt nạ cuối cùng |
| `multisigcheck.m` | Hợp nhất các tín hiệu trùng |
| `nosigcorr.m` | Loại bỏ tín hiệu tương quan cao |
| `comparemask.m` | So sánh mặt nạ đầu ra và lý tưởng |
| `calcELNR.m` | Tính toán chỉ số cải thiện SNR |


# 🧠 Báo cáo cải tiến và giải thích các Figure trong chương trình tách tín hiệu (22/10/2025)

## I. Các cải tiến trong phiên bản hiện tại

Dưới đây là các thay đổi và tối ưu đã được thực hiện so với phiên bản gốc:

| Nhóm cải tiến | Mô tả chi tiết | Mục đích |
|----------------|----------------|-----------|
| **1️⃣ Chuẩn hóa dữ liệu đầu vào (`Xn`, `Xm`)** | Thêm kiểm tra kích thước ma trận trước khi trừ trung bình, xử lý tự động theo chiều hàng/cột. | Tránh lỗi “Matrix dimensions must agree”. |
| **2️⃣ Sửa lỗi trong `nosigcorr.m`** | - Kiểm tra và đồng bộ độ dài tín hiệu.<br>- Đảm bảo chỉ lấy 2 kênh trái/phải nếu có (`s(:,1)` và `s(:,2)`). | Giải quyết lỗi “Subscripted assignment dimension mismatch”. |
| **3️⃣ Kiểm tra file âm thanh trước khi đọc** | Thêm `exist(str, 'file')` trước `audioread`. Nếu không tồn tại thì bỏ qua hoặc cảnh báo. | Tránh lỗi “filename specified was not found”. |
| **4️⃣ Tự động tạo thư mục kết quả (`result_folder`)** | Tự động `mkdir(result_folder)` nếu chưa có. | Đảm bảo lưu kết quả an toàn. |
| **5️⃣ Lưu toàn bộ Figure ra file ảnh (PNG)** | Sau khi vẽ xong, tự động lưu bằng:<br>`saveas(gcf, fullfile(result_folder, sprintf('figure%d.png', figIndex)));` | Giúp xem lại kết quả mà không cần giữ session MATLAB. |
| **6️⃣ Dọn dẹp bộ nhớ** | Sử dụng `close all` sau khi lưu figure. | Giảm chiếm dụng RAM và tránh lỗi khi chạy nhiều lần. |

---

## II. Giải thích Figure 1 và Figure 2

Hai hình này là **biểu đồ phổ tần (Spectrogram)** thể hiện năng lượng tín hiệu theo **thời gian – tần số**, dùng để quan sát hiệu quả tách nguồn âm.

---

### 🎧 Figure 1 – Phổ năng lượng của tín hiệu gốc (Original Spectrogram)

#### 🔹 Nội dung hiển thị:
- **Trục hoành (X-axis)**: Thời gian (giây)  
- **Trục tung (Y-axis)**: Tần số (Hz)  
- **Màu sắc (Color)**: Cường độ năng lượng (biên độ) – càng sáng, năng lượng càng mạnh.

#### 🔹 Ý nghĩa:
- Biểu diễn tín hiệu **hỗn hợp đầu vào**, chứa cả hai (hoặc nhiều) nguồn âm trộn lẫn.
- Các **vệt năng lượng chồng lấn** cho thấy vùng giao thoa của các nguồn.
- Giúp ta hình dung **mức độ trộn** trước khi áp dụng thuật toán tách.

#### 🔹 Mục tiêu:
Dùng để **so sánh trực quan** với Figure 2 – nhằm đánh giá hiệu quả tách tín hiệu.

---

### 🎚️ Figure 2 – Phổ năng lượng sau tách (Separated Spectrogram)

#### 🔹 Nội dung hiển thị:
- Cấu trúc trục tương tự Figure 1.  
- Hiển thị dữ liệu sau khi áp dụng **thuật toán tách nguồn âm** (masking, ICA, hoặc năng lượng tương quan).

#### 🔹 Ý nghĩa:
- Thể hiện **kết quả tách**: năng lượng của từng nguồn âm được tách riêng rõ ràng hơn.
- Vùng tần số chồng lấn đã giảm, phổ năng lượng gọn và rõ ràng hơn.
- Nếu thuật toán hiệu quả, phổ của mỗi nguồn sẽ **ít nhiễu, rõ biên độ đặc trưng**.

#### 🔹 Mục tiêu:
Giúp **đánh giá trực quan hiệu quả của thuật toán** bằng cách so sánh với Figure 1.

---

## III. Cách lưu Figure thủ công (nếu muốn)

Nếu bạn muốn lưu thủ công thay vì tự động:

```matlab
figure(1);
% ... code vẽ ...
title('Spectrogram of Original Signal');
saveas(gcf, 'result/figure1_original.png');

figure(2);
% ... code vẽ ...
title('Spectrogram of Separated Signal');
saveas(gcf, 'result/figure2_separated.png');
