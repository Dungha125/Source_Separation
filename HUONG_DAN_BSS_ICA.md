# Hướng Dẫn Sử Dụng BSS-ICA để Tách Nguồn Âm Thanh

## Mô tả
**Mục tiêu: Tách được nguồn nói của mỗi người từ tín hiệu hỗn hợp**

Chương trình sử dụng nhiều phương pháp để tách các nguồn âm thanh (đặc biệt tối ưu cho giọng nói) từ tín hiệu hỗn hợp thu được từ 2 microphone (mic1.wav và mic2.wav).

**Xem file `NGUYEN_LY_TACH_NGUON.md` để hiểu chi tiết cách các thuật toán nhận biết và tách nguồn âm từ 2 người khác nhau.**

Các phương pháp được sử dụng:

- **BSS-ICA** (Blind Source Separation - Independent Component Analysis) - Dựa trên tính độc lập thống kê
- **GSC Beamformer** (Generalized Sidelobe Canceller) - Dựa trên hướng không gian
- **Clustering** (dựa trên IPD/ILD) - Dựa trên đặc trưng binaural
- **Delay-and-Sum Beamformer** - Dựa trên hướng không gian đơn giản
- **Differential Microphone Array** - Dựa trên sự khác biệt giữa các microphone
- **MVDR Beamformer** - Dựa trên covariance matrix
- **LCMV Beamformer** - Dựa trên constraints tuyến tính

**Đặc điểm tối ưu cho giọng nói:**
- Pre-emphasis filter để tăng cường tần số cao của giọng nói
- Bandpass filter (300-3400 Hz) cho dải tần giọng nói
- Tự động chọn hướng tốt nhất cho GSC
- Hậu xử lý để loại bỏ nhiễu và cải thiện chất lượng

## Cách sử dụng

### Phương pháp 1: Thuật toán lai nâng cao (Khuyến nghị nhất)
Chạy file `advanced_hybrid_separation.m` trong MATLAB:

```matlab
advanced_hybrid_separation
```

**Thuật toán lai KẾT HỢP TẤT CẢ các phương pháp beamforming!**

**Các bước chính**:
1. **Ước tính hướng từ IPD/ILD**: Clustering để tìm hướng của 2 người
2. **Áp dụng TẤT CẢ beamformers**: 
   - Delay-and-Sum
   - GSC
   - MVDR
   - LCMV
   - Differential
   - ICA
3. **Chọn cặp tốt nhất**: Correlation thấp VÀ khác input
4. **Ideal Binary Mask**: Tạo mask từ cặp đã chọn
5. **Iterative refinement**: Lặp 5 lần với Wiener + Spatial mask
6. **Hậu xử lý mạnh**: Orthogonalization lặp 5 lần + VAD

**Ưu điểm**:
- Kết hợp TẤT CẢ các phương pháp → chọn tốt nhất
- Binary masking mạnh (mask^2.5) → loại bỏ dính tốt
- Iterative refinement → cải thiện dần
- Đảm bảo output KHÁC input
- Hiển thị chi tiết từng bước

**Output**: `nguoi_1.wav` và `nguoi_2.wav` trong thư mục `output_advanced_hybrid/`

### Phương pháp 2: Script kết hợp đơn giản
Chạy file `tach_2_nguoi.m` trong MATLAB:

```matlab
tach_2_nguoi
```

Script này:
- Chạy 3 phương pháp riêng biệt (ICA, Clustering, Beamforming)
- Chọn 2 nguồn có correlation thấp nhất
- Loại bỏ phần dính bằng Wiener Masking
- **Output**: `nguoi_1.wav` và `nguoi_2.wav`

### Phương pháp 3: Script kết hợp nhiều thuật toán (Nâng cao)
Chạy file `speech_separation_combined.m` trong MATLAB:

```matlab
speech_separation_combined
```

**Đây là script chính để tách tiếng 2 người riêng biệt!**

Script này:
- **Kết hợp tất cả các thuật toán**: BSS-ICA, GSC Beamformer, và Clustering
- **Tự động chọn kết quả tốt nhất**: Tính correlation giữa các kết quả và chọn 2 kết quả khác nhau nhất
- **Output cuối cùng**: 2 file riêng biệt `nguoi_1.wav` và `nguoi_2.wav`
- Tiền xử lý tối ưu (pre-emphasis, bandpass)
- Hậu xử lý để cải thiện chất lượng
- Hiển thị kết quả chi tiết

### Phương pháp 4: Script tối ưu cho giọng nói
Chạy file `speech_separation_optimized.m` trong MATLAB:

```matlab
speech_separation_optimized
```

Script này được tối ưu hóa đặc biệt để tách giọng nói, bao gồm:
- Tiền xử lý tối ưu (pre-emphasis, bandpass)
- Tự động chọn hướng tốt nhất
- Hậu xử lý để cải thiện chất lượng
- Hiển thị kết quả chi tiết

### Phương pháp 5: Script đơn giản
Chạy file `bss_ica_separate.m` trong MATLAB:

```matlab
bss_ica_separate
```

Script này sẽ:
1. Tự động tải `mic1.wav` và `mic2.wav`
2. Tiền xử lý tín hiệu (lọc high-pass, chuẩn hóa)
3. Áp dụng BSS-ICA để tách nguồn
4. Áp dụng GSC Beamformer để tách nguồn
5. Lưu kết quả vào thư mục `result_bss_ica/`

### Phương pháp 6: Script đầy đủ
Chạy file `main.m` trong MATLAB:

```matlab
main
```

Script này bao gồm:
- Tách nguồn bằng Clustering
- Tách nguồn bằng BSS-ICA
- Tách nguồn bằng GSC Beamformer
- Kết quả được lưu trong thư mục `result_enhanced/`

## Yêu cầu

1. **File đầu vào:**
   - `mic1.wav` - Tín hiệu từ microphone 1
   - `mic2.wav` - Tín hiệu từ microphone 2
   - Cả hai file phải có trong thư mục gốc

2. **File hỗ trợ cần thiết:**
   - `fastica_robust.m` - Hàm thực hiện FastICA
   - `gsc_beamformer.m` - Hàm thực hiện GSC Beamformer
   - `my_istft.m` - Hàm Inverse STFT (chỉ cần cho main.m)

## Kết quả

Sau khi chạy, bạn sẽ nhận được:

### Với `advanced_hybrid_separation.m` (Khuyến nghị):
- **`nguoi_1.wav`** - Tiếng người thứ nhất (kết quả tốt nhất)
- **`nguoi_2.wav`** - Tiếng người thứ hai (kết quả tốt nhất)
- `ket_qua_advanced.fig` - Figure chi tiết

Script này kết hợp tất cả các beamformer và tự động chọn kết quả tốt nhất!

### Với `tach_2_nguoi.m`:
- **`nguoi_1.wav`** - Tiếng người thứ nhất
- **`nguoi_2.wav`** - Tiếng người thứ hai

Script này kết hợp ICA + Clustering + Beamforming, chọn cặp tốt nhất.

### Với `speech_separation_combined.m`:
- **`nguoi_1.wav`** - Tiếng người thứ nhất (kết quả cuối cùng)
- **`nguoi_2.wav`** - Tiếng người thứ hai (kết quả cuối cùng)

Đây là 2 file output chính, được chọn tự động từ kết quả của tất cả các phương pháp!

### Với `speech_separation_optimized.m`:
- `speaker_1_ica.wav` - Người 1 (BSS-ICA)
- `speaker_2_ica.wav` - Người 2 (BSS-ICA)
- `speaker_1_ica_enhanced.wav` - Người 1 (BSS-ICA + hậu xử lý)
- `speaker_2_ica_enhanced.wav` - Người 2 (BSS-ICA + hậu xử lý)
- `speaker_X_gsc_Xdeg.wav` - Người X (GSC Beamformer, tự động chọn hướng tốt nhất)

### Với `bss_ica_separate.m`:
- `result_bss_ica/separated_source_1.wav` - Nguồn âm thanh thứ nhất (BSS-ICA)
- `result_bss_ica/separated_source_2.wav` - Nguồn âm thanh thứ hai (BSS-ICA)
- `result_bss_ica/separated_gsc_source_1.wav` - Nguồn 1 (GSC Beamformer)
- `result_bss_ica/separated_gsc_source_2.wav` - Nguồn 2 (GSC Beamformer)
- `result_bss_ica/input_mic1.wav` - Tín hiệu mic1 đã xử lý
- `result_bss_ica/input_mic2.wav` - Tín hiệu mic2 đã xử lý

### Với `main.m`:
- `result_enhanced/output_ica_src1.wav` - Nguồn 1 (BSS-ICA)
- `result_enhanced/output_ica_src2.wav` - Nguồn 2 (BSS-ICA)
- `result_enhanced/output_cluster_src1.wav` - Nguồn 1 (Clustering)
- `result_enhanced/output_cluster_src2.wav` - Nguồn 2 (Clustering)
- `result_enhanced/output_gsc_src1.wav` - Nguồn 1 (GSC Beamformer)
- `result_enhanced/output_gsc_src2.wav` - Nguồn 2 (GSC Beamformer)

## Lưu ý

1. **Chất lượng tách nguồn** phụ thuộc vào:
   - Vị trí tương đối của các nguồn âm so với microphone
   - Mức độ độc lập của các nguồn âm
   - Mức độ nhiễu trong tín hiệu

2. **BSS-ICA hoạt động tốt nhất** khi:
   - Số lượng nguồn ≤ số lượng microphone
   - Các nguồn âm độc lập thống kê
   - Tín hiệu không có nhiễu quá lớn

3. **GSC Beamformer hoạt động tốt nhất** khi:
   - Biết được hướng của nguồn âm cần tách
   - Nguồn âm ở các hướng khác nhau
   - Khoảng cách giữa các microphone phù hợp (thường 4-5cm)

4. Nếu gặp lỗi, kiểm tra:
   - File `mic1.wav` và `mic2.wav` có tồn tại không
   - File hỗ trợ (`fastica_robust.m`, `gsc_beamformer.m`, `my_istft.m`) có trong thư mục không
   - Tín hiệu có đủ dài và có nội dung không

## Ví dụ sử dụng

```matlab
% Trong MATLAB Command Window:
cd 'E:\paper\thesis-bss-master\thesis-bss-master'
bss_ica_separate
```

Sau khi chạy xong, mở các file WAV trong thư mục `result_bss_ica/` để nghe kết quả.

