# Hướng Dẫn Tách Tiếng 2 Người từ 2 Microphone

## Vấn đề cần giải quyết

**Đầu vào**: 2 file microphone (mic1.wav và mic2.wav) thu cùng lúc 2 người nói  
**Đầu ra**: 2 file riêng biệt (nguoi_1.wav và nguoi_2.wav) chứa tiếng của từng người

## Cách các thuật toán nhận biết 2 người khác nhau

### 1. BSS-ICA (Phương pháp chính - Khuyến nghị)

**Nguyên lý**: Dựa trên tính độc lập thống kê

**Cách nhận biết**:
- Mỗi người có **pattern nói** khác nhau
- **Phổ tần số** khác nhau (formant khác nhau)
- **Thời điểm nói** khác nhau (khi người này nói, người kia có thể im)
- ICA tìm các thành phần độc lập trong tín hiệu hỗn hợp

**Ví dụ**:
```
Mic 1: Người A + Người B (hỗn hợp 1)
Mic 2: Người A + Người B (hỗn hợp 2, với tỷ lệ khác)

ICA tìm cách biến đổi:
Output 1: Chủ yếu là Người A
Output 2: Chủ yếu là Người B
```

**Điều kiện tốt**:
- 2 người nói **không đồng thời** (ít overlap)
- 2 người có **giọng khác nhau** (cao/thấp, nam/nữ)
- **Vị trí** 2 microphone khác nhau (để có 2 hỗn hợp khác nhau)

### 2. GSC/Delay-and-Sum Beamformer

**Nguyên lý**: Dựa trên hướng không gian

**Cách nhận biết**:
- Mỗi người ở một **hướng/vị trí** khác nhau
- Sóng âm đến 2 microphone với **thời gian trễ** khác nhau
- Công thức: `delay = (khoảng_cách_mic * sin(góc)) / 343m/s`

**Ví dụ**:
```
Người A ở bên trái → sóng đến mic 1 trước, mic 2 sau
Người B ở bên phải → sóng đến mic 2 trước, mic 1 sau

Beamforming hướng trái → tăng cường Người A
Beamforming hướng phải → tăng cường Người B
```

**Điều kiện tốt**:
- 2 người ở **các hướng khác nhau** (góc > 30°)
- Khoảng cách giữa 2 mic phù hợp (4-5cm)
- Biết hoặc ước tính được hướng của mỗi người

### 3. Clustering (IPD/ILD)

**Nguyên lý**: Dựa trên đặc trưng binaural

**Cách nhận biết**:
- **IPD** (Độ lệch pha): Người ở bên trái có IPD khác người ở bên phải
- **ILD** (Độ lệch biên độ): Người gần mic 1 sẽ to hơn ở mic 1
- Clustering nhóm các điểm có IPD/ILD tương tự

**Ví dụ**:
```
Người A ở trái → IPD < 0, ILD > 0
Người B ở phải → IPD > 0, ILD < 0

Clustering tách thành 2 nhóm → 2 người
```

**Điều kiện tốt**:
- 2 người ở **các hướng khác nhau**
- Khoảng cách mic phù hợp

---

## Script nên dùng

### **Script đơn giản nhất: `tach_2_nguoi.m`** (Khuyến nghị)

```matlab
tach_2_nguoi
```

**Script này kết hợp 3 phương pháp để tách tốt nhất**:

#### Bước 1-2: Tải và tiền xử lý
- Tải mic1.wav và mic2.wav
- Pre-emphasis, bandpass (300-3400 Hz)

#### Bước 3: Tách bằng nhiều phương pháp
1. **BSS-ICA**: Tách dựa trên tính độc lập thống kê
2. **Clustering (IPD/ILD)**: Tách dựa trên đặc trưng binaural
3. **Beamforming**: Tách dựa trên hướng không gian

#### Bước 4: Chọn kết quả tốt nhất
- Thử tất cả các nguồn từ 3 phương pháp
- Tính correlation giữa từng cặp
- Chọn 2 nguồn có correlation thấp nhất (khác nhau nhất)

#### Bước 5: Loại bỏ phần dính (Wiener Masking)
- Nếu correlation > 0.3: áp dụng Wiener masking
- Tạo mask dựa trên công suất: `mask1 = P1 / (P1 + P2)`
- Loại bỏ phần chung giữa 2 nguồn

#### Bước 6: Voice Activity Detection (VAD)
- Phát hiện phần có giọng nói
- Loại bỏ phần im lặng và nhiễu

#### Bước 7: Orthogonalization
- Làm 2 nguồn vuông góc với nhau
- Loại bỏ hoàn toàn phần dính còn lại

#### Bước 8-10: Lưu và đánh giá
- Lưu nguoi_1.wav và nguoi_2.wav
- Hiển thị kết quả chi tiết
- Đánh giá chất lượng tách

**Ưu điểm**:
- Kết hợp nhiều phương pháp → tăng khả năng tách được
- Tự động chọn kết quả tốt nhất
- Loại bỏ phần dính bằng nhiều phương pháp
- Có VAD để loại bỏ im lặng
- Đánh giá chi tiết chất lượng

---

## Tại sao có thể chưa tách được?

### 1. **2 người nói đồng thời quá nhiều**
- ICA cần các nguồn độc lập
- Nếu 2 người nói cùng lúc, khó tách

**Giải pháp**: 
- Sử dụng Clustering hoặc Beamforming
- Cải thiện tiền xử lý

### 2. **2 microphone thu tín hiệu quá giống nhau**
- Nếu 2 mic ở quá gần hoặc cùng hướng
- ICA cần 2 hỗn hợp khác nhau

**Giải pháp**: 
- Đặt microphone ở vị trí khác nhau hơn
- Sử dụng phương pháp khác ngoài ICA

### 3. **Chất lượng tín hiệu kém**
- Nhiễu quá lớn
- Volume quá nhỏ

**Giải pháp**:
- Cải thiện tiền xử lý
- Điều chỉnh tham số lọc

### 4. **Tham số không phù hợp**
- Bandpass filter quá hẹp/rộng
- Noise gate quá cao/thấp

**Giải pháp**:
- Điều chỉnh tham số trong script

---

## So sánh kết quả

Sau khi chạy, kiểm tra:

1. **Correlation giữa 2 output**:
   - < 0.3: Rất tốt
   - 0.3 - 0.5: Khá tốt
   - > 0.5: Chưa tốt, còn dính

2. **Xem spectrogram**:
   - 2 spectrogram có pattern khác nhau → tốt
   - 2 spectrogram giống nhau → chưa tách được

3. **Nghe thử**:
   - Nghe rõ từng người riêng → tốt
   - Nghe cả 2 người trong cùng file → chưa tách được

---

## Nếu kết quả chưa tốt

### Thử script khác:

1. **Script kết hợp nhiều phương pháp**:
```matlab
speech_separation_combined
```
- Chạy 7 phương pháp
- Tự động chọn kết quả tốt nhất

2. **Điều chỉnh tham số trong `tach_2_nguoi.m`**:
- Dòng 59: Thay đổi bandpass filter (300-3400 Hz → 200-4000 Hz)
- Dòng 114: Thay đổi VAD threshold
- Dòng 152: Thay đổi noise gate

3. **Kiểm tra file đầu vào**:
- Đảm bảo mic1.wav và mic2.wav thực sự là 2 microphone khác nhau
- Kiểm tra có thực sự có 2 người nói không
- Đảm bảo chất lượng thu tốt

