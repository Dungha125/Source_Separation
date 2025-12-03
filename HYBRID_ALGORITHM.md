# Thuật Toán Lai: Hybrid ICA + Spatial Masking

## Tổng Quan

**Đây là thuật toán KẾT HỢP thực sự, không phải chạy riêng rồi chọn!**

Thuật toán lai phát triển từ BSS-ICA, kết hợp với spatial information từ beamforming/clustering để cải thiện chất lượng tách nguồn.

---

## Ý Tưởng Chính

### Vấn đề của ICA thuần túy:
- ICA tách dựa trên độc lập thống kê
- Nhưng nếu 2 người nói đồng thời nhiều → khó tách
- Kết quả có thể bị "dính" (crosstalk)

### Giải pháp - Thuật toán lai:
**Kết hợp ICA với spatial information** để tách tốt hơn:

```
ICA (độc lập thống kê) + Spatial info (vị trí không gian) = Tách tốt hơn
```

---

## Các Bước Chi Tiết

### Bước 1: ICA Khởi Tạo
```matlab
[S_ica, A, W] = fastica_robust(X_bp);
```
- Tách nguồn ban đầu bằng ICA
- Ước tính **mixing matrix A**: `X = A * S`
- Ma trận A chứa thông tin về vị trí tương đối của các nguồn

**Kết quả**: 2 nguồn khởi tạo (có thể còn dính)

---

### Bước 2: Trích Xuất Spatial Information

#### 2.1. Từ Mixing Matrix A:
```matlab
A = [a11, a12]
    [a21, a22]
```
- Nếu `a11 > a21`: Nguồn 1 gần mic 1 hơn → bên trái
- Nếu `a22 > a12`: Nguồn 2 gần mic 2 hơn → bên phải

#### 2.2. Từ IPD/ILD (Clustering):
```matlab
IPD = angle(S2_tf / S1_tf)  % Độ lệch pha
ILD = 20*log10(|S2_tf| / |S1_tf|)  % Độ lệch biên độ
```
- Clustering IPD/ILD thành 2 nhóm
- Mỗi nhóm = 1 vị trí không gian = 1 người

**Kết quả**: 2 spatial masks (mask_spatial_1, mask_spatial_2)

---

### Bước 3: Tính Wiener Mask từ ICA

```matlab
P1 = |S1_ica|^2  # Công suất nguồn 1
P2 = |S2_ica|^2  # Công suất nguồn 2

mask_wiener_1 = P1 / (P1 + P2)
mask_wiener_2 = P2 / (P1 + P2)
```

- Chia năng lượng theo tỷ lệ
- Mask cao (gần 1) = thuộc về nguồn đó
- Mask thấp (gần 0) = thuộc về nguồn khác

**Kết quả**: 2 Wiener masks

---

### Bước 4: Kết Hợp Masks (KÉT HỢP CHÍNH)

```matlab
mask_combined_1 = mask_wiener_1 * mask_spatial_1
mask_combined_2 = mask_wiener_2 * mask_spatial_2
```

**Ý nghĩa**:
- Điểm có mask cao ở **CẢ HAI** → chắc chắn thuộc nguồn đó
- Điểm chỉ cao ở 1 mask → không chắc chắn
- Nhân 2 mask → chỉ giữ phần chắc chắn

**Soft masking** (làm mềm):
```matlab
mask = mask^alpha  # alpha > 1: cứng hơn, alpha < 1: mềm hơn
```
- Alpha = 1.5: Mask vừa phải
- Giữ được chi tiết tín hiệu

**Normalize**:
```matlab
mask_sum = mask1 + mask2
mask1 = mask1 / mask_sum
mask2 = mask2 / mask_sum
```
- Đảm bảo tổng = 1
- Không mất năng lượng

---

### Bước 5: Áp Dụng Combined Mask

```matlab
S1_separated = S_input * mask_combined_1
S2_separated = S_input * mask_combined_2
```

- Áp dụng mask lên **spectrogram gốc** (không phải ICA source)
- Tách từ signal gốc → giữ chất lượng tốt hơn

**Kết quả**: 2 nguồn đã tách tốt hơn

---

### Bước 6: Hậu Xử Lý

#### 6.1. De-emphasis:
- Đảo ngược pre-emphasis
- Khôi phục phổ tần số tự nhiên

#### 6.2. Voice Activity Detection (VAD):
- Phát hiện phần có giọng nói
- Loại bỏ im lặng

#### 6.3. Orthogonalization:
```matlab
nguoi_2 = nguoi_2 - projection(nguoi_2 onto nguoi_1)
```
- Làm 2 tín hiệu vuông góc
- Loại bỏ hoàn toàn phần dính còn lại

---

## So Sánh với Các Phương Pháp Khác

### ICA thuần túy:
```
Input → ICA → Output
```
- Ưu: Đơn giản
- Nhược: Có thể còn dính

### Beamforming thuần túy:
```
Input → Beamforming → Output
```
- Ưu: Tách theo hướng
- Nhược: Cần biết hướng chính xác

### **Thuật toán lai (Hybrid)**:
```
Input → ICA (khởi tạo) → Spatial mask (từ IPD/ILD)
      ↓                      ↓
   Wiener mask          Combined mask → Apply → Output
```
- Ưu: **Kết hợp điểm mạnh của cả 2**
- Kết quả: Tách tốt hơn cả ICA và Beamforming riêng lẻ

---

## Tại Sao Thuật Toán Lai Tốt Hơn?

### 1. **ICA cung cấp**:
- Phân tích thống kê ban đầu
- Wiener mask dựa trên công suất

### 2. **Spatial information cung cấp**:
- Thông tin vị trí không gian
- Spatial mask dựa trên IPD/ILD

### 3. **Kết hợp**:
- **ICA tốt** ở vùng 2 người nói **không đồng thời**
- **Spatial tốt** ở vùng 2 người nói **đồng thời** (phân biệt theo vị trí)
- Combined mask: `mask1 * mask2` → Chỉ giữ phần **CẢ HAI đều chắc chắn**

### 4. **Kết quả**:
- Giảm crosstalk (phần dính)
- Giữ chất lượng tín hiệu tốt hơn
- Tách được cả khi 2 người nói đồng thời

---

## Cách Sử Dụng

```matlab
hybrid_ica_beamforming
```

Script sẽ hiển thị:
- ICA khởi tạo
- Combined mask
- Kết quả cuối cùng
- Đánh giá chi tiết

---

## Đánh Giá Kết Quả

### Xem Figure:
- **Row 2**: ICA khởi tạo (có thể còn dính)
- **Row 3-4**: Kết quả cuối (sau khi kết hợp)
- **Column 3**: Mask và spectrogram

### Xem Correlation:
- < 0.15: Xuất sắc
- 0.15-0.30: Tốt
- 0.30-0.50: Khá
- > 0.50: Cần cải thiện

### Xem Spectrogram:
- 2 spectrogram **khác nhau rõ ràng** → Tốt
- Vùng màu vàng/cam ở **các vị trí khác nhau** → Tốt
- Vùng màu tối (xanh dương) nhiều → Đã loại bỏ nhiễu tốt

---

## Tùy Chỉnh

Có thể điều chỉnh các tham số trong script:

1. **Dòng 156**: `mask_alpha = 1.5`
   - Tăng (2.0): Mask cứng hơn (binary), loại bỏ dính tốt nhưng mất chi tiết
   - Giảm (1.0): Mask mềm hơn, giữ chi tiết nhưng có thể còn dính

2. **Dòng 82-83**: Bandpass filter (300-3400 Hz)
   - Có thể mở rộng (200-4000 Hz) nếu cần

3. **Dòng 224**: VAD threshold
   - Tăng: Loại bỏ nhiễu nhiều hơn
   - Giảm: Giữ lại nhiều tín hiệu hơn

