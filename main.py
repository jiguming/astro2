import streamlit as st
import numpy as np
import matplotlib.pyplot as plt

# Streamlit 앱 제목
st.title("외계 행성 탐사: 항성 광도 변화 시뮬레이션")

# 설명
st.write("""
이 앱은 외계 행성이 항성을 통과할 때 발생하는 광도 변화를 시뮬레이션합니다.
항성과 행성의 반지름을 조정하여 광도 변화 곡선을 확인하세요.
""")

# 입력 슬라이더
st.header("입력 매개변수")
star_radius = st.slider("항성 반지름 (태양 반지름 단위, R☉)", 
                        min_value=0.1, max_value=2.0, value=1.0, step=0.1)
planet_radius = st.slider("행성 반지름 (목성 반지름 단위, R_J)", 
                          min_value=0.1, max_value=2.0, value=1.0, step=0.1)

# 반지름 단위 변환 (태양 반지름과 목성 반지름 비율)
# 1 목성 반지름 ≈ 0.10045 태양 반지름
planet_radius_solar = planet_radius * 0.10045

# 광도 변화 계산 함수
def transit_light_curve(star_radius, planet_radius, time):
    """
    행성 통과에 따른 상대 광도 계산
    time: 정규화된 시간 배열 (-1 to 1)
    star_radius: 항성 반지름 (태양 반지름 단위)
    planet_radius: 행성 반지름 (태양 반지름 단위)
    """
    flux = np.ones_like(time)  # 기본 광도 = 1
    # 행성이 항성을 가리는 구간
    transit_mask = np.abs(time) <= 1
    # 광도 감소: (행성 반지름 / 항성 반지름)^2
    flux[transit_mask] = 1 - (planet_radius / star_radius) ** 2
    return flux

# 시간 배열 생성 (정규화된 시간, -1.5 ~ 1.5)
time = np.linspace(-1.5, 1.5, 1000)

# 광도 변화 계산
flux = transit_light_curve(star_radius, planet_radius_solar, time)

# 최대 광도 감소 비율 계산
max_flux_drop = (planet_radius_solar / star_radius) ** 2 * 100  # 퍼센트 단위

# 그래프 생성
fig, ax = plt.subplots()
ax.plot(time, flux, color='blue', label='상대 광도')
ax.set_xlabel('정규화된 시간')
ax.set_ylabel('상대 광도 (F/F₀)')
ax.set_title('행성 통과에 따른 항성 광도 변화')
ax.grid(True)
ax.legend()

# 그래프 표시
st.pyplot(fig)

# 결과 출력
st.header("결과")
st.write(f"**항성 반지름**: {star_radius:.2f} R☉")
st.write(f"**행성 반지름**: {planet_radius:.2f} R_J ({planet_radius_solar:.3f} R☉)")
st.write(f"**최대 광도 감소**: {max_flux_drop:.3f}%")

# 추가 정보
st.write("""
### 참고
- 광도 감소는 행성이 항성을 완전히 가리는 구간에서 (R_p / R_s)² 비율로 계산됩니다.
- 시간은 정규화된 단위로, 실제 통과 시간은 궤도 주기와 항성 크기에 따라 달라집니다.
- 이 모델은 단순화를 위해 림 다크닝(limb darkening) 효과를 무시합니다.
""")
