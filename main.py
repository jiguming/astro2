import streamlit as st
import numpy as np
from astropy.io import fits
from PIL import Image
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import Time
from datetime import datetime

# --- 앱 페이지 설정 ---
st.set_page_config(page_title="천문 이미지 분석기", layout="wide", initial_sidebar_state="expanded")
st.title("🔭 천문 이미지 처리 앱")
st.caption("💡 좌측 사이드바에서 다른 분석 기능으로 이동할 수 있어요.")

# --- 파일 업로더 ---
uploaded_file = st.file_uploader("분석할 FITS 파일을 선택하세요.", type=['fits', 'fit', 'fz'])

# --- 서울 위치 (고정) ---
seoul_location = EarthLocation(lat=37.5665, lon=126.9780, height=50)
now = datetime.utcnow()
now_astropy = Time(now)

# --- 파일 처리 로직 ---
if uploaded_file:
    try:
        with fits.open(uploaded_file) as hdul:
            image_hdu = next((h for h in hdul if h.data is not None and h.is_image), None)
            if image_hdu is None:
                st.error("유효한 이미지 데이터를 찾을 수 없습니다.")
            else:
                header = image_hdu.header
                data = np.nan_to_num(image_hdu.data)
                st.success(f"**'{uploaded_file.name}'** 파일을 성공적으로 처리했습니다.")
                
                # 좌우 컬럼
                col1, col2 = st.columns(2)
                with col1:
                    st.header("이미지 정보")
                    st.text(f"크기: {data.shape[1]} x {data.shape[0]} 픽셀")
                    if 'OBJECT' in header:
                        st.text(f"관측 대상: {header['OBJECT']}")
                    if 'EXPTIME' in header:
                        st.text(f"노출 시간: {header['EXPTIME']} 초")
                    st.header("물리량")
                    mean_brightness = np.mean(data)
                    st.metric("이미지 전체 평균 밝기", f"{mean_brightness:.2f}")

                with col2:
                    st.header("이미지 미리보기")
                    if data.max() == data.min():
                        norm_data = np.zeros_like(data, dtype=np.uint8)
                    else:
                        scale_min, scale_max = np.percentile(data, 5), np.percentile(data, 99.5)
                        data_clipped = np.clip(data, scale_min, scale_max)
                        norm_data = ((data_clipped - scale_min) / (scale_max - scale_min) * 255).astype(np.uint8)
                    img = Image.fromarray(norm_data)
                    st.image(img, caption="FITS 이미지", use_container_width=True)

                # --- 사이드바: 천체 위치 계산 ---
                st.sidebar.header("🧭 현재 천체 위치 (서울 기준)")
                if 'RA' in header and 'DEC' in header:
                    try:
                        coord = SkyCoord(ra=header['RA'], dec=header['DEC'], unit=('hourangle', 'deg'))
                        altaz = coord.transform_to(AltAz(obstime=now_astropy, location=seoul_location))
                        st.sidebar.metric("고도 (°)", f"{altaz.alt.degree:.2f}")
                        st.sidebar.metric("방위각 (°)", f"{altaz.az.degree:.2f}")
                    except Exception as e:
                        st.sidebar.warning(f"천체 위치 계산 실패: {e}")
                else:
                    st.sidebar.info("RA/DEC 정보가 없습니다.")
    except Exception as e:
        st.error(f"파일 처리 중 오류: {e}")
        st.warning("FITS 형식이 유효하지 않을 수 있습니다.")
else:
    st.info("시작하려면 FITS 파일을 업로드해주세요.")

# --- 💬 댓글 기능 ---
st.divider()
st.header("💬 의견 남기기")
if "comments" not in st.session_state:
    st.session_state.comments = []

with st.form(key="comment_form"):
    name = st.text_input("이름 입력", key="name_input")
    comment = st.text_area("댓글 입력", key="comment_input")
    if st.form_submit_button("댓글 남기기"):
        if name.strip() and comment.strip():
            st.session_state.comments.append((name.strip(), comment.strip()))
            st.success("댓글이 저장되었습니다.")
        else:
            st.warning("이름과 댓글을 모두 입력해주세요.")

if st.session_state.comments:
    st.subheader("📋 전체 댓글")
    for i, (n, c) in enumerate(reversed(st.session_state.comments), 1):
        st.markdown(f"**{i}. {n}**: {c}")
else:
    st.info("아직 댓글이 없습니다. 첫 댓글을 남겨보세요!")
