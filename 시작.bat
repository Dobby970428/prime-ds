@echo off
chcp 65001 >nul
title PRIME-DS

echo.
echo  ==============================================================
echo   PRIME-DS  ·  불순물 예측 및 ICH M7 유전독성 평가 플랫폼
echo  ==============================================================
echo.
echo   브라우저가 자동으로 열립니다.
echo   창을 닫거나 Ctrl+C 를 누르면 종료됩니다.
echo.

cd /d "%~dp0"

where py >nul 2>nul
if %errorlevel%==0 (
    py -3.11 p00_flask_app.py 2>nul || py p00_flask_app.py
) else (
    python p00_flask_app.py
)

if errorlevel 1 (
    echo.
    echo  ------------------------------------------------------------
    echo   실행에 실패했습니다.
    echo.
    echo   필요한 패키지가 설치되어 있는지 확인하세요:
    echo     pip install flask rdkit scikit-learn pandas openpyxl
    echo     pip install python-docx fpdf2 rdchiral requests
    echo  ------------------------------------------------------------
    echo.
    pause
)
