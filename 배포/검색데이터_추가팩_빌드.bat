@echo off
chcp 65001 >nul
setlocal enabledelayedexpansion
title PRIME-DS 검색 데이터 추가팩 빌드

rem ====================================================================
rem  검색(Search) 모듈 데이터를 별도 ZIP 으로 묶습니다.
rem
rem  묶는 것 : drug_index.db  +  data\CID*.csv (16,000개 이상)
rem  받는 쪽 : PRIME-DS 포터블 폴더에 그대로 풀고 재색인 한 번
rem
rem  먼저 «배포_빌드.bat» 을 돌려 두면 그 런타임을 씁니다.
rem  아직 안 돌렸다면 이 PC 의 Python 을 씁니다.
rem
rem  실행 방법은 배포_빌드.bat 과 같습니다.
rem    배포_빌드.bat "<프로젝트 폴더 경로>"  처럼 경로를 줄 수도 있습니다.
rem ====================================================================

set "KIT=%~dp0"
if "%KIT:~-1%"=="\" set "KIT=%KIT:~0,-1%"

if not exist "%KIT%\_build_searchpack.py" (
    echo   [실패] !KIT!\_build_searchpack.py 가 없습니다.
    echo          빌드 키트 7개 파일이 같은 폴더에 있어야 합니다.
    pause
    exit /b 1
)

rem --- 프로젝트 폴더 찾기 (p00_flask_app.py 기준) ---------------------
set "PROJ="
if not "%~1"=="" (
    if exist "%~1\p00_flask_app.py" set "PROJ=%~1"
)
if not defined PROJ (
    for %%P in ("%KIT%\..") do if exist "%%~fP\p00_flask_app.py" set "PROJ=%%~fP"
)
if not defined PROJ (
    if exist "%KIT%\p00_flask_app.py" set "PROJ=%KIT%"
)
if not defined PROJ (
    echo   [실패] PRIME-DS 프로젝트 폴더를 찾지 못했습니다.
    echo          프로젝트 폴더를 인자로 주세요. 예:
    echo            검색데이터_추가팩_빌드.bat "<프로젝트 폴더 경로>"
    pause
    exit /b 1
)

set "OUT=%PROJ%\_dist"
set "PY=%OUT%\PRIME-DS\runtime\python.exe"

if not exist "%PY%" (
    where py >nul 2>nul
    if errorlevel 1 (
        set "PY=python"
    ) else (
        set "PY=py"
    )
    echo   [안내] 배포본 런타임이 없어 이 PC 의 Python 을 사용합니다.
)

echo.
echo  ==============================================================
echo   검색 데이터 추가팩 빌드
echo  ==============================================================
echo   프로젝트 : %PROJ%
echo   출력     : %OUT%
echo.
echo   파일이 많아 10 ~ 30 분쯤 걸립니다.
echo  ==============================================================
echo.

"!PY!" "%KIT%\_build_searchpack.py" "%PROJ%" "%OUT%"
if errorlevel 1 (
    echo.
    echo   [실패] 위 메시지를 확인해 주세요.
    pause
    exit /b 1
)

echo.
echo  ==============================================================
echo   완료 : %OUT%\PRIME-DS_검색데이터팩.zip
echo  ==============================================================
echo.
pause
