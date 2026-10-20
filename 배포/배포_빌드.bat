@echo off
chcp 65001 >nul
setlocal enabledelayedexpansion
title PRIME-DS 배포본 빌드

rem ====================================================================
rem  PRIME-DS 포터블 배포본 빌드
rem
rem  이 PC 에 Python 이 설치되어 있지 않아도 됩니다.
rem  python.org 에서 "embeddable" 런타임을 받아 배포본 안에 넣고,
rem  거기에 필요한 패키지를 설치하는 방식입니다.
rem
rem  필요한 것 : 인터넷 (이 빌드 한 번만) · 여유 공간 약 3 GB
rem  걸리는 시간 : 10 ~ 20 분
rem  결과물 : <프로젝트>\_dist\PRIME-DS\  와  PRIME-DS_포터블.zip
rem
rem  실행 방법
rem    (1) 이 파일이 프로젝트 안 배포\ 폴더에 있으면 그냥 더블클릭
rem    (2) 프로젝트가 다른 곳에 있으면 폴더 경로를 끌어다 놓거나
rem          배포_빌드.bat "<프로젝트 폴더 경로>"
rem ====================================================================

set "KIT=%~dp0"
if "%KIT:~-1%"=="\" set "KIT=%KIT:~0,-1%"

echo.
echo  ==============================================================
echo   PRIME-DS 포터블 배포본 빌드
echo  ==============================================================
echo.

rem --- 0) 키트 파일이 전부 있는지 먼저 본다 ---------------------------
rem     예전 판은 이 확인을 패키지 설치 단계에서야 했습니다. 런타임을
rem     10 분 내려받은 뒤에 «파일 없음» 으로 멈추므로 확인을 앞으로 옮겼습니다.
set "LACK="
for %%F in (_build.py _verify.py requirements-dist.txt) do (
    if not exist "%KIT%\%%F" set "LACK=!LACK! %%F"
)
if defined LACK (
    echo   [실패] 이 폴더에 빌드 키트 파일이 빠져 있습니다.
    echo.
    echo          폴더   : !KIT!
    echo          없는 것:!LACK!
    echo.
    echo          배포_빌드.bat 하나만 있으면 동작하지 않습니다.
    echo          아래 7개가 **같은 폴더에** 있어야 합니다.
    echo.
    echo            배포_빌드.bat
    echo            검색데이터_추가팩_빌드.bat
    echo            _build.py
    echo            _verify.py
    echo            _build_searchpack.py
    echo            requirements-dist.txt
    echo            README_배포.md
    echo.
    echo          프로젝트의 «배포» 폴더에 전부 들어 있습니다. 그 폴더의
    echo          배포_빌드.bat 을 실행하시면 됩니다.
    echo.
    pause
    exit /b 1
)

rem --- 0b) 프로젝트 폴더 찾기 -----------------------------------------
rem     기준은 p00_flask_app.py 가 있는 폴더입니다.
rem       ① 명령줄 인자  ② 키트의 부모 폴더  ③ 키트 폴더 자신
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
    echo.
    echo          p00_flask_app.py 가 있는 폴더를 찾습니다. 다음을 봤습니다:
    if not "%~1"=="" echo            · %~1
    echo            · !KIT!\..
    echo            · !KIT!
    echo.
    echo          프로젝트 폴더를 인자로 주시면 됩니다. 예:
    echo            배포_빌드.bat "<프로젝트 폴더 경로>"
    echo          (탐색기에서 폴더를 이 bat 파일 위로 끌어다 놓아도 됩니다)
    echo.
    pause
    exit /b 1
)

set "OUT=%PROJ%\_dist"
set "PKG=%OUT%\PRIME-DS"
set "RT=%PKG%\runtime"
set "WORK=%OUT%\_work"

echo   키트     : %KIT%
echo   프로젝트 : %PROJ%
echo   출력     : %OUT%
echo  ==============================================================
echo.

rem --- 0c) 다운로드 도구 확인 -----------------------------------------
where curl.exe >nul 2>nul
if errorlevel 1 (
    set "DL=powershell"
) else (
    set "DL=curl"
)
echo  [0/6] 키트 확인 완료 · 다운로드 도구 : !DL!

rem --- 1) 지난 빌드 정리 ---------------------------------------------
echo  [1/6] 지난 빌드 정리
if exist "%PKG%" rmdir /s /q "%PKG%"
mkdir "%PKG%" 2>nul
mkdir "%RT%" 2>nul
mkdir "%WORK%" 2>nul

rem --- 2) Python embeddable 런타임 내려받기 ---------------------------
rem     패치 버전은 수시로 바뀌므로 최신 것부터 차례로 시도합니다.
rem     이미 받아 둔 파일이 있으면 그것을 씁니다(프록시 우회용).
echo  [2/6] Python 3.12 embeddable 런타임 받는 중
set "PYZIP=%WORK%\python-embed.zip"
set "GOTPY="
if exist "%PYZIP%" (
    for %%S in ("%PYZIP%") do if %%~zS GTR 5000000 (
        set "GOTPY=미리 받아 둔 파일"
        echo        - !PYZIP! 를 사용합니다 ^(%%~zS bytes^)
    )
)
for %%V in (3.12.14 3.12.13 3.12.12 3.12.11 3.12.10 3.12.9 3.12.8 3.12.7) do (
    if not defined GOTPY (
        if exist "%PYZIP%" del /q "%PYZIP%"
        set "URL=https://www.python.org/ftp/python/%%V/python-%%V-embed-amd64.zip"
        if "!DL!"=="curl" (
            curl.exe -sS -f -L -o "%PYZIP%" "!URL!" >nul 2>nul
        ) else (
            powershell -NoProfile -Command "try{Invoke-WebRequest -Uri '!URL!' -OutFile '%PYZIP%' -UseBasicParsing}catch{exit 1}" >nul 2>nul
        )
        rem  오류 응답이 빈 파일로 남는 경우가 있어 크기로 확인합니다 (정상은 10 MB 쯤)
        if exist "%PYZIP%" (
            for %%S in ("%PYZIP%") do if %%~zS GTR 5000000 (
                set "GOTPY=%%V"
                echo        - Python %%V  ^(%%~zS bytes^)
            )
        )
    )
)
if not defined GOTPY (
    echo.
    echo   [실패] python.org 에서 런타임을 받지 못했습니다.
    echo          인터넷 연결 또는 사내 프록시를 확인해 주세요.
    echo          수동으로 받으시려면:
    echo            https://www.python.org/downloads/windows/
    echo            "Windows embeddable package (64-bit)" 를 받아
    echo            !PYZIP!
    echo            로 저장한 뒤 이 창을 다시 실행하세요.
    echo.
    pause
    exit /b 1
)

powershell -NoProfile -Command "Expand-Archive -Force -LiteralPath '%PYZIP%' -DestinationPath '%RT%'"
if errorlevel 1 (
    echo   [실패] 런타임 압축을 풀지 못했습니다.
    pause
    exit /b 1
)

rem --- 3) ._pth 패치 : site-packages 와 app 폴더를 경로에 추가 ---------
echo  [3/6] 런타임 경로 설정
rem  표준 라이브러리 zip 의 이름(python312.zip 등)을 실제 파일명으로 읽어 씁니다.
set "STDZIP="
for %%Z in ("%RT%\python*.zip") do set "STDZIP=%%~nxZ"
if not defined STDZIP (
    echo   [실패] 런타임 안에 표준 라이브러리 zip 이 없습니다. 압축 해제를 확인하세요.
    pause
    exit /b 1
)
for %%F in ("%RT%\python*._pth") do (
    > "%%F" echo !STDZIP!
    >>"%%F" echo .
    >>"%%F" echo Lib\site-packages
    >>"%%F" echo ..\app
    >>"%%F" echo import site
)

rem --- 4) pip 설치 -----------------------------------------------------
echo  [4/6] pip 설치
set "GETPIP=%WORK%\get-pip.py"
if not exist "%GETPIP%" (
    if "!DL!"=="curl" (
        curl.exe -sS -f -L -o "%GETPIP%" "https://bootstrap.pypa.io/pip/get-pip.py"
    ) else (
        powershell -NoProfile -Command "Invoke-WebRequest -Uri 'https://bootstrap.pypa.io/pip/get-pip.py' -OutFile '%GETPIP%' -UseBasicParsing"
    )
)
if not exist "%GETPIP%" (
    echo   [실패] get-pip.py 를 받지 못했습니다.
    echo          https://bootstrap.pypa.io/pip/get-pip.py 를 브라우저로 받아
    echo          !GETPIP! 로 저장한 뒤 다시 실행하세요.
    pause
    exit /b 1
)
"%RT%\python.exe" "%GETPIP%" --no-warn-script-location -q
if errorlevel 1 (
    echo   [실패] pip 설치에 실패했습니다.
    pause
    exit /b 1
)

rem --- 5) 패키지 설치 ---------------------------------------------------
echo  [5/6] 패키지 설치 (가장 오래 걸립니다 · RDKit 만 150 MB 쯤)
"%RT%\python.exe" -m pip install --no-warn-script-location --only-binary=:all: -r "%KIT%\requirements-dist.txt"
if errorlevel 1 (
    echo.
    echo   [실패] 패키지 설치에 실패했습니다. 위 메시지를 확인해 주세요.
    echo          사내망이라면 프록시 설정이 필요할 수 있습니다:
    echo            set HTTPS_PROXY=http://프록시주소:포트
    echo.
    pause
    exit /b 1
)

rem --- 6) 앱 파일 복사 · 검증 · 압축 ------------------------------------
echo  [6/6] 앱 파일 복사 · 검증 · 압축
"%RT%\python.exe" "%KIT%\_build.py" "%PROJ%" "%PKG%" "%KIT%"
if errorlevel 1 (
    echo.
    echo   [실패] 빌드 단계에서 멈췄습니다. 위 메시지를 확인해 주세요.
    pause
    exit /b 1
)

rmdir /s /q "%WORK%" 2>nul

echo.
echo  ==============================================================
echo   완료
echo  ==============================================================
echo   폴더 : %PKG%
echo   압축 : %OUT%\PRIME-DS_포터블.zip
echo.
echo   받는 쪽에서는 압축을 풀고
echo   «PRIME-DS 실행.bat» 을 더블클릭하면 됩니다.
echo   Python 설치는 필요 없습니다.
echo  ==============================================================
echo.
pause
