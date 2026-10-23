@echo off
setlocal EnableExtensions
cd /d "%~dp0..\.."

set "OUTDIR=data\external_rice\Robertson2026"
set "OUT=%OUTDIR%\GSE313726_OsLeafStressIntegrated.rds"
set "URL=https://ftp.ncbi.nlm.nih.gov/geo/series/GSE313nnn/GSE313726/suppl/GSE313726_OsLeafStressIntegrated.rds"
set "EXPECTED=897951595"

if not exist "%OUTDIR%" mkdir "%OUTDIR%"

echo ============================================================
echo Robertson 2026 rice stress processed Seurat object
echo Supports resume with curl -C -
echo Expected bytes: %EXPECTED%
echo Output: %OUT%
echo ============================================================
echo.

curl.exe -L --fail --retry 5 --retry-delay 5 -C - "%URL%" -o "%OUT%"
if errorlevel 1 (
  echo.
  echo DOWNLOAD_FAILED_OR_INTERRUPTED
  echo Re-run this same BAT to resume.
  exit /b 1
)

for %%A in ("%OUT%") do set "ACTUAL=%%~zA"
echo.
echo Actual bytes: %ACTUAL%
if not "%ACTUAL%"=="%EXPECTED%" (
  echo SIZE_MISMATCH
  echo Expected %EXPECTED% bytes but got %ACTUAL%.
  echo Re-run this same BAT to resume or repair the download.
  exit /b 2
)

echo DOWNLOAD_COMPLETE_SIZE_OK
endlocal
