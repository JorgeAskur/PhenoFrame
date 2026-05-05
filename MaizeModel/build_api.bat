@echo off
setlocal

set "ROOT=%~dp0"
pushd "%ROOT%" >nul 2>nul
if errorlevel 1 (
    echo Build failed: cannot enter %ROOT%
    exit /b 1
)

set "OUT_DLL=%CD%\maize_c_api.dll"

if exist "a.exe" del /f /q "a.exe" >nul 2>nul

echo Building %OUT_DLL%
g++ -shared -O2 -std=c++17 -DMAIZE_C_API_BUILD -I. maize_c_api.cpp Maize.cpp Descriptor.cpp vect3d.cpp tinyxml2.cpp -o "%OUT_DLL%"
if errorlevel 1 (
    echo Build failed during compilation.
    popd >nul
    exit /b 1
)

if not exist "%OUT_DLL%" (
    echo Build failed: expected output not created: %OUT_DLL%
    popd >nul
    exit /b 1
)

echo Build complete: %OUT_DLL%
popd >nul
exit /b 0
