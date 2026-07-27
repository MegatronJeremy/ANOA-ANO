@echo off
REM setup.bat -- cmd.exe shim around setup.ps1 (for double-click / plain cmd use).
REM Forwards all arguments and preserves the exit code. Not a second code path.
cd /d "%~dp0"
powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0setup.ps1" %*
exit /b %ERRORLEVEL%
