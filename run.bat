@echo off
REM run.bat -- cmd.exe shim around run.ps1 (for double-click / plain cmd use).
REM Forwards all arguments and preserves the exit code. Not a second code path.
REM Examples:  run           (menu)      run check      run qc -Smoke -Debug
REM            run all        run data    run test
cd /d "%~dp0"
powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0run.ps1" %*
exit /b %ERRORLEVEL%
