@echo off
setlocal

:: Navigate to the parent directory of this batch file
cd /d "%~dp0.."

python -m src.gui planetary
pause