@echo off
REM Installe l indexation automatique dans GSAS-II (Windows)
REM Installs the auto-indexing plugin into GSAS-II (Windows)
set "DEST=%USERPROFILE%\.GSASII\exports"
if not exist "%DEST%" mkdir "%DEST%"
copy /Y "%~dp0G2export_AutoIndex.py" "%DEST%\" >nul
copy /Y "%~dp0g2_autoindex.py" "%DEST%\" >nul
echo.
echo Fichiers copies dans / Files copied to: %DEST%
echo Conograph : le plugin cherche "Conograph CUI.exe" dans Program Files (sinon indiquez son chemin).
echo Conograph: the plugin looks for "Conograph CUI.exe" in Program Files (otherwise give its path).
echo Redemarrez GSAS-II / Restart GSAS-II:
echo   Export - Entire project as - (user) Auto-index (find unit cell)
echo.
pause
