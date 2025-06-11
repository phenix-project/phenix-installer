REM unset problematic environment variables
set "LIBTBX_BUILD="

REM set paths
set PATH=%PREFIX%;%PREFIX%\Library\mingw-w64\bin;%PREFIX%\Library\usr\bin;%PREFIX%\Library\bin;%PREFIX%\Scripts;%PREFIX%\bin;%PATH%

REM compile Python files
echo ""
echo "Compiling Python files"
echo "======================"
pushd %PREFIX%\lib
call %PREFIX%\Library\bin\libtbx.py_compile_all.bat -i
popd
echo "done"

REM create custom bin directory
echo ""
echo "Creating phenix_bin directory"
echo "============================="
set SCRIPT_LOCATION=%PREFIX%\Lib\site-packages\libtbx\auto_build\conda_build\create_custom_bin.py
%PREFIX%\python.exe %SCRIPT_LOCATION% ^
  --prefix %PREFIX% ^
  --custom_bin %INSTALLER_NAME%_bin ^
  --packages phenix dials dxtbx
if %ERRORLEVEL% neq 0 exit 1
echo "done"

REM write phenix_env files
echo ""
echo "Writing phenix_env files"
echo "========================"
set SCRIPT_LOCATION=%PREFIX%\Lib\site-packages\libtbx\auto_build\conda_build\write_env_files.py
%PREFIX%\python.exe %SCRIPT_LOCATION% ^
  --program %INSTALLER_NAME% ^
  --prefix %PREFIX% ^
  --bin-dir %INSTALLER_NAME%_bin ^
  --version %INSTALLER_VER% ^
  --destination %PREFIX%
if %ERRORLEVEL% neq 0 exit 1
if exist "%PREFIX%\\%INSTALLER_NAME%_env.bat" echo "%PREFIX%\\%INSTALLER_NAME%.bat created"
if exist "%PREFIX%\\%INSTALLER_NAME%_env.bat" copy "%PREFIX%\\%INSTALLER_NAME%.bat" "%PREFIX%\\%INSTALLER_NAME%.sh"
if exist "%PREFIX%\\%INSTALLER_NAME%_env.sh" echo "%PREFIX%\\%INSTALLER_NAME%.sh created"
echo "done"

REM rebuild rotarama and cablam caches
echo ""
echo "Rebuilding rotarama and cablam caches"
echo "====================================="
call %PREFIX%\Library\bin\mmtbx.rebuild_rotarama_cache.bat
if %ERRORLEVEL% neq 0 exit 1
call %PREFIX%\Library\bin\mmtbx.rebuild_cablam_cache.bat
if %ERRORLEVEL% neq 0 exit 1
echo "done"

REM build documentation
echo ""
echo "Building Phenix documentation"
echo "============================="
set PYTHONIOENCODING="utf8"
call %PREFIX%\Library\bin\phenix_html.rebuild_docs.bat
if %ERRORLEVEL% neq 0 exit 1
echo "done"

REM print copy of conclusion.txt without unix text
echo ""
echo "=========================================================================="
echo ""
echo "                      Phenix installation complete"
echo "                      ----------------------------"
echo ""
echo "You can begin using Phenix now by setting your environment with the"
echo "command:"
echo ""
echo "  Windows users:"
echo "    <installation directory>\phenix_env.bat"
echo ""
echo "To use Phenix, go to a work directory and type:"
echo ""
echo "  phenix"
echo ""
echo "On Windows, you can add a shortcut to the batch file to your Startup"
echo "folder."
echo ""
echo "You have successfully installed Phenix."
echo ""
echo "=========================================================================="
echo ""
