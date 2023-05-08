REM set paths
set PATH=%PREFIX%;%PREFIX%\Library\mingw-w64\bin;%PREFIX%\Library\usr\bin;%PREFIX%\Library\bin;%PREFIX%\Scripts;%PREFIX%\bin;%PATH%

REM compile Python files
echo ""
echo "Compiling Python files"
echo "======================"
pushd ${PREFIX}\lib
${PREFIX}\bin\libtbx.py_compile_all -i
popd
echo "done"

REM write phenix_env files
echo ""
echo "Writing phenix_env files"
echo "========================"
SCRIPT_LOCATION=`${PREFIX}\bin\libtbx.find_in_repositories libtbx\auto_build\conda_build\write_env_files.py`
${PREFIX}\bin\python ${SCRIPT_LOCATION} ^
  --program ${INSTALLER_NAME} ^
  --prefix ${PREFIX} ^
  --bin-dir Library\bin ^
  --version ${INSTALLER_VER} ^
  --destination ${PREFIX}
if exist "${PREFIX}/${INSTALLER_NAME}_env.bat" echo "${PREFIX}/${INSTALLER_NAME}_env.bat created"
echo "done"

REM rebuild rotarama and cablam caches
echo ""
echo "Rebuilding rotarama and cablam caches"
echo "====================================="
${PREFIX}/bin/mmtbx.rebuild_rotarama_cache
${PREFIX}/bin/mmtbx.rebuild_cablam_cache
echo "done"

REM build documentation
echo ""
echo "Building Phenix documentation"
echo "============================="
${PREFIX}/bin/phenix_html.rebuild_docs
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
