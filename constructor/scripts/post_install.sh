# compile Python files
echo ""
echo "Compiling Python files"
echo "======================"
pushd ${PREFIX}/lib
${PREFIX}/bin/libtbx.py_compile_all -i
popd

# write phenix_env files
echo ""
echo "Writing phenix_env files"
echo "========================"
SCRIPT_LOCATION=`${PREFIX}/bin/libtbx.find_in_repositories libtbx/auto_build/conda_build/write_env_files.py`

${PREFIX}/bin/python ${SCRIPT_LOCATION} \
  --program ${INSTALLER_NAME} \
  --prefix ${PREFIX} \
  --bin-dir bin \
  --version ${INSTALLER_VER} \
  --destination ${PREFIX}

# rebuild rotarama and cablam caches
echo ""
echo "Rebuilding rotarama and cablam caches"
echo "====================================="
${PREFIX}/bin/mmtbx.rebuild_rotarama_cache
${PREFIX}/bin/mmtbx.rebuild_cablam_cache

echo ""
