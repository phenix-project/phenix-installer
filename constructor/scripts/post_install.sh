SCRIPT_LOCATION=`${PREFIX}/bin/libtbx.find_in_repositories libtbx/auto_build/conda_build/write_env_files.py`

${PREFIX}/bin/python ${SCRIPT_LOCATION} \
  --program ${INSTALLER_NAME} \
  --prefix ${PREFIX} \
  --bin-dir bin \
  --version ${INSTALLER_VER} \
  --destination ${PREFIX}
