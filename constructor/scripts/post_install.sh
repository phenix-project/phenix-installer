# compile Python files
echo ""
echo "Compiling Python files"
echo "======================"
pushd ${PREFIX}/lib
${PREFIX}/bin/libtbx.py_compile_all -i
popd
echo "done"

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
if [ -f "${PREFIX}/${INSTALLER_NAME}_env.sh" ]; then
  echo "${PREFIX}/${INSTALLER_NAME}_env.sh exists"
fi
if [ -f "${PREFIX}/${INSTALLER_NAME}_env.csh" ]; then
  echo "${PREFIX}/${INSTALLER_NAME}_env.csh exists"
fi
echo "done"

# rebuild rotarama and cablam caches
echo ""
echo "Rebuilding rotarama and cablam caches"
echo "====================================="
${PREFIX}/bin/mmtbx.rebuild_rotarama_cache
${PREFIX}/bin/mmtbx.rebuild_cablam_cache
echo "done"

# build documentation
echo ""
echo "Building Phenix documentation"
echo "============================="
${PREFIX}/bin/phenix.rebuild_docs
echo "done"

# print copy of conclusion.txt without Windows text
echo "
==========================================================================

                      Phenix installation complete
                      ----------------------------

You can begin using Phenix now by setting your environment with the
'source' command:

  csh users:
    source ${PREFIX}/phenix_env.csh

  bash (zsh) users:
    source ${PREFIX}/phenix_env.sh

To use Phenix, go to a work directory and type:

  phenix

On linux and macOS, you may wish to put the appopriate source statement
in your shell startup script (e.g. .cshrc or .bashrc or .zshrc).

You have successfully installed Phenix.

==========================================================================
"

echo ""
