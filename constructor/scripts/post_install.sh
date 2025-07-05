# stop for any failure
set -e

# store current directory
current_dir=`pwd`

# unset problematic environment variables
unset LIBTBX_BUILD

# compile Python files
echo ""
echo "Compiling Python files"
echo "======================"
cd ${PREFIX}/lib
${PREFIX}/bin/libtbx.py_compile_all -i
cd ${current_dir}
echo "done"

# create custom bin directory
echo ""
echo "Creating phenix_bin directory"
echo "============================="
SCRIPT_LOCATION=`${PREFIX}/bin/libtbx.find_in_repositories libtbx/auto_build/conda_build/create_custom_bin.py`
${PREFIX}/bin/python ${SCRIPT_LOCATION} \
  --prefix ${PREFIX} \
  --custom_bin phenix_bin \
  --packages phenix dials dxtbx
echo "done"

# write phenix_env files
echo ""
echo "Writing phenix_env files"
echo "========================"
SCRIPT_LOCATION=`${PREFIX}/bin/libtbx.find_in_repositories libtbx/auto_build/conda_build/write_env_files.py`
${PREFIX}/bin/python ${SCRIPT_LOCATION} \
  --program ${INSTALLER_NAME} \
  --prefix ${PREFIX} \
  --bin-dir phenix_bin \
  --version ${INSTALLER_VER} \
  --destination ${PREFIX}
if [ -f "${PREFIX}/${INSTALLER_NAME}_env.sh" ]; then
  echo "${PREFIX}/${INSTALLER_NAME}_env.sh created"
fi
if [ -f "${PREFIX}/${INSTALLER_NAME}_env.csh" ]; then
  echo "${PREFIX}/${INSTALLER_NAME}_env.csh created"
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
${PREFIX}/bin/phenix_html.rebuild_docs
echo "done"

# build .app on macOS
if [[ "$OSTYPE" == "darwin"* ]]; then
  echo ""
  echo "Creating macOS phenix.app"
  echo "========================="
  cd ${PREFIX}
  EXTRA_LINES=`echo "os.environ['PHENIX'] = '${PREFIX}'\nos.environ['PHENIX_PREFIX'] = '${PREFIX}'\nos.environ['PHENIX_VERSION'] = '${INSTALLER_VER}'"`
  ${PREFIX}/bin/libtbx.create_mac_app \
    phenix \
    --app_name=${INSTALLER_NAME}-${INSTALLER_VER} \
    --dest=${PREFIX} \
    --alias_build \
    --extra_lines="${EXTRA_LINES}"
  cd ${current_dir}
  echo "done"
fi

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

On macOS, you can also start Phenix by using the app bundle in the
${PREFIX} directory.

You have successfully installed Phenix.

==========================================================================
"

echo ""
