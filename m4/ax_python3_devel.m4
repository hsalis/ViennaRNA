AC_DEFUN([AX_PYTHON3_DEVEL],[

    dnl Allow the user to override the interpreter path
    AC_ARG_VAR([PYTHON3], [Path to Python3 interpreter (e.g.: /usr/bin/python3)])
    AC_ARG_VAR([PYTHON3_DIR], [Directory to install python3 scripts in])
    AC_ARG_VAR([PYTHON3_EXECDIR], [Directory to install architecture dependent python3 things in])

    dnl (AM_PATH_PYTHON) cannot be used for multiple Python versions at once
    if test -z "$PYTHON3" ; then
      AC_PATH_PROGS([PYTHON3], [python3 python312 python3.12 python311 python3.11 python310 python3.10 python39 python3.9 python38 python3.8 python37 python3.7 python36 python3.6 python35 python3.5 python34 python3.4], [no])
    fi
    AC_SUBST([PYTHON3])

    dnl Initialize outputs so they are always defined
    PYTHON3_INC=""
    PYTHON3_LDFLAGS=""
    PYTHON3_LD=""
    PYTHON3_SO=""
    PYTHON3_CACHE_TAG=""
    PYTHON3_CACHE_OPT1_EXT=""
    PYTHON3_INCLUDE_DIR=""
    AC_SUBST([PYTHON3_INC])
    AC_SUBST([PYTHON3_LDFLAGS])
    AC_SUBST([PYTHON3_LD])
    AC_SUBST([PYTHON3_SO])
    AC_SUBST([PYTHON3_CACHE_TAG])
    AC_SUBST([PYTHON3_CACHE_OPT1_EXT])

    if test "x${PYTHON3}" != "xno" ; then
      if $PYTHON3 -c 'import sysconfig' >/dev/null 2>&1
      then
        AC_MSG_CHECKING([for Python3 include path])
        PYTHON3_INCLUDE_DIR=`$PYTHON3 -c "import sysconfig; print(sysconfig.get_path('include') or '')" 2>/dev/null`
        AC_MSG_RESULT([$PYTHON3_INCLUDE_DIR])

        dnl Keep PYTHON3_INC as compiler flags, but test Python.h using the raw directory path
        AC_MSG_CHECKING([for $PYTHON3_INCLUDE_DIR/Python.h])
        if test -n "$PYTHON3_INCLUDE_DIR" && test -f "$PYTHON3_INCLUDE_DIR/Python.h"
        then
          AC_MSG_RESULT([yes])
          PYTHON3_INC="-I$PYTHON3_INCLUDE_DIR"
          AC_SUBST([PYTHON3_INC])

          AC_MSG_CHECKING([for Python3 ldflags])
          case "$host" in
            dnl Handle OSX Python extensions differently
            dnl see: http://blog.tim-smith.us/2015/09/python-extension-modules-os-x/
            *-darwin* | *-macos10*)
              PYTHON3_LDFLAGS="-bundle -undefined dynamic_lookup"
              ;;
            *)
              PYTHON3_LDFLAGS=`$PYTHON3 -c "import sysconfig; \
libdir = sysconfig.get_config_var('LIBDIR') or ''; \
ldlibrary = sysconfig.get_config_var('LDLIBRARY') or ''; \
libs = sysconfig.get_config_var('LIBS') or ''; \
syslibs = sysconfig.get_config_var('SYSLIBS') or ''; \
linkforshared = sysconfig.get_config_var('LINKFORSHARED') or ''; \
parts = []; \
if linkforshared: parts.extend(linkforshared.split()); \
if libdir: parts.append('-L' + libdir); \
if ldlibrary.startswith('lib') and '.so' in ldlibrary: \
    parts.append('-l' + ldlibrary[3:].split('.so')[0]); \
elif ldlibrary.startswith('lib') and '.a' in ldlibrary: \
    parts.append('-l' + ldlibrary[3:].split('.a')[0]); \
elif ldlibrary.startswith('lib') and '.dylib' in ldlibrary: \
    parts.append('-l' + ldlibrary[3:].split('.dylib')[0]); \
parts.extend(libs.split()); \
parts.extend(syslibs.split()); \
print(' '.join(parts))" 2>/dev/null`
              ;;
          esac
          AC_SUBST([PYTHON3_LDFLAGS])
          AC_MSG_RESULT([$PYTHON3_LDFLAGS])

          AC_MSG_CHECKING([for Python3 extension linker])
          PYTHON3_LD=`$PYTHON3 -c "import sysconfig; print(sysconfig.get_config_var('LDSHARED') or '')" 2>/dev/null`
          AC_SUBST([PYTHON3_LD])
          AC_MSG_RESULT([$PYTHON3_LD])

          AC_MSG_CHECKING([for directory to install Python3 scripts in])
          if test -z "$PYTHON3_DIR" ; then
            dnl Prevent premature substitution of ${prefix}
            PYTHON3_DIR=`$PYTHON3 -c "import sysconfig; \
print(sysconfig.get_path('purelib', vars={'base':'$' '{prefix}'}))" 2>/dev/null`
          fi
          AC_SUBST([PYTHON3_DIR])
          AC_SUBST([python3dir], [$PYTHON3_DIR])
          AC_SUBST([python3dir_expanded], [$(eval printf "%s" ${python3dir})])
          AC_MSG_RESULT([$PYTHON3_DIR])

          AC_MSG_CHECKING([for directory to install architecture dependent python3 things in])
          if test -z "$PYTHON3_EXECDIR" ; then
            dnl Prevent premature substitution of ${exec_prefix}
            PYTHON3_EXECDIR=`$PYTHON3 -c "import sysconfig; \
print(sysconfig.get_path('platlib', vars={'platbase':'$' '{exec_prefix}'}))" 2>/dev/null`
          fi
          AC_SUBST([PYTHON3_EXECDIR])
          AC_SUBST([py3execdir], [$PYTHON3_EXECDIR])
          AC_SUBST([py3execdir_expanded], [$(eval printf "%s" ${py3execdir})])
          AC_MSG_RESULT([$PYTHON3_EXECDIR])

          AC_MSG_CHECKING([for Python3 module extension])
          dnl Usually ".so", but for example, Mac OS X may differ.
          PYTHON3_SO=`$PYTHON3 -c "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX') or '')" 2>/dev/null`
          AC_SUBST([PYTHON3_SO])
          AC_MSG_RESULT([$PYTHON3_SO])

          AC_MSG_CHECKING([for Python3 tag for cached compiled scripts])
          PYTHON3_CACHE_TAG=`$PYTHON3 -c "import sys; print(getattr(sys.implementation, 'cache_tag', ''))" 2>/dev/null`
          AC_SUBST([PYTHON3_CACHE_TAG])
          AC_MSG_RESULT([$PYTHON3_CACHE_TAG])

          AC_MSG_CHECKING([for Python3 extension of cached and optimized bytecode])
          PYTHON3_CACHE_OPT1_EXT=`$PYTHON3 -c "import importlib.util, sys; \
print('%s.pyo' % sys.implementation.cache_tag) if sys.version_info.minor < 5 \
else print('{1}{2}'.format(*importlib.util.cache_from_source('', optimization=1).rpartition(sys.implementation.cache_tag)))" 2>/dev/null`
          AC_SUBST([PYTHON3_CACHE_OPT1_EXT])
          AC_MSG_RESULT([$PYTHON3_CACHE_OPT1_EXT])

        else
          AC_MSG_RESULT([no])
          AC_MSG_WARN([
**********************************************************************
Python.h not found!
You probably need to install python-dev, python-devel, or something similar.
**********************************************************************
])
          python3_enabled_but_failed="Python.h missing"
        fi
      else
        AC_MSG_WARN([
**********************************************************************
Failed to import sysconfig!
**********************************************************************
])
        python3_enabled_but_failed="Can't import sysconfig"
      fi
    else
      python3_enabled_but_failed="python3 executable missing"
    fi
])