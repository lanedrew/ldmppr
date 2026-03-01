# Package startup hook.
#
# Avoid changing process-wide preferences (such as OMP_THREAD_LIMIT) at load time.
.onLoad <- function(libname, pkgname) {
  invisible(NULL)
}
