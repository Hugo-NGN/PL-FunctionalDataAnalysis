library(fda)

derive_fd <- function(fd_obj, fine_grid, order =1) {
  fd_obj_deriv1 <- eval.fd(fine_grid, fd_obj, Lfdobj = order)
  return(fd_obj_deriv)
}