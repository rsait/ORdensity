#' @title 
#' IindexRobust
#'
#' @description
#' IindexRobust
#' 
#' @name IindexRobust
#' 
#' @author Itziar Irigoien, Concepcion Arenas, Jose Maria Martinez-Otzeta
#' 
#' @export

IindexRobust <- function(di, vg){
# Depth RObust index of i in C: IRobust(i) = [1 + proxiRobust(i,C)/2VRobust(C)]^{-1}
# Input:
#      di: distances from individual i to the rest 
#      vg: half geometric variability of C, estimated by means of medians
# Output:
#       index: IRobust(i)
#--------------------------------------------------------------------
  
  # di <- di[di > 0]
#   phi2 <- median(di)^2 - vg  #median(di^2) - vg  #median(di)^2 - vg  Estas dos cosas no tienen porque coincidir para n par!!!
   I <- vg/Rfast::med(di^2)
   I <- as.numeric(I)

   return(I)
   
 }
