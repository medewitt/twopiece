########################################################################################################
#
# DEFINITION OF GENERIC FUNCTIONS
#
########################################################################################################

setGeneric("clusterprobs", function(fit, x, iter) standardGeneric("clusterprobs"))

setGeneric("fixLabelSwitch", function(fit,x,z,method='ECR') standardGeneric("fixLabelSwitch"))
