rm(list=ls())
load("data/example_data.rda")
data = example_data

mult_model = case ~ tumor_1 + tumor_2 + tumor_3* tumor_2
add_model = ~ geno_1 + geno_2
unobserved = 888

# mvpoly = function(data, ymodel = y ~ tumor, x_stuff = ~ x,
# unobserved = 888) {

col_vars = all.vars(mult_model)
sub_data = data[, col_vars]
sub_data = as.matrix(sub_data)
checks = all(sub_data %in% c(0, 1, NA, unobserved))
if (!checks) {
  stop(paste0("Not all data within 0, 1, NA, ", unobserved))
}
sub_data[ sub_data %in% unobserved] = NA
sub_data = as.data.frame(sub_data)
   mm = model.matrix(
     object = mult_model,
     data = sub_data)

   cn = colnames(mm)
   mm = mm[ , !cn %in% "(Intercept)"]
   bad = array(mm %in% unobserved, dim = dim(mm))
   keep = mm[rowSums(bad) == 0,]
   um = unique(keep)
   rownames(um) = NULL

   mm[ mm %in% unobserved] = NA
   delta0 = colMeans(mm, na.rm = TRUE)

   add_model = delete.response(add_model)
   add_data = model.matrix(object = add_model,
                            data = data)
   # model.frame may be better


   # }
