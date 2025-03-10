#' Power analysis for a Cormack-Jolly-Seber Survival and Detection Model
#' @author Ryan N. Kinzer
#' @param mark_size 
#' @param phi1
#' @param p
#' @param iterations 
#'
#' @return
#' @export
#'
#' @examples
cjs_power_analysis <- function(mark_size, phi, p, iterations = 100){
  df <- NULL
  df_mark <- NULL
  for(n in 1:length(mark_size)){
    for(i in 1:iterations){
      dat <- simulate_cjs_data(mark_size[n], phi, p)
      dat.proc <- process.data(dat, model = "cjs", begin.time = 1)
      dat.ddl <- make.design.data(dat.proc)
      tryCatch({
        cjs.model <- crm(dat.proc, dat.ddl, hessian = TRUE, model.parameters = list(Phi = list(formula = ~time), p = list(formula = ~time)))
        tmp <- bind_rows(cjs.model$results$reals, .id = 'param') %>%
          mutate(n = mark_size[n],
                 iter = i)
        df_mark <- bind_rows(df_mark, tmp)
      },
      error = function(e){
        return(NULL)
      }
      )
    }
    df <- bind_rows(df, df_mark)
  }
  return(df)
}
