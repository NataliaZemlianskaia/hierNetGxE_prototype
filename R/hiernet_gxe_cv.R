library(parallel)
library(dplyr)

hierNet.gxe.cv = function(G, E, Y, GxE, nFolds=5, ncores=5,
                          grid, max_iter=10000, tol=1e-3,
                          target_b_gxe_non_zero=NULL, 
                          working_set_min_size=100,
                          seed=2020){
  
  n = nrow(G)
  set.seed(seed)
  ncores = min(ncores, nFolds)
  folds = rep_len(1:nFolds, n)
  folds = sample(folds, n)
  
  df = do.call(rbind, mclapply(
    1:nFolds,
    function(i, ...) {
      fold = which(folds == i)
      
      G_train = G[-fold,]
      E_train = E[-fold]
      GxE_train = GxE[-fold,]
      Y_train = Y[-fold]
      
      G_valid = G[fold,]
      E_valid = E[fold]
      GxE_valid = GxE[fold,]
      Y_valid = Y[fold]
      
      fit = hierNet.gxe.fit(G_train, E_train, Y_train, GxE_train,
                            G_valid, E_valid, Y_valid, GxE_valid,
                            grid=grid, max_iter=max_iter, tol=tol,
                            target_lambdas=NULL,
                            working_set_min_size=working_set_min_size)
        
      path = fit$path
      path$fold = rep(i, nrow(path))
      return(path)
    },
    mc.cores=ncores
  ))
  
  if (is.null(target_b_gxe_non_zero)) {
    best_lambdas = df %>% 
      group_by(lambda_1,lambda_2 ) %>%
      summarize_all(mean) %>%
      ungroup() %>%
      slice(which.min(valid_loss))
  } else {
    best_lambdas = df %>%
      group_by(lambda_1,lambda_2 ) %>%
      summarize_all(mean) %>%
      filter(b_gxe_non_zero >= target_b_gxe_non_zero) %>% # && b_gxe_non_zero <= 10 * target_b_gxe_non_zero
      ungroup() %>%
      slice(which.min(valid_loss))
  }
  
  cat('-- Best lambdas lambda_1=', best_lambdas$lambda_1, ', lambda_2=', best_lambdas$lambda_2, "\n")

  fit = hierNet.gxe.fit(G, E, Y, GxE, 
                        NULL, NULL, NULL, NULL,
                        grid=grid, max_iter=max_iter, tol=tol,
                        target_lambdas=best_lambdas,
                        working_set_min_size=working_set_min_size)
    
  return(list(path=df, lambda_1=fit$lambda_1, lambda_2=fit$lambda_2,
              train_loss=fit$train_loss, result=fit$result, rules_stats=fit$rules_stat))
}
