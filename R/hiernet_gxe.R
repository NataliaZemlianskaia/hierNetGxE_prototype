library(Rcpp)
library(RcppEigen)

sourceCpp("/home/zemlians@PREVMED.USC.EDU/FIGI/hierNetGxE/src/hiernet_gxe_bcd.cpp",
          verbose=TRUE, rebuild=TRUE)
#sourceCpp("/Users/nataliazemlianskaia/Desktop/hierNetGxE/src/hiernet_gxe_bcd.cpp",
#          verbose=TRUE, rebuild=TRUE)

linear.predictor = function(G, E, b_0, beta_G, beta_E, beta_GxE){
  f = b_0 + G %*% beta_G + beta_E * E + (G * E) %*% beta_GxE
  return(f)
}


hierNet.gxe.fit = function(G, E, Y, GxE, 
                           G_valid, E_valid, Y_valid, GxE_valid,
                           grid=NULL, grid_size=40, max_iter=10000, tol=1e-3,
                           target_lambdas=NULL,
                           working_set_min_size=100){
  active_set_tol = tol
  n = dim(G)[1]
  p = dim(G)[2]
  path = data.frame(lambda_1=double(),
                    lambda_2=double(),
                    valid_loss=double(),
                    train_loss=double(),
                    b_g_non_zero=integer(),
                    b_gxe_non_zero=integer())
  rules_stat = data.frame(
    iteration=integer(),
    safe_g=integer(),
    safe_gxe=integer(),
    working_set_size=integer(),
    b_g_nonzero=integer(),
    b_gxe_nonzero=integer(),
    lambda_1=double(),
    lambda_2=double(),
    stringsAsFactors=FALSE)   
  
  target_result = NULL
  target_train_loss = NULL
  target_lambda_1 = NULL
  target_lambda_2 = NULL
  
  set.seed(1)
  
  sum_E = sum(E)
  norm2_E = (E %*% E)[1,1]
  denominator_E = n * norm2_E - sum_E^2
  norm2_GxE = colSums(GxE^2); length(norm2_GxE)
  norm_GxE = sqrt(norm2_GxE)
  norm2_G = colSums(G^2); length(norm2_G)
  norm_G = sqrt(norm2_G)
  G_by_GxE = colSums(G * GxE); length(G_by_GxE)
  n_1 = 1 / n
  norm2_GxE_inv = 1 / norm2_GxE
  norm2_G_div_n = n_1 * norm2_G
  norm2_GxE_div_n = n_1 * norm2_GxE
  G_by_GxE_div_n = n_1 * G_by_GxE
  norm2_Y_div_n2 = (Y %*% Y)[1,1] / (n^2)
  b_gxe_denumerator = n_1 * norm2_GxE
  b_gxe_denumerator_inv = 1 / b_gxe_denumerator
  case1_detA = norm2_G * norm2_GxE - G_by_GxE * G_by_GxE
  case1_A22_div_detA = norm2_GxE / case1_detA
  case1_A12_div_detA = G_by_GxE / case1_detA
  case_3_A = norm2_G_div_n + b_gxe_denumerator
  case_3_B = 2 * G_by_GxE_div_n
  case5_detA = norm2_G_div_n * b_gxe_denumerator - G_by_GxE_div_n * G_by_GxE_div_n
  case5_A22_div_detA = b_gxe_denumerator / case5_detA
  case5_A12_div_detA = G_by_GxE_div_n / case5_detA
  active_set = as.integer(rep(1, p))
  b_0 = 0; b_e = 0
  b_g = rep(0, p)
  b_gxe = rep(0, p) 
  delta = rep(0, p)
  res = Y
  best_valid_loss = Inf
  GxE_by_Yn_abs = abs(Y %*% GxE)[1,] / n
  G_by_Yn_abs = abs(Y %*% G)[1,] / n

  if (is.null(grid)) {
    lambda_max = max(c(G_by_Yn_abs, GxE_by_Yn_abs))
    lambda_min = 1e-4 * lambda_max
    grid = 10^seq(log10(lambda_min), log10(lambda_max), length.out=grid_size)
  } 
  grid = sort(grid, decreasing=TRUE)
  snake = TRUE  
  
  abs_res_by_G_are_uptodate = FALSE
  lambda_iter = 0
  
  for (lambda_1 in grid){
    if (snake) {
      lambda_2_grid = grid
    } else {
      lambda_2_grid = sort(grid)
    }
    snake = !(snake)
    for (lambda_2 in lambda_2_grid){
      lambda_iter = lambda_iter + 1
      inner_nu = NULL
      SAFE_set_g = rep(TRUE, p)
      SAFE_set_gxe = rep(TRUE, p)
      working_set_size = sum((b_g != 0) | (b_gxe != 0))
      if (working_set_size <= 0) {
        working_set_size = working_set_min_size
      }
      
      dual_objective = -Inf
      nu = NULL      
      
      lambda_2_plus_1 = lambda_2 + lambda_1
      case_3_E = G_by_GxE_div_n * (lambda_1 - lambda_2)
      case_3_F = lambda_1 * b_gxe_denumerator - lambda_2 * norm2_G_div_n
      
      for (i_outter in 1:max_iter) {
        res = res + b_0 + E*b_e; b_0_old = b_0; b_e_old = b_e; sum_res = sum(res)
        b_e = (n * (E %*% res)[1,1] - sum_E * sum_res) / denominator_E
        b_0 = (sum_res - sum_E * b_e) / n
        res = res - b_0 - E*b_e
        if (b_0_old != b_0 || b_e_old != b_e) {
          abs_res_by_G_are_uptodate = FALSE
        }        
        
        nu_prev = res / n
        if (!abs_res_by_G_are_uptodate) {
          abs_res_by_G_div_n = absVectorByMatrix(n, res, G)
          abs_res_by_GxE_div_n = absVectorByMatrix(n, res, GxE)
          abs_res_by_G_are_uptodate = TRUE
        }    
        a = abs_res_by_G_div_n
        b = abs_res_by_GxE_div_n
        delta_safe_rule = ComputeDeltaSafeRule(lambda_1, lambda_2, abs_res_by_G_div_n, abs_res_by_GxE_div_n)
        b_tmp = (lambda_2 + delta_safe_rule)/b
        b_tmp[is.nan(b_tmp)] = Inf
        a_tmp = (lambda_1 - delta_safe_rule)/a
        a_tmp[is.nan(a_tmp)] = Inf      
        M = min(b_tmp, a_tmp)
        x_hat = ((Y %*% nu_prev)/(n * nu_prev %*% nu_prev))[1,1]
        x_opt = ifelse(abs(x_hat) <= M, x_hat, sign(x_hat) * M)
        nu_res = x_opt * nu_prev 
        dual_objective_res = (n/2) * (norm2_Y_div_n2 - sum((Y / n - nu_res)^2))
        
        if (dual_objective_res > dual_objective) {
          dual_objective = dual_objective_res
          nu = nu_res
        } else {
          a = absVectorByMatrix(1, nu, G)
          b = absVectorByMatrix(1, nu, GxE)
          x_opt = 1
        }
        
        primal_objective = (res %*% res)[1,1] / (2 * n) + lambda_1 * sum(pmax(abs(b_g), abs(b_gxe))) + lambda_2 * sum(abs(b_gxe))
        dual_gap = primal_objective - dual_objective  
        if (dual_gap < tol) {
          break
        }        

        r = sqrt(2 * dual_gap / n)
        nu_by_GxE = abs(x_opt) * b + r * norm_GxE
        nu_by_G = abs(x_opt) * a + r * norm_G
        d_j = (lambda_1 - lambda_2 - abs(x_opt) * a - pmax(lambda_2 - r * norm_GxE, abs(x_opt) * b)) / (norm_GxE + norm_G)
        
        safe_set_zero = pmax(0, nu_by_GxE - lambda_2) < lambda_1 - nu_by_G
        SAFE_set_gxe = SAFE_set_gxe & !safe_set_zero & (nu_by_GxE >= lambda_2)
        SAFE_set_g = SAFE_set_g & !safe_set_zero & ((nu_by_G >= lambda_1) | SAFE_set_gxe)
        
        update_to_zero_b_gxe = (1:p)[(b_gxe != 0) & (!SAFE_set_gxe)]
        if (length(update_to_zero_b_gxe) > 0) {
          if (length(update_to_zero_b_gxe) > 1) {
            res = res + (GxE[,update_to_zero_b_gxe] %*% b_gxe[update_to_zero_b_gxe])[,1]
          } else {
            res = res + GxE[,update_to_zero_b_gxe] * b_gxe[update_to_zero_b_gxe]
          }        
          b_gxe[update_to_zero_b_gxe] = 0     
          abs_res_by_G_are_uptodate = FALSE
        }
        update_to_zero_b_g = (1:p)[(b_g != 0) & (!SAFE_set_g)]
        if (length(update_to_zero_b_g) > 0) {
          if (length(update_to_zero_b_g) > 1) {
            res = res + (G[,update_to_zero_b_g] %*% b_g[update_to_zero_b_g])[,1]
          } else {
            res = res + G[,update_to_zero_b_g] * b_g[update_to_zero_b_g]
          }        
          b_g[update_to_zero_b_g] = 0          
          abs_res_by_G_are_uptodate = FALSE
        }
        
        if (i_outter >= 2) {
          working_set_size = 2 * working_set_size
        }
        d_j[(b_g != 0) | (b_gxe != 0)] = -Inf
        d_j[!SAFE_set_g] = Inf
        working_set = sort(d_j, index.return=TRUE)$ix
        working_set = working_set[SAFE_set_g[working_set]]
        if (length(working_set) <= working_set_size) {
          working_set_size = length(working_set)
        } else {
          working_set = working_set[1:working_set_size]
        }

        inner_dual_objective = -Inf
        inner_nu = NULL

        G_working_set = G[,working_set]
        GxE_working_set = GxE[,working_set]

        not_working_set = setdiff(1:p, working_set)
        stopifnot(sum(abs(b_g[not_working_set])) == 0)
        stopifnot(sum(abs(b_gxe[not_working_set])) == 0)

        abs_res_by_G_div_n_working_set = rep(0, working_set_size)
        abs_res_by_GxE_div_n_working_set = rep(0, working_set_size)
        
        active_set = as.integer(rep(0, p))
        
        current_active_set_tol = active_set_tol

        for (i_inner in 1:max_iter) {
          res = res + b_0 + E*b_e; b_0_old = b_0; b_e_old = b_e; sum_res = sum(res)
          b_e = (n * (E %*% res)[1,1] - sum_E * sum_res) / denominator_E
          b_0 = (sum_res - sum_E * b_e) / n
          res = res - b_0 - E*b_e
          b_0_vector = c(b_0); b_e_vector = c(b_e)
          max_diff = max((b_0_old - b_0)^2, n_1 * norm2_E * ((b_e_old - b_e)^2)) 
          
          if (b_0_old != b_0 || b_e_old != b_e) {
            abs_res_by_G_are_uptodate = FALSE
          }
          
          inner_nu_prev = res / n
         
          abs_res_by_G_div_n_working_set = absVectorByMatrix(n, res, G_working_set)
          abs_res_by_GxE_div_n_working_set = absVectorByMatrix(n, res, GxE_working_set)
          delta_safe_rule = ComputeDeltaSafeRule(lambda_1, lambda_2, abs_res_by_G_div_n_working_set, abs_res_by_GxE_div_n_working_set)
          
          b_tmp = (lambda_2 + delta_safe_rule)/abs_res_by_GxE_div_n_working_set
          b_tmp[is.nan(b_tmp)] = Inf
          a_tmp = (lambda_1 - delta_safe_rule)/abs_res_by_G_div_n_working_set
          a_tmp[is.nan(a_tmp)] = Inf      
          M = min(b_tmp, a_tmp)
          x_hat = ((Y %*% inner_nu_prev)/(n * inner_nu_prev %*% inner_nu_prev))[1,1]
          x_opt = ifelse(abs(x_hat) <= M, x_hat, sign(x_hat) * M)
          inner_nu_res = x_opt * inner_nu_prev
          inner_dual_objective_res = (n/2) * (norm2_Y_div_n2 - sum((Y / n - inner_nu_res)^2))

          if (inner_dual_objective_res > inner_dual_objective) {
            inner_dual_objective = inner_dual_objective_res
            inner_nu = inner_nu_res
          }
          
          primal_objective = (res %*% res)[1,1] / (2 * n) + lambda_1 * sum(pmax(abs(b_g[working_set]), abs(b_gxe[working_set]))) + lambda_2 * sum(abs(b_gxe[working_set]))
          inner_dual_gap = primal_objective - inner_dual_objective
          if (current_active_set_tol < 1e-29) {
            browser()
          }

          if (inner_dual_gap < tol) {
            break
          } else {
            if (i_inner > 1) {
              current_active_set_tol = current_active_set_tol / 10
            }
          }
          
          current_diff = CoordDescendStep(
            n, lambda_1, lambda_2,
            G, E, GxE, res,
            G_by_GxE, norm2_G, norm2_E, norm2_GxE,
            sum_E, denominator_E,
            case1_A12_div_detA, case1_A22_div_detA,
            case_3_A, case_3_E, case_3_F,
            case5_A12_div_detA, case5_A22_div_detA,
            SAFE_set_gxe, 
            b_0_vector, b_e_vector, b_g, b_gxe, delta,
            working_set, active_set, current_active_set_tol)
          max_diff = max(max_diff, current_diff)
          b_0 = b_0_vector[1]; b_e = b_e_vector[1]; 

          abs_res_by_G_are_uptodate = abs_res_by_G_are_uptodate && (max_diff > 0)
        }
        if (i_inner >= max_iter) {
          cat("-- WARNING: Reached max (inner) iterations for lambda_1=", lambda_1, ", lambda_2=", lambda_2, "\n")
        }
      }
      rules_stat[nrow(rules_stat) + 1,] = list(lambda_iter,
                                               sum(!SAFE_set_g),
                                               sum(!SAFE_set_gxe),
                                               working_set_size,
                                               sum(b_g != 0),
                                               sum(b_gxe != 0),
                                               lambda_1,
                                               lambda_2)
      if (i_outter >= max_iter) {
        cat("-- WARNING: Reached max (outter) iterations for lambda_1=", lambda_1, ", lambda_2=", lambda_2, "\n")
      }
      result = list(beta_0_hat=b_0, beta_G_hat=b_g, beta_E_hat=b_e, beta_GxE_hat=b_gxe)
      if (!is.null(Y_valid)) {
        valid_loss = mean((Y_valid - linear.predictor(G_valid, E_valid, 
                                                      b_0, b_g, b_e, b_gxe))^2) / 2
      } else {
        valid_loss = 0
      }
      
      train_loss = (res %*% res)[1,1] / (2 * n) + lambda_1 * sum(pmax(abs(b_g), abs(b_gxe))) + lambda_2 * sum(abs(b_gxe))
      path[nrow(path) + 1,] = list(lambda_1=lambda_1,
                                   lambda_2=lambda_2,
                                   valid_loss=valid_loss,
                                   train_loss=train_loss,
                                   b_g_non_zero=sum(b_g != 0),
                                   b_gxe_non_zero=sum(b_gxe != 0))
      
      if (!is.null(target_lambdas)) {
        if (lambda_1 == target_lambdas$lambda_1 && lambda_2 == target_lambdas$lambda_2) {
          target_result = result
          target_train_loss = train_loss
          target_lambda_1 = lambda_1
          target_lambda_2 = lambda_2
          return(list(path=path, result=target_result, train_loss=target_train_loss,
                      lambda_1=target_lambda_1, lambda_2=target_lambda_2,
                      rules_stat=rules_stat))
        }
      }
    }
  }
  return(list(path=path, result=target_result, train_loss=target_train_loss,
              lambda_1=target_lambda_1, lambda_2=target_lambda_2, rules_stat=rules_stat))
}
