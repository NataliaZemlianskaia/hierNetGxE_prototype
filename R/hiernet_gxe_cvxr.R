library("CVXR")

hier_linear_regression_e = function(p, q, xx, y, lambda_1, lambda_3, penalty, hierarchy="strong"){
  stopifnot(hierarchy == "strong")
  n = dim(xx)[1]
  beta_0 = Variable(1)
  beta_x = Variable(p)
  beta_e = Variable(q)
  theta = Variable(p, q)
  
  linear_predictor = function(xx, beta_0, beta_x, beta_e, theta){
    f = beta_0 + xx[,1:p] %*% beta_x + xx[,(p + 1):(p + q)] %*% beta_e + xx[,(p + q + 1):(p*q + p + q)] %*% reshape_expr(t(theta), c(p*q, 1))
    return(f)
  }
  
  loss = ((sum((y - linear_predictor(xx, beta_0, beta_x, beta_e, theta))^2))) / (2 * n)
  
  #if (hierarchy == "strong"){
  #  constraint1 = list(theta == t(theta), diag(theta) == 0)
  #}
  #if (hierarchy == "weak"){
  #  constraint1 = list(diag(theta) == 0)
  #}
  
  lasso_reg = function(theta, lambda) {
    lasso =  p_norm(theta, 1)
    lambda * (lasso)
  }
  if (penalty == "reg"){
    obj = loss
  }
  if (penalty == "group_lasso"){
    group_reg = function(beta_x, beta_e, theta, lambda) {
      #group_ = sum(p_norm(hstack(beta, theta), p=2, axis=1))
      group_x = sum(p_norm(hstack(beta_x, theta), p=2, axis=1))
      group_e = sum(p_norm(vstack(t(beta_e), theta), p=2, axis=2))
      lambda * (group_x + group_e)
    }
    obj = loss + lasso_reg(theta, lambda_3) + group_reg(beta_x, beta_e, theta, lambda_1)# + lasso_reg(beta, lambda_beta)
  }
  if (penalty == "l_inf"){
    l_inf_reg = function(beta_x, beta_e, theta, lambda) {
      group_x = sum(max_entries(abs(hstack(beta_x, theta)), axis=1))
      group_e = sum(max_entries(abs(vstack(t(beta_e), theta)), axis=2))
      lambda * (group_x + group_e)
    }
    obj = loss + lasso_reg(theta, lambda_3) + l_inf_reg(beta_x, beta_e, theta, lambda_1)# + lasso_reg(beta, lambda_beta)
  }
  max_reg = function(beta_x, beta_e, theta, lambda) {
    #max_e = sum(max_elemwise(sum_entries(abs(theta), axis=2), abs(t(beta_e))))
    max_x = sum(max_elemwise(sum_entries(abs(theta), axis=1), abs(beta_x)))
    lambda * (max_x)# + max_e)
  }
  if (penalty == "hierNet"){
    #epsilon = (1e-8)*lambda_1
    #ridge_reg = function(beta_x, beta_e, theta, epsilon) {
    #  ridge = p_norm(beta_x, 2) + 0.5 * p_norm(theta, 2)
    #  0.5 * epsilon * (ridge)
    #}
    
    obj = loss + lasso_reg(theta, lambda_1) + max_reg(beta_x, beta_e, theta, lambda_1) #- lasso_reg(beta, lambda_beta)
  }
  if (penalty == "hierNet2"){
    #epsilon = (1e-8)*lambda_1
    #ridge_reg = function(beta, theta, epsilon){
    #  ridge = p_norm(beta, 2) + 0.5 * p_norm(theta, 2)
    #  0.5 * epsilon * (ridge)
    #}
    obj = loss + lasso_reg(theta, lambda_1) + max_reg(beta_x, beta_e, theta, lambda_3) #- lasso_reg(beta, lambda_beta)
  }
  if (penalty == "all_pairs_lasso"){
    obj = loss + lasso_reg(theta, lambda_1) + lasso_reg(beta_x, lambda_1) #+ lasso_reg(beta_e, lambda_1)
  }
  if (penalty == "l_1/l_2"){
    l12_reg = function(beta_x, beta_e, theta, lambda){
      l12_x = sum_entries(p_norm(vstack(beta_x, sum_entries(abs(theta), axis=1)), p=2, axis=2), axis=1)
      #l12_e = sum_entries(p_norm(hstack(t(beta_e), sum_entries(abs(theta), axis=2)), p=2, axis=1), axis=2)
      lambda * (l12_x)# + l12_e)
    }
    obj = loss + lasso_reg(theta, lambda_1) + l12_reg(beta_x, beta_e, theta, lambda_3)
  }
  if (penalty == "l_1/l_2_1"){
    l12_reg = function(beta_x, beta_e, theta, lambda){
      l12_x = sum_entries(p_norm(vstack(beta_x, sum_entries(abs(theta), axis=1)), p=2, axis=2), axis=1)
      #l12_e = sum_entries(p_norm(hstack(t(beta_e), sum_entries(abs(theta), axis=2)), p=2, axis=1), axis=2)
      lambda * (l12_x)# + l12_e)
    }
    obj = loss + lasso_reg(theta, lambda_1) + l12_reg(beta_x, beta_e, theta, lambda_1)
  }
  
  #prob = Problem(Minimize(obj))#, constraints = constraint1)
  constraint1 = list(beta_x[16] == 0.01353009, theta[16, 1] == 0)
  prob = Problem(Minimize(obj), constraints = constraint1)
  
  
  #result = solve(prob, "SCS", ignore_dcp=TRUE)
  #result = solve(prob, ignore_dcp=TRUE)
  result = solve(prob, ignore_dcp=TRUE, verbose=TRUE, abstol=1e-8, reltol=1e-8, feastol=1e-8, max_iter=20000)
  beta_x = result$getValue(beta_x)
  beta_e = result$getValue(beta_e)
  beta_0 = result$getValue(beta_0)
  beta_main = c(beta_x, beta_e)
  beta_interaction = result$getValue(theta)
  #if (hierarchy == "strong"){
  #  beta_interaction = result$getValue(theta)
  #}
  #if (hierarchy == "weak"){
  #  beta_interaction = (result$getValue(theta) + t(result$getValue(theta)))/2
  #}
  return(list(beta_0=beta_0, beta_main=beta_main, beta_x=beta_x, beta_e=beta_e, beta_interaction=beta_interaction, value=result$value))
}

cvxr_res = hier_linear_regression_e(p, 1,
                                    cbind(dataset$G_train, dataset$E_train, dataset$GxE_train),
                                    dataset$Y_train,
                                    target_lambdas$lambda_2, target_lambdas$lambda_1, penalty="hierNet2", hierarchy="strong")
