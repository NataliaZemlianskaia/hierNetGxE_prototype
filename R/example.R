library(ggplot2)

seed = 2020
grid = NULL
#grid = 10^seq(-4, log10(1), length.out = 40)
grid_size = 40
p = 1000
sample_size = 200
family = "gaussian"
mode = "strong_hierarchical"
n_g_non_zero = 15
n_gxe_non_zero = 10
normalize = TRUE

max_iter = 10000
tol = 1e-4
working_set_min_size = 100
target_lambdas = NULL

dataset = data.gen(seed, sample_size, p, n_g_non_zero, n_gxe_non_zero,
                   family=family, normalize=normalize, mode=mode)

start = Sys.time()
fit = hierNet.gxe.fit(dataset$G_train, dataset$E_train, dataset$Y_train, dataset$GxE_train, 
                      dataset$G_valid, dataset$E_valid, dataset$Y_valid, dataset$GxE_valid,
                      grid=grid, grid_size=grid_size, max_iter=max_iter, tol=tol,
                      target_lambdas=NULL, working_set_min_size=working_set_min_size)
stop_cd = Sys.time() - start; stop_cd

df = fit$rules_stat
df$lambda_1_factor = factor(df$lambda_1)
df$lambda_2_factor = factor(df$lambda_2)

ggplot(df, aes(lambda_1_factor, lambda_2_factor, fill=100.0 * working_set_size / p)) + 
  scale_fill_distiller(palette = "RdBu") +
  scale_x_discrete("lambda_1", breaks=c(1)) +
  scale_y_discrete("lambda_2", breaks=c(1)) +
  labs(fill='Working Set Size %') +
  geom_tile()


start = Sys.time()
cv = hierNet.gxe.cv(dataset$G_train, dataset$E_train, dataset$Y_train, dataset$GxE_train, 
                    nFolds=5, ncores=5,
                    grid=grid, grid_size=grid_size, max_iter=max_iter, tol=tol, 
                    working_set_min_size=working_set_min_size,
                    target_b_gxe_non_zero=NULL, 
                    seed=seed)
stop_cd = Sys.time() - start; stop_cd

avg_valid_loss = cv$path %>% 
  group_by(lambda_1, lambda_2) %>%
  summarize_all(mean) %>%
  ungroup()
avg_valid_loss$lambda_1_factor = factor(avg_valid_loss$lambda_1)
avg_valid_loss$lambda_2_factor = factor(avg_valid_loss$lambda_2)

ggplot(avg_valid_loss, aes(lambda_1_factor, lambda_2_factor, fill=valid_loss)) + 
  scale_fill_distiller(palette = "RdBu") +
  scale_x_discrete("lambda_1", breaks=c()) +
  scale_y_discrete("lambda_2", breaks=c()) +
  labs(fill='Average Valid Loss') +
  geom_tile() +
  geom_point(aes(x=as.factor(cv$lambda_1), y=as.factor(cv$lambda_2)), colour="green")

selection.metrics(dataset, cv$result)
