// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]

#include <cmath>
#include <algorithm>
#include <vector>
#include <omp.h>
#include <Rcpp.h>
#include <RcppEigen.h>

inline double sqr(double x) {
  return x * x;
}

inline double soft_threshold(double x, double lambda) {
  if (x > lambda) {
    return x - lambda;
  }
  if (x < - lambda){
    return x + lambda;
  }
  return 0;
}

inline void UpdateRes(
  const Eigen::Map<Eigen::MatrixXd>& G,
  const Eigen::Map<Eigen::MatrixXd>& GxE,
  const Eigen::Map<Eigen::VectorXd>& b_g,
  const Eigen::Map<Eigen::VectorXd>& b_gxe,
  double sign,
  int index,
  Eigen::Map<Eigen::VectorXd> res) {
  if (b_gxe[index] != 0) {
    res += sign * b_gxe[index] * GxE.col(index);
  }
  if (b_g[index] != 0) {
    res += sign * b_g[index] * G.col(index);
  }
}

// [[Rcpp::export]]
Eigen::VectorXd absVectorByMatrix(double n, const Eigen::Map<Eigen::VectorXd> v, Eigen::Map<Eigen::MatrixXd> M){
  return (v.transpose() * M).cwiseAbs() / n;
}
  
  
// [[Rcpp::export]]
Eigen::VectorXd ComputeDeltaSafeRule(
    double lambda_1,
    double lambda_2, 
    const Eigen::Map<Eigen::VectorXd>& abs_res_by_G_div_n,
    const Eigen::Map<Eigen::VectorXd>& abs_res_by_GxE_div_n
) {
  return (lambda_1 * abs_res_by_GxE_div_n - lambda_2 * abs_res_by_G_div_n).cwiseQuotient(abs_res_by_GxE_div_n + abs_res_by_G_div_n).cwiseMax(0).cwiseMin(lambda_1);
}  
  
// [[Rcpp::export]]
double CoordDescendStep(
  double n,    
  double lambda_1,
  double lambda_2,
  const Eigen::Map<Eigen::MatrixXd>& G,
  const Eigen::Map<Eigen::VectorXd>& E,
  const Eigen::Map<Eigen::MatrixXd>& GxE,
  Eigen::Map<Eigen::VectorXd> res,
  const Eigen::Map<Eigen::VectorXd>& G_by_GxE,
  const Eigen::Map<Eigen::VectorXd>& norm2_G,
  const double norm2_E,
  const Eigen::Map<Eigen::VectorXd>& norm2_GxE,
  const double sum_E,
  const double denominator_E,
  const Eigen::Map<Eigen::VectorXd>& case1_A12_div_detA,
  const Eigen::Map<Eigen::VectorXd>& case1_A22_div_detA,
  const Eigen::Map<Eigen::VectorXd>& case_3_A,
  const Eigen::Map<Eigen::VectorXd>& case_3_E,
  const Eigen::Map<Eigen::VectorXd>& case_3_F,
  const Eigen::Map<Eigen::VectorXd>& case5_A12_div_detA,
  const Eigen::Map<Eigen::VectorXd>& case5_A22_div_detA,
  const Rcpp::LogicalVector& SAFE_set_gxe,
  Eigen::Map<Eigen::VectorXd> b_0,
  Eigen::Map<Eigen::VectorXd> b_e,
  Eigen::Map<Eigen::VectorXd> b_g,
  Eigen::Map<Eigen::VectorXd> b_gxe,
  Eigen::Map<Eigen::VectorXd> delta,
  const Eigen::Map<Eigen::VectorXi>& indices,
  Eigen::Map<Eigen::VectorXi> active_set,
  double active_set_tol
) {
  double max_diff, curr_diff;
  double sum_res;
  bool has_solved;
  int index;
  double b_0_old, b_e_old, b_g_old, b_gxe_old;
  bool active_set_iteration = false, first_iteration = true;
  
  while (true) {
    max_diff = 0;
    
    if (first_iteration) {
      first_iteration = false;
    } else {
      res += E * b_e[0];
      res = res.array() + b_0[0];
      b_0_old = b_0[0]; b_e_old = b_e[0];
      sum_res = res.sum();
      b_e[0] = (n * E.dot(res) - sum_E * sum_res) / denominator_E;
      b_0[0] = (sum_res - sum_E * b_e[0]) / n;
      res = res.array() - b_0[0];
      res -= E*b_e;
      curr_diff = std::max(sqr(b_0[0] - b_0_old), norm2_E * sqr(b_e_old - b_e[0]) / n);
      max_diff = std::max(max_diff, curr_diff);
    }
    
    for (int k = 0; k < indices.size(); ++k) {
      index = indices[k] - 1;
      if (active_set_iteration && !active_set[index]) {
        continue;
      }
      b_g_old = b_g[index];
      b_gxe_old = b_gxe[index];
      
      UpdateRes(G, GxE, b_g, b_gxe, 1, index, res);
      
      const double G_by_res = G.col(index).dot(res);
      const double GxE_by_res = GxE.col(index).dot(res);
      
      const double GxE_by_res_div_n = GxE_by_res / n;
      const double norm2_G_div_n = norm2_G[index] / n;
      const double norm2_GxE_div_n = norm2_GxE[index] / n;
      
      double delta_upperbound, delta_lowerbound;
      double b_g_new;
      
      if (norm2_GxE[index] == 0.0 || !SAFE_set_gxe[index]) {
        delta_upperbound = lambda_1 - std::abs(G_by_res / n);
        delta_lowerbound = std::max(-lambda_2 + std::abs(GxE_by_res_div_n), 0.0);
        if (delta_lowerbound <= delta_upperbound) {
          b_g[index] = 0; b_gxe[index] = 0; delta[index] = delta_upperbound;
          UpdateRes(G, GxE, b_g, b_gxe, -1, index, res);
          curr_diff = norm2_G_div_n * sqr(b_g_old);
          max_diff = std::max(max_diff, curr_diff);
          if (!active_set_iteration && !active_set[index]) {
            active_set[index] = int(curr_diff > 0);
          }
          continue;
        } else {
          b_g_new = soft_threshold(G_by_res / n, lambda_1) / norm2_G_div_n;
          b_g[index] = b_g_new; b_gxe[index] = 0; delta[index] = 0;
          UpdateRes(G, GxE, b_g, b_gxe, -1, index, res);
          curr_diff = norm2_G_div_n * sqr(b_g_new - b_g_old);
          max_diff = std::max(max_diff, curr_diff);
          if (!active_set_iteration && !active_set[index]) {
            active_set[index] = int(curr_diff > 0);
          }        
          continue;
        }
      }  
      
      // Case 2
      const double G_by_res_div_n = G_by_res / n;
      delta_upperbound = lambda_1 - std::abs(G_by_res_div_n);
      delta_lowerbound = std::max(-lambda_2 + std::abs(GxE_by_res_div_n), 0.0);
      if (delta_lowerbound <= delta_upperbound) {
        b_g[index] = 0; b_gxe[index] = 0; delta[index] = delta_upperbound;
        UpdateRes(G, GxE, b_g, b_gxe, -1, index, res);
        curr_diff = std::max(norm2_G_div_n * sqr(b_g_old), 
                             norm2_GxE_div_n * sqr(b_gxe_old));
        max_diff = std::max(max_diff, curr_diff);      
        if (!active_set_iteration && !active_set[index]) {
          active_set[index] = int(curr_diff > 0);
        }      
        continue;
      }
      
      const double plus_minus_one[] = {-1.0, 1.0};
      const double lambda_2_plus_1 = lambda_1 + lambda_2;  
      const double G_by_GxE_div_n = G_by_GxE[index] / n;
      has_solved = false;
      
      // Case 1
      const double case1_B1_A22_div_detA = G_by_res * case1_A22_div_detA[index];
      double case1_B2, s, root_1, root_2, b_gxe_numerator;
      for (int i = 0; i < 2; ++i) {
        s = plus_minus_one[i];
        case1_B2 = GxE_by_res - s * lambda_2_plus_1 * n;
        root_1 = case1_B1_A22_div_detA - case1_B2 * case1_A12_div_detA[index];
        root_2 = (case1_B2 - root_1 * G_by_GxE[index]) / norm2_GxE[index];
        if (std::abs(root_2) > std::abs(root_1)) {
          b_gxe_numerator = GxE_by_res_div_n - G_by_GxE_div_n * root_1;
          if (s * b_gxe_numerator > lambda_2_plus_1) {
            b_g[index] = root_1; b_gxe[index] = root_2; delta[index] = lambda_1;
            UpdateRes(G, GxE, b_g, b_gxe, -1, index, res);
            curr_diff = std::max(norm2_G_div_n * sqr(b_g_old - root_1), 
                                 norm2_GxE_div_n * sqr(b_gxe_old - root_2));
            max_diff = std::max(max_diff, curr_diff);              
            if (!active_set_iteration && !active_set[index]) {
              active_set[index] = int(curr_diff > 0);
            }          
            has_solved = true;
            break;
          }
        }
      }
      if (has_solved) {
        continue;
      }
      
      // Case 3
      //const double case_3_B_1 = G_by_res_div_n;
      //const double case_3_B_2 = GxE_by_res_div_n;
      const double case_3_C = GxE_by_res_div_n * G_by_GxE_div_n - G_by_res_div_n * norm2_GxE[index] / n;
      const double case_3_D = GxE_by_res_div_n * norm2_G_div_n - G_by_res_div_n * G_by_GxE_div_n;
      double s_g, s_gxe, case_3_E_D, case_3_C_F, case_3_B_s_g, root_3, b_g_numerator;
      for (int i = 0; i < 2; ++i) {
        s_g = plus_minus_one[i];
        case_3_E_D = s_g * case_3_E[index] + case_3_D;
        case_3_C_F = s_g * case_3_C + case_3_F[index];
        case_3_B_s_g = s_g * 2 * G_by_GxE_div_n;
        for (int j = 0; j < 2; ++j) {
          s_gxe = plus_minus_one[j];
          root_3 = (s_gxe * case_3_E_D + case_3_C_F) / (case_3_A[index] + s_gxe * case_3_B_s_g);
          if ((root_3 >= 0) && (root_3 < lambda_1)) {
            root_1 = (G_by_res_div_n - s_g * (lambda_1 - root_3)) / (norm2_G_div_n + s_g * s_gxe * G_by_GxE_div_n);
            root_2 = s_g * s_gxe * root_1;
            b_gxe_numerator = GxE_by_res_div_n - G_by_GxE_div_n * root_1;
            b_g_numerator = (G_by_res - root_2 * G_by_GxE[index]) / n;
            if ((s_gxe * b_gxe_numerator > lambda_2 + root_3) &&
                (s_g * b_g_numerator > lambda_1 - root_3)) {
              b_g[index] = root_1; b_gxe[index] = root_2; delta[index] = root_3;
              UpdateRes(G, GxE, b_g, b_gxe, -1, index, res);
              curr_diff = std::max(norm2_G_div_n * sqr(b_g_old - root_1), 
                                   norm2_GxE_div_n * sqr(b_gxe_old - root_2));
              max_diff = std::max(max_diff, curr_diff);               
              if (!active_set_iteration && !active_set[index]) {
                active_set[index] = int(curr_diff > 0);
              }            
              has_solved = true;
              break;
            }          
          }
        }
        if (has_solved) {
          break;
        }      
      }
      if (has_solved) {
        continue;
      }    
      
      // Case 4
      b_g_new = soft_threshold(G_by_res_div_n, lambda_1) / norm2_G_div_n;
      b_gxe_numerator = GxE_by_res_div_n - G_by_GxE_div_n * b_g_new;
      if (std::abs(b_gxe_numerator) <= lambda_2) {
        b_g[index] = b_g_new; b_gxe[index] = 0; delta[index] = 0;
        UpdateRes(G, GxE, b_g, b_gxe, -1, index, res);
        curr_diff = std::max(norm2_G_div_n * sqr(b_g_old - b_g_new), 
                             norm2_GxE_div_n * sqr(b_gxe_old));
        max_diff = std::max(max_diff, curr_diff);         
        if (!active_set_iteration && !active_set[index]) {
          active_set[index] = int(curr_diff > 0);
        }      
        continue;
      }    
      
      // Case 5
      double case5_B2;
      for (int i = 0; i < 2; ++i) {
        s_g = plus_minus_one[i];
        for (int j = 0; j < 2; ++j) {
          s_gxe = plus_minus_one[j];
          case5_B2 = GxE_by_res_div_n - s_gxe * lambda_2;
          root_1 = (G_by_res_div_n - s_g * lambda_1) * case5_A22_div_detA[index] - case5_B2 * case5_A12_div_detA[index];
          b_gxe_numerator = GxE_by_res_div_n - G_by_GxE_div_n * root_1;
          if  (s_gxe * b_gxe_numerator > lambda_2) {
            root_2 = n * (case5_B2 - root_1 * G_by_GxE_div_n) / norm2_GxE[index];
            b_g_numerator = (G_by_res - root_2 * G_by_GxE[index]) / n;
            if (s_g * b_g_numerator > lambda_1) {
              b_g[index] = root_1; b_gxe[index] = root_2; delta[index] = 0;
              UpdateRes(G, GxE, b_g, b_gxe, -1, index, res);
              curr_diff = std::max(norm2_G_div_n * sqr(b_g_old - root_1), 
                                   norm2_GxE_div_n * sqr(b_gxe_old - root_2));
              max_diff = std::max(max_diff, curr_diff);                     
              if (!active_set_iteration && !active_set[index]) {
                active_set[index] = int(curr_diff > 0);
              }            
              has_solved = true;
              break;
            }
          }
        }
        if (has_solved) {
          break;
        }      
      }
    }
    
    if (active_set_iteration) {
      if (max_diff < active_set_tol) {
        active_set_iteration = false;
      }      
    } else {
      if (max_diff < active_set_tol) {
        return max_diff;
      } else {
        active_set_iteration = true;
      }
    }
  }
}