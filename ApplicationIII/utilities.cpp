#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

////////////////////// BC function //////////////////////
// [[Rcpp::export]]
Rcpp::List bc_fun_cpp(const arma::vec &p_value, double alpha = 0.05)
{
  // initialization
  int n = p_value.n_elem;
  arma::uvec index_select;
  double e_tmp = 0;
  double thres = 0;
  double n_num = 1.0;
  arma::vec evalue = arma::zeros(n);
  // order p-value
  double p_tmp, n_dom, hat_fdp;
  arma::vec ordered_p = arma::sort(p_value, "descend");
  // BC algorithm
  for (int iter = 0; iter <= n; ++iter)
  {
    if (iter == n)
    {
      index_select.reset();
    }
    else
    {
      p_tmp = ordered_p(iter);
      n_num = 1 + arma::sum(p_value >= (1 - p_tmp));
      n_dom = n - iter;
      hat_fdp = (double)n_num / n_dom;
      if (hat_fdp <= alpha)
      {
        thres = p_tmp;
        index_select = arma::find(p_value <= p_tmp);
        e_tmp = (double)n / n_num;
        evalue.elem(index_select).fill(e_tmp);
        break;
      }
    }
  }
  // return results
  return Rcpp::List::create(Named("evalue") = evalue, Named("e_tmp") = e_tmp,
                            Named("thres") = thres, Named("n_num") = n_num,
                            Named("index_select") = index_select);
}

////////////////////// function to be optim //////////////////////
// [[Rcpp::export]]
double func_optim(arma::vec x, int len1, int len2, const arma::vec &q0,
                  const arma::mat &pi0_var, const arma::mat &f1_var,
                  const arma::vec &lpval)
{
  arma::vec theta = x.subvec(0, len1 - 1);
  arma::vec beta = x.subvec(len1, len1 + len2 - 1);
  arma::vec q1 = 1 - q0;
  arma::vec exp_eta_theta = exp(pi0_var * theta);
  arma::vec pi0_temp = exp_eta_theta / (1 + exp_eta_theta);
  arma::vec exp_eta_beta = exp(f1_var * beta);
  arma::vec k_temp = exp_eta_beta / (1 + exp_eta_beta);
  double res = -arma::sum(q0 % log(pi0_temp) + q1 % log(1 - pi0_temp) +
                          q1 % (-k_temp % lpval + log(1 - k_temp)));
  return res;
}

////////////////////// calculate weight //////////////////////
// [[Rcpp::export]]
Rcpp::List calculate_weight(const arma::vec &pvals, const arma::mat &pi0_var,
                            const arma::mat &f1_var, double pi0_init, double k_init,
                            double tol, double iterlim, double nlm_iter)
{
  // initialization
  int len_pi0 = pi0_var.n_cols;
  int len_f1 = f1_var.n_cols;
  int iter = 0;
  double pi0_est = pi0_init;
  arma::vec lpval = log(pvals);
  // No covariate effect - initialization
  arma::vec theta1 = arma::zeros(len_pi0);
  theta1(0) = log(pi0_est / (1 - pi0_est));
  arma::vec beta1 = arma::zeros(len_f1);
  beta1(0) = log(k_init / (1 - k_init));
  // EM iteration
  Environment stats = Environment::namespace_env("stats");
  Function nlm = stats["nlm"];
  Rcpp::Function func_optim("func_optim");
  double loglik0 = 0;
  double loglik1 = 0;
  arma::vec exp_eta_theta, exp_eta_beta, pi0, k, f1, f01, q0, st;
  arma::vec theta2, beta2;
  Rcpp::List obj;
  while (iter < iterlim)
  {
    // E-step
    exp_eta_theta = exp(pi0_var * theta1);
    exp_eta_beta = exp(f1_var * beta1);
    pi0 = exp_eta_theta / (1 + exp_eta_theta);
    k = exp_eta_beta / (1 + exp_eta_beta);
    f1 = (1 - k) % pow(pvals, -k);
    f01 = (1 - pi0) % f1 + pi0;
    q0 = pi0 / f01;
    // M step
    st = arma::join_cols(theta1, beta1);
    obj = nlm(Named("f") = func_optim, Named("p") = st, Named("len1") = len_pi0,
              Named("len2") = len_f1, Named("q0") = q0, Named("pi0_var") = pi0_var,
              Named("f1_var") = f1_var, Named("lpval") = lpval, Named("iterlim") = nlm_iter);

    // update new likelihood
    arma::vec est_tmp = obj["estimate"];
    theta2 = est_tmp.subvec(0, len_pi0 - 1);
    beta2 = est_tmp.subvec(len_pi0, len_pi0 + len_f1 - 1);
    loglik1 = arma::sum(log(f01));

    if (abs((loglik1 - loglik0) / loglik0) < tol)
    {
      break;
    }
    else
    {
      theta1 = theta2;
      beta1 = beta2;
      iter = iter + 1;
      loglik0 = loglik1;
    }
  }
  // return results
  return Rcpp::List::create(Named("pi0") = pi0, Named("k") = k,
                            Named("pi0_coef") = theta1, Named("k_coef") = beta1,
                            Named("loglik") = loglik1, Named("EM_iter") = iter);
}

//////////////////////  local FDR function //////////////////////
// [[Rcpp::export]]
arma::vec lfdr_fun_cpp(const arma::vec &pvalue, const arma::vec &kappa_fit,
                       const arma::vec &pi0_fit)
{
  arma::vec f1_fit = (1 - kappa_fit) % pow(pvalue, -kappa_fit);
  arma::vec lfdr_fit = pi0_fit / (pi0_fit + (1 - pi0_fit) % f1_fit);
  return lfdr_fit;
}

//////////////////////  GBC function //////////////////////
// [[Rcpp::export]]
Rcpp::List gbc_fun_cpp(const arma::vec &lfdr_fit, const arma::vec &lfdr_sym,
                       const arma::vec &kappa_fit, const arma::vec &pi0_fit,
                       double alpha_bc)
{
  // initialization
  int n = lfdr_fit.n_elem;
  arma::uvec index_select;
  double thres = 0;
  double num = 1.0;
  arma::vec evalue = arma::zeros(n);
  // order lfdr
  double lfdr_tmp, dom, hat_fdp, e_tmp;
  arma::vec ordered_lfdr = arma::sort(lfdr_fit, "descend");
  // BC algorithm
  for (int iter = 0; iter <= n; ++iter)
  {
    if (iter == n)
    {
      index_select.reset();
    }
    else
    {
      lfdr_tmp = ordered_lfdr(iter);
      num = 1 + arma::sum(lfdr_sym <= lfdr_tmp);
      dom = n - iter;
      hat_fdp = (double)num / dom;
      if (hat_fdp <= alpha_bc)
      {
        thres = lfdr_tmp;
        index_select = arma::find(lfdr_fit <= thres);
        e_tmp = (double)n / num;
        evalue.elem(index_select).fill(e_tmp);
        break;
      }
    }
  }
  // return results
  return Rcpp::List::create(Named("evalue") = evalue, Named("thres") = thres,
                            Named("num") = num, Named("index_select") = index_select);
}

////////////////////// model fit //////////////////////
// [[Rcpp::export]]
Rcpp::List modelfit_fun(const arma::vec &pvalue_train, const arma::mat &pi0_var_train,
                        const arma::mat &f1_var_train, const arma::vec &pvalue_test,
                        const arma::mat &pi0_var_test, const arma::mat &f1_var_test,
                        double alpha_bc, double pi0_init, double k_init, double tol,
                        double iterlim, double nlm_iter, bool weight_bool)
{
  // fit the model
  List res_fit = calculate_weight(pvalue_train, pi0_var_train, f1_var_train,
                                  pi0_init, k_init, tol, iterlim, nlm_iter);
  // pi0
  arma::vec pi0_coef = res_fit["pi0_coef"];
  arma::vec exp_eta_theta = exp(pi0_var_test * pi0_coef);
  arma::vec pi0_fit = exp_eta_theta / (1 + exp_eta_theta);
  pi0_fit.elem(arma::find(pi0_fit < 0.1)).fill(0.1);
  pi0_fit.elem(arma::find(pi0_fit > 1 - 1e-5)).fill(1 - 1e-5);

  // kappa
  arma::vec k_coef = res_fit["k_coef"];
  arma::vec exp_eta_beta = exp(f1_var_test * k_coef);
  arma::vec kappa_fit = exp_eta_beta / (1 + exp_eta_beta);

  // calculate lfdr
  arma::vec tilde_pvalue = 1 - pvalue_test;
  arma::vec lfdr_fit = lfdr_fun_cpp(pvalue_test, kappa_fit, pi0_fit);
  arma::vec lfdr_sym = lfdr_fun_cpp(tilde_pvalue, kappa_fit, pi0_fit);

  // GBH
  List gbc_fit = gbc_fun_cpp(lfdr_fit, lfdr_sym, kappa_fit, pi0_fit, alpha_bc);
  arma::vec evalue_vec = gbc_fit["evalue"];
  double weight_num = gbc_fit["num"];
  arma::uvec index_select = gbc_fit["index_select"];
  double thres = gbc_fit["thres"];
  double weight_dom;
  if (weight_bool) {
    // save lfdr for weight
    int n_test = pvalue_test.n_elem;
    arma::vec weight_dom_vec = arma::zeros(n_test);
    arma::vec lfdr_fit_tmp = arma::zeros(n_test);
    arma::vec lfdr_sym_tmp = arma::zeros(n_test);
    double thres_tmp;
    Rcpp::List gbc_fit_tmp;
    for (int iter_test = 0; iter_test < n_test; ++iter_test)
    {
      if (pvalue_test(iter_test) <= 0.5)
      {
        thres_tmp = thres;
      }
      else
      {
        // re-calculate LFDR
        lfdr_fit_tmp = lfdr_fit;
        lfdr_sym_tmp = lfdr_sym;
        lfdr_fit_tmp(iter_test) = lfdr_sym(iter_test);
        lfdr_sym_tmp(iter_test) = lfdr_fit(iter_test);
        gbc_fit_tmp = gbc_fun_cpp(lfdr_fit_tmp, lfdr_sym_tmp, kappa_fit, pi0_fit, alpha_bc);
        thres_tmp = gbc_fit_tmp["thres"];
      }
      if (lfdr_sym(iter_test) <= thres_tmp)
      {
        weight_dom_vec(iter_test) = 1;
      }
    }
    weight_dom = arma::sum(weight_dom_vec);
  } else {
    weight_dom = 0;
    weight_num = 1;
  }
  
  // return results
  return Rcpp::List::create(Named("weight_dom") = weight_dom,
                            Named("weight_num") = weight_num,
                            Named("evalue_vec") = evalue_vec,
                            Named("index_select") = index_select);
}

////////////////////////  get the index //////////////////////
// [[Rcpp::export]]
List get_train_test_indices(const List &split_groups, int iter_split)
{
  // Getting the test indices from the list
  arma::uvec index_test = split_groups[iter_split];
  // Getting the training indices by removing the test indices from the list
  arma::uvec index_train;
  for (int i = 0; i < split_groups.size(); ++i)
  {
    if (i != iter_split)
    {
      arma::uvec temp = split_groups[i];
      index_train.insert_rows(index_train.n_elem, temp);
    }
  }
  // return results
  return List::create(Named("index_test") = index_test, Named("index_train") = index_train);
}

////////////////////// cross fit //////////////////////
// [[Rcpp::export]]
arma::vec crossfit_fun(const arma::vec &pvalue_vec, const arma::mat &pi0_var,
                       const arma::mat &f1_var, double alpha_bc, double pi0_init,
                       double k_init, double tol, double iterlim, double nlm_iter,
                       int n_split, List split_groups, bool weight_bool)
{
  // Shuffle the index vector randomly
  int n = pvalue_vec.n_elem;
  List evalue_vec_list(n_split);
  arma::vec weight_num_vec = arma::zeros(n_split);
  arma::vec weight_dom_vec = arma::zeros(n_split);
  arma::vec evalue_vec = arma::zeros(n);
  // cross fit
  arma::vec pvalue_test, pvalue_train;
  arma::mat pi0_var_test, f1_var_test, pi0_var_train, f1_var_train;
  List model_fit;
  for (int iter_split = 0; iter_split < n_split; ++iter_split)
  {
    // spllit index
    List index_list = get_train_test_indices(split_groups, iter_split);
    arma::uvec index_test = index_list["index_test"];
    arma::uvec index_train = index_list["index_train"];

    // train data
    pvalue_train = pvalue_vec.elem(index_train);
    pi0_var_train = pi0_var.rows(index_train);
    f1_var_train = f1_var.rows(index_train);

    // test data
    pvalue_test = pvalue_vec.elem(index_test);
    pi0_var_test = pi0_var.rows(index_test);
    f1_var_test = f1_var.rows(index_test);

    // fit the model
    model_fit = modelfit_fun(pvalue_train, pi0_var_train, f1_var_train,
                             pvalue_test, pi0_var_test, f1_var_test,
                             alpha_bc, pi0_init, k_init, tol, iterlim,
                             nlm_iter, weight_bool);
    evalue_vec_list[iter_split] = model_fit["evalue_vec"];
    weight_num_vec(iter_split) = model_fit["weight_num"];
    weight_dom_vec(iter_split) = model_fit["weight_dom"];
  }
  double weight_dom_sum = arma::sum(weight_dom_vec);
  // get e-value
  if (weight_bool) {
    for (int iter_split = 0; iter_split < n_split; ++iter_split)
    {
      arma::vec evalue_vec_tmp = evalue_vec_list[iter_split];
      if (arma::sum(evalue_vec_tmp) > 0)
      {
        arma::uvec index_train = split_groups(iter_split);
        int n_train = index_train.n_elem;
        double num_tmp = n * weight_num_vec(iter_split) / n_train;
        double dom_tmp = weight_num_vec(iter_split) + weight_dom_sum - weight_dom_vec(iter_split);
        double weight = num_tmp / dom_tmp;
        evalue_vec_tmp = weight * evalue_vec_tmp;
        evalue_vec.elem(index_train) = evalue_vec_tmp;
      }
    }
  } else {
    double weight = 1.0;
    for (int iter_split = 0; iter_split < n_split; ++iter_split)
    {
      arma::vec evalue_vec_tmp = evalue_vec_list[iter_split];
      if (arma::sum(evalue_vec_tmp) > 0)
      {
        arma::uvec index_train = split_groups(iter_split);
        evalue_vec_tmp = weight * evalue_vec_tmp;
        evalue_vec.elem(index_train) = evalue_vec_tmp;
      }
    }
  }
  return evalue_vec;
}

////////////////////// shuffle the index //////////////////////
// [[Rcpp::export]]
List shuffle_fun(int n, int n_split)
{
  // Generate a shuffled index vector
  arma::uvec shuffled_index = arma::randperm(n);
  // Split the shuffled index vector into groups
  List split_groups(n_split);
  int group_size = ceil(static_cast<double>(n) / n_split);

  for (int i = 0; i < n_split; ++i)
  {
    int start_idx = i * group_size;
    int end_idx = std::min(start_idx + group_size - 1, n - 1);
    split_groups[i] = shuffled_index.subvec(start_idx, end_idx);
  }

  return split_groups;
}

////////////////////// eBH method //////////////////////
// [[Rcpp::export]]
arma::uvec ebh_fun_cpp(const arma::vec &e_value, double alpha = 0.05)
{
  // initialization
  int n = e_value.n_elem;
  arma::vec ordered_evalue = arma::sort(e_value, "ascend");
  arma::uvec index_select;
  double e_tmp, thres;
  // eBH algorithm
  for (int iter = 0; iter <= n; ++iter)
  {
    if (iter == n)
    {
      index_select.reset();
    }
    else
    {
      e_tmp = ordered_evalue(iter);
      thres = (double)n / (alpha * (n - iter));
      if (e_tmp + 1e-10 >= thres)
      {
        index_select = arma::find(e_value >= e_tmp);
        break;
      }
    }
  }
  // return results
  return index_select;
}

////////////////////// ebh-gbc function //////////////////////
// [[Rcpp::export]]
List ebh_gbc_cpp(const arma::vec &pvalue_vec, const arma::mat &pi0_var, arma::mat &f1_var,
                 int n_try, int n_split, double alpha_ebh, double alpha_bc,
                 double pi0_init, double k_init, double tol, double iterlim,
                 double nlm_iter, bool weight_bool)
{
  // initialization
  int n = pvalue_vec.n_elem;
  arma::mat evalue_mat(n_try, n);
  arma::vec evalue_tmp(n);
  arma::vec evalue_ave(n);
  // random split
  List index_list(n_try);
  for (int iter_try = 0; iter_try < n_try; ++iter_try)
  {
    index_list(iter_try) = shuffle_fun(n, n_split);
  }
  // cross fit
  for (int iter_try = 0; iter_try < n_try; ++iter_try)
  {
    List split_groups = index_list(iter_try);
    evalue_tmp = crossfit_fun(pvalue_vec, pi0_var, f1_var, alpha_bc, pi0_init,
                              k_init, tol, iterlim, nlm_iter, n_split,
                              split_groups, weight_bool);
    evalue_mat.row(iter_try) = evalue_tmp.t();
  }
  // average
  evalue_ave = arma::mean(evalue_mat, 0).t();
  // eBH
  arma::uvec index_select = ebh_fun_cpp(evalue_ave, alpha_ebh);
  // return results
  return List::create(Named("evalue_ave") = evalue_ave,
                      Named("index_select") = index_select,
                      Named("evalue_mat") = evalue_mat);
}