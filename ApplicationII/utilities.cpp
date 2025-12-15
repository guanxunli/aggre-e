#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// BH function
// [[Rcpp::export]]
Rcpp::List bh_fun_cpp(const arma::vec &p_value, double alpha = 0.05)
{
  // initialization
  int n = p_value.n_elem;
  arma::uvec index_select;
  double e_tmp = 0;
  double tau = 0;
  arma::vec evalue = arma::zeros(n);
  // order p-value
  double p_tmp, thres;
  arma::vec ordered_p = arma::sort(p_value, "descend");
  // BH algorithm
  for (int iter = 0; iter <= n; ++iter)
  {
    if (iter == n)
    {
      index_select.reset();
    }
    else
    {
      p_tmp = ordered_p(iter);
      thres = alpha * (n - iter) / n;
      if (p_tmp <= thres)
      {
        tau = thres;
        index_select = arma::find(p_value <= tau);
        e_tmp = 1.0 / tau;
        evalue.elem(index_select).fill(e_tmp);
        break;
      }
    }
  }
  // return results
  return Rcpp::List::create(Named("evalue") = evalue, Named("e_tmp") = e_tmp,
                            Named("tau_thres") = tau,
                            Named("index_select") = index_select);
}

// Storey'smethod
// [[Rcpp::export]]
Rcpp::List storey_fun_cpp(const arma::vec &p_value, double alpha = 0.05,
                          double lambda= 0.5)
{
  // initialization
  int n = p_value.n_elem;
  arma::uvec index_select, index_select_lambda;
  double e_tmp = 0;
  double tau = 0;
  double pi_0 = 1;
  arma::vec evalue = arma::zeros(n);
  // estimate pi_0
  int nreject_lambda = arma::sum(p_value <= lambda);
  pi_0 = (1.0 + n - nreject_lambda) / ((1 - lambda) * n);
  // order p-value
  double p_tmp, thres;
  arma::vec ordered_p = arma::sort(p_value, "descend");
  // BH algorithm
  for (int iter = 0; iter <= n; ++iter)
  {
    if (iter == n)
    {
      index_select.reset();
    }
    else
    {
      p_tmp = ordered_p(iter);
      thres = alpha * (n - iter) / (n * pi_0);
      if (p_tmp <= thres)
      {
        tau = thres;
        index_select = arma::find(p_value <= tau);
        e_tmp = 1.0 / tau;
        evalue.elem(index_select).fill(e_tmp);
        break;
      }
    }
  }
  // return results
  return Rcpp::List::create(Named("evalue") = evalue, Named("e_tmp") = e_tmp,
                            Named("tau_thres") = tau, Named("pi_0") = pi_0,
                            Named("nreject_lambda") = nreject_lambda,
                            Named("index_select") = index_select);
}

// BC function
// [[Rcpp::export]]
Rcpp::List bc_fun_cpp(const arma::vec &p_value, double alpha = 0.05)
{
  // initialization
  int n = p_value.n_elem;
  arma::uvec index_select;
  double e_tmp = 0;
  double tau = 0;
  arma::vec evalue = arma::zeros(n);
  // order p-value
  double p_tmp, n_num, n_dom, hat_fdp;
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
        tau = p_tmp;
        index_select = arma::find(p_value <= p_tmp);
        e_tmp = (double)n / n_num;
        evalue.elem(index_select).fill(e_tmp);
        break;
      }
    }
  }
  // return results
  return Rcpp::List::create(Named("evalue") = evalue, Named("e_tmp") = e_tmp,
                            Named("tau_thres") = tau, Named("n_num") = n_num,
                            Named("index_select") = index_select);
}

// eBH method
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

// BH BC average
// [[Rcpp::export]]
Rcpp::List bhbc_fix(const arma::vec &p_value, double alpha_bh = 0.025,
                    double alpha_ebh = 0.05, double wbh = 0.5, double wbc = 0.5)
{
  // BH method
  Rcpp::List bh_res = bh_fun_cpp(p_value, alpha_bh);
  arma::vec evalue_bh = bh_res["evalue"];
  // BC method
  Rcpp::List bc_res = bc_fun_cpp(p_value, alpha_bh);
  arma::vec evalue_bc = bc_res["evalue"];
  // generate evalue
  arma::vec evalue = wbh * evalue_bh + wbc * evalue_bc;
  // eBH method
  arma::uvec index_select = ebh_fun_cpp(evalue, alpha_ebh);
  return Rcpp::List::create(Named("evalue") = evalue,
                            Named("index_select") = index_select);
}

// BH BC ada
// [[Rcpp::export]]
Rcpp::List bhbc_ada(const arma::vec &p_value, double alpha_bh = 0.05,
                    double alpha_ebh = 0.05)
{
  //// initialization
  int n = p_value.n_elem;
  // BH method
  Rcpp::List bh_res = bh_fun_cpp(p_value, alpha_bh);
  arma::vec evalue_bh = bh_res["evalue"];
  arma::uvec index_bh = bh_res["index_select"];
  // BC method
  Rcpp::List bc_res = bc_fun_cpp(p_value, alpha_bh);
  arma::vec evalue_bc = bc_res["evalue"];
  arma::uvec index_bc = bc_res["index_select"];
  double n_num = bc_res["n_num"];
  double ave_wbc = (double)n_num / n;
  // double tau_bc = bc_res["tau_thres"];

  //// calculate weight
  // tilde p vector
  arma::vec tilde_p = arma::zeros(n);
  for (int iter_n = 0; iter_n < n; iter_n++)
  {
    tilde_p(iter_n) = std::min(p_value(iter_n), 1.0 - p_value(iter_n));
  }
  // calculate threshold BH
  arma::vec tau_bh_weight(n);
  arma::vec tildep_tmp;
  Rcpp::List bhres_tmp;
  for (int iter_i = 0; iter_i < n; iter_i++)
  {
    tildep_tmp = tilde_p;
    tildep_tmp(iter_i) = 0;
    bhres_tmp = bh_fun_cpp(tildep_tmp, alpha_bh);
    tau_bh_weight(iter_i) = bhres_tmp["tau_thres"];
  }
  double tau_bh_max = arma::max(tau_bh_weight);
  // weight for BH
  int n_bh = index_bh.size();
  arma::vec w_bh = arma::zeros(n);
  if (n_bh > 0)
  {
    Rcpp::List bcres_tmp;
    double ave_wbh_tmp;
    double w_bh_tmp;
    for (int iter_bh = 0; iter_bh < n_bh; iter_bh++)
    {
      int iter_i = index_bh(iter_bh);
      arma::vec tau_bc_weight = arma::zeros(n);
      for (int iter_j = 0; iter_j < n; iter_j++)
      {
        if (iter_j == iter_i)
        {
          tau_bc_weight(iter_j) = 0;
        }
        else
        {
          tildep_tmp = p_value;
          tildep_tmp(iter_j) = tilde_p(iter_j);
          tildep_tmp(iter_i) = 0;
          bcres_tmp = bc_fun_cpp(tildep_tmp, alpha_bh);
          tau_bc_weight(iter_j) = bcres_tmp["tau_thres"];
        }
      }
      // BH weight
      ave_wbh_tmp = (1 + sum(p_value >= 1 - tau_bc_weight)) / n;
      w_bh_tmp = tau_bh_weight(iter_i) / (tau_bh_weight(iter_i) + ave_wbh_tmp);
      if (ISNAN(w_bh_tmp))
      {
        w_bh_tmp = 1;
      }
      w_bh(iter_i) = w_bh_tmp;
    }
  }
  // weight for BC
  int n_bc = index_bc.size();
  double w_bc = 0;
  if (n_bc > 0)
  {
    w_bc = ave_wbc / (ave_wbc + tau_bh_max);
  }
  // generate evalue
  arma::vec evalue = w_bh % evalue_bh + w_bc * evalue_bc;
  // eBH method
  arma::uvec index_select = ebh_fun_cpp(evalue, alpha_ebh);
  return Rcpp::List::create(Named("evalue") = evalue,
                            Named("index_select") = index_select);
}

// fast BH BC ada
// [[Rcpp::export]]
Rcpp::List fastbhbc_ada(const arma::vec &p_value, double alpha_bh = 0.05,
                        double alpha_ebh = 0.05)
{
  //// initialization
  int n = p_value.n_elem;
  // BH method
  Rcpp::List bh_res = bh_fun_cpp(p_value, alpha_bh);
  arma::vec evalue_bh = bh_res["evalue"];
  arma::uvec index_bh = bh_res["index_select"];
  // BC method
  Rcpp::List bc_res = bc_fun_cpp(p_value, alpha_bh);
  arma::vec evalue_bc = bc_res["evalue"];
  arma::uvec index_bc = bc_res["index_select"];
  // double tau_bc = bc_res["tau_thres"];
  double n_num = bc_res["n_num"];
  double ave_wbc = (double)n_num / n;

  //// calculate weight
  // tilde p vector
  arma::vec tilde_p = arma::zeros(n);
  for (int iter_n = 0; iter_n < n; iter_n++)
  {
    tilde_p(iter_n) = std::min(p_value(iter_n), 1.0 - p_value(iter_n));
  }
  // calculate threshold BH
  arma::vec tau_bh_weight(n);
  arma::vec tildep_tmp;
  Rcpp::List bhres_tmp;
  for (int iter_i = 0; iter_i < n; iter_i++)
  {
    tildep_tmp = tilde_p;
    tildep_tmp(iter_i) = 0;
    bhres_tmp = bh_fun_cpp(tildep_tmp, alpha_bh);
    tau_bh_weight(iter_i) = bhres_tmp["tau_thres"];
  }
  double tau_bh_max = arma::max(tau_bh_weight);
  // calculate threshold BC
  arma::vec tau_bc_weight(n);
  Rcpp::List bcres_tmp;
  for (int iter_j = 0; iter_j < n; iter_j++)
  {
    tildep_tmp = p_value;
    tildep_tmp(iter_j) = tilde_p(iter_j);
    bcres_tmp = bc_fun_cpp(tildep_tmp, alpha_bh);
    tau_bc_weight(iter_j) = bcres_tmp["tau_thres"];
  }

  // weight for BH
  int n_bh = index_bh.size();
  arma::vec w_bh = arma::zeros(n);
  if (n_bh > 0)
  {
    double ave_wbh_tmp;
    double w_bh_tmp;
    arma::vec tau_bc_weight_tmp = arma::zeros(n);
    for (int iter_bh = 0; iter_bh < n_bh; iter_bh++)
    {
      int iter_i = index_bh(iter_bh);
      tau_bc_weight_tmp = tau_bc_weight;
      tau_bc_weight_tmp(iter_i) = 1;
      // BH weight
      ave_wbh_tmp = sum(p_value >= 1 - tau_bc_weight_tmp) / n;
      w_bh_tmp = tau_bh_weight(iter_i) / (tau_bh_weight(iter_i) + ave_wbh_tmp);
      if (ISNAN(w_bh_tmp))
      {
        w_bh_tmp = 1;
      }
      w_bh(iter_i) = w_bh_tmp;
    }
  }
  // weight for BC
  int n_bc = index_bc.size();
  double w_bc = 0;
  if (n_bc > 0)
  {
    w_bc = ave_wbc / (ave_wbc + tau_bh_max);
  }
  // generate evalue
  arma::vec evalue = w_bh % evalue_bh + w_bc * evalue_bc;
  // eBH method
  arma::uvec index_select = ebh_fun_cpp(evalue, alpha_ebh);
  return Rcpp::List::create(Named("evalue") = evalue,
                            Named("index_select") = index_select);
}

// // fast BH BC rand
// // [[Rcpp::export]]
// Rcpp::List bhbc_rand(const arma::vec &p_value, double alpha_ebh = 0.05)
// {
//   //// initialization
//   int n = p_value.n_elem;
//   // BH method
//   Rcpp::List bh_res = bh_fun_cpp(p_value, alpha_ebh);
//   arma::vec evalue_bh = bh_res["evalue"];
//   arma::uvec index_bh = bh_res["index_select"];
//   // BC method
//   Rcpp::List bc_res = bc_fun_cpp(p_value, alpha_ebh);
//   arma::vec evalue_bc = bc_res["evalue"];
//   arma::uvec index_bc = bc_res["index_select"];
//   // random select
//   double prob = 0.5;
//   double w_bh = Rcpp::rbinom(1, 1, prob)(0);
//   arma::vec evalue = w_bh * evalue_bh + (1.0 - w_bh) * evalue_bc;
//   // eBH method
//   arma::uvec index_select = ebh_fun_cpp(evalue, alpha_ebh);
//   return Rcpp::List::create(Named("evalue") = evalue,
//                             Named("index_select") = index_select);
// }