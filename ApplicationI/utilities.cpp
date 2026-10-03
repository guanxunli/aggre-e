#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

////////////////////// BC function //////////////////////
// [[Rcpp::export]]
Rcpp::List bc_fun_cpp(const arma::vec &p_value, double alpha = 0.05)
{
  // initialization
  int n = p_value.n_elem;
  if (!p_value.is_finite() || arma::any(p_value < 0.0) || arma::any(p_value > 1.0))
  {
    Rcpp::stop("p-values must be finite and lie in [0, 1]");
  }
  if (!(alpha > 0.0 && alpha < 1.0))
  {
    Rcpp::stop("alpha must lie in (0, 1)");
  }
  arma::uvec index_select;
  double e_tmp = 0;
  double thres = 0;
  double num = 1.0 + arma::sum(p_value >= 1.0);
  arma::vec evalue = arma::zeros(n);
  // order p-value
  double p_tmp, dom, hat_fdp;
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
      if (p_tmp >= 0.5)
      {
        continue;
      }
      num = 1 + arma::sum(p_value >= (1 - p_tmp));
      dom = arma::sum(p_value <= p_tmp);
      hat_fdp = (double)num / dom;
      if (hat_fdp <= alpha)
      {
        thres = p_tmp;
        index_select = arma::find(p_value <= p_tmp);
        e_tmp = (double)n / num;
        evalue.elem(index_select).fill(e_tmp);
        break;
      }
    }
  }
  // return results
  return Rcpp::List::create(Named("evalue") = evalue, Named("e_tmp") = e_tmp,
                            Named("thres") = thres, Named("num") = num,
                            Named("index_select") = index_select);
}

//////////////////////// eBH method //////////////////////
// [[Rcpp::export]]
arma::uvec ebh_fun_cpp(const arma::vec &e_value, double alpha = 0.05)
{
  // initialization
  int n = e_value.n_elem;
  if (!e_value.is_finite() || arma::any(e_value < 0.0))
  {
    Rcpp::stop("e-values must be finite and nonnegative");
  }
  if (!(alpha > 0.0 && alpha < 1.0))
  {
    Rcpp::stop("alpha must lie in (0, 1)");
  }
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
      if (e_tmp >= thres)
      {
        index_select = arma::find(e_value >= thres);
        break;
      }
    }
  }
  // return results
  return index_select;
}

//// data dependent ebh
// [[Rcpp::export]]
Rcpp::List ebh_apa_cpp(Rcpp::List pvalue_list,
                       double alpha_bc = 0.05, double alpha_ebh = 0.05)
{
  //// initialization
  int n_group = pvalue_list.length();
  arma::vec n_vec(n_group), num_vec(n_group), n_rej_vec(n_group), thres_vec(n_group);
  Rcpp::List evalue_list(n_group), index_select_list(n_group);
  // separate BC method
  Rcpp::List bcres;
  for (int iter_group = 0; iter_group < n_group; ++iter_group)
  {
    arma::vec pvalue_tmp = pvalue_list[iter_group];
    n_vec(iter_group) = pvalue_tmp.n_elem;
    bcres = bc_fun_cpp(pvalue_tmp, alpha_bc);
    num_vec(iter_group) = bcres["num"];
    evalue_list[iter_group] = bcres["evalue"];
    arma::uvec index_select_tmp = bcres["index_select"];
    thres_vec(iter_group) = bcres["thres"];
    index_select_list[iter_group] = index_select_tmp;
    n_rej_vec(iter_group) = index_select_tmp.size();
  }
  double n = arma::sum(n_vec);

  //// calculate weight vector
  arma::vec dom_vec(n_group);
  Rcpp::List bcres_tmp;
  double tau_tmp, n_tmp;
  for (int iter_group = 0; iter_group < n_group; ++iter_group)
  {
    n_tmp = n_vec(iter_group);
    arma::vec tauvec_tmp = arma::zeros(n_tmp);
    arma::vec pvalue_tmp = pvalue_list[iter_group];
    arma::vec puse_tmp(n_tmp);
    for (int iter_n = 0; iter_n < n_tmp; ++iter_n)
    {
      if (pvalue_tmp(iter_n) <= 0.5)
      {
        tau_tmp = thres_vec(iter_group);
      }
      else
      {
        puse_tmp = pvalue_tmp;
        puse_tmp(iter_n) = 1 - pvalue_tmp(iter_n);
        bcres_tmp = bc_fun_cpp(puse_tmp, alpha_bc);
        tau_tmp = bcres_tmp["thres"];
      }
      if (pvalue_tmp(iter_n) >= 1 - tau_tmp)
      {
        tauvec_tmp(iter_n) = 1;
      }
    }
    dom_vec(iter_group) = arma::sum(tauvec_tmp);
  }
  double dom_sum = arma::sum(dom_vec);
  // update evalue
  double w, num, dom, num_tmp;
  for (int iter_group = 0; iter_group < n_group; ++iter_group)
  {
    if (n_rej_vec(iter_group) > 0)
    {
      n_tmp = n_vec(iter_group);
      num_tmp = num_vec(iter_group);
      num = num_tmp * n / n_tmp;
      dom = num_tmp + dom_sum - dom_vec(iter_group);
      w = num / dom;
      arma::vec evalue_tmp = evalue_list[iter_group];
      evalue_list[iter_group] = w * evalue_tmp;
    }
  }

  // generate evalue
  arma::vec evalue(n);
  int index_begin = 0;
  int index_end;
  for (int iter_group = 0; iter_group < n_group; ++iter_group)
  {
    index_end = index_begin + n_vec(iter_group) - 1;
    arma::vec evalue_tmp = evalue_list[iter_group];
    arma::uvec index_use = arma::regspace<arma::uvec>(index_begin, index_end);
    evalue.elem(index_use) = evalue_tmp;
    index_begin = index_end + 1;
  }
  // eBH method
  arma::uvec index_select = ebh_fun_cpp(evalue, alpha_ebh);
  return Rcpp::List::create(Named("evalue") = evalue,
                            Named("index_select") = index_select);
}
