

plot_diagnostics = function(fit){


  table = flextable::flextable(data.frame(fit$diagnostic_summary()))
  traces = list(mcmc_trace(fit$draws("gm")),mcmc_trace(fit$draws("tau_u")))

  return(list(table,traces))
}
