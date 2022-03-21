function find_periodogram_peak(ode_solution, t_eval, f_sampling)
  nodes = size(ode_solution, 1)
  peaks = fill(Float64[], nodes)
  for i in 1:nodes
    # Remove mean from solution to avoid peak at low frequencies
    signal = ode_solution(t_eval)[i,:] .- mean(ode_solution(t_eval)[i,:])
    pgram = periodogram(signal; fs=f_sampling)
    pks, vals = findmaxima(pgram.power, 3)
    # Save largest peak in the periodogram
    peaks[i] = [pks[argmax(vals)], vals[argmax(vals)]]
  end
  return peaks
end
