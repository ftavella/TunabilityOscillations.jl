function find_periodogram_peak(ode_solution, t_eval)
  nodes = size(ode_solution, 1)
  tpoints = size(t_eval, 1)
  peaks = zeros(nodes)
  for i in nodes
    # Remove mean from solution to avoid peak at low frequencies
    # TODO
    #=
    signal = ode_solution(t_eval)[i,:] .- mean(ode_solution(t_eval)[i,:])
    pks, vals = findmaxima(signal, 5)
    pgram = periodogram(signal; fs=tpoints)
    pks, vals = findmaxima(pgram.power, 3)
    # Save largest peak
    peaks[i] = [pks[argmax(vals)], vals[argmax(vals)]]
    =#
  end
  return peaks
end
