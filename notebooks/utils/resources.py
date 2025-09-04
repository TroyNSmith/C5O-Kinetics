def optimization_pars(num_heavy_ats: int, label: str):
    # [0 heavy atoms, 1 heavy atom, 2 heavy atoms, ...]
    ccsdt_cpus = [0, 4, 4, 12, 12, 16, 20, 20, 20]
    ccsdt_mem = [0, 1, 1, 2, 4, 4, 4, 6, 12]
    ccsdt_lscratch = [0, 50, 50, 100, 200, 200, 200, 200, 200]
    ccsdt_time = [
        "00:00:00",
        "00:30:00",
        "00:30:00",
        "00:30:00",
        "00:60:00",
        "02:00:00",
        "02:00:00",
        "02:00:00",
        "06:00:00",
    ]

    if "XTB" in label:
        pars = {
            "num_cpus": 8,
            "mem_per_cpu": 1,
            "lscratch_size": 10,
            "time": "00:30:00",
        }
    elif "REV" in label:
        pars = {
            "num_cpus": 16,
            "mem_per_cpu": 1,
            "lscratch_size": 20,
            "time": "02:00:00",
        }
    elif "CCS" in label:
        pars = {
            "num_cpus": ccsdt_cpus[num_heavy_ats],
            "mem_per_cpu": ccsdt_mem[num_heavy_ats],
            "lscratch_size": ccsdt_lscratch[num_heavy_ats],
            "time": ccsdt_time[num_heavy_ats],
        }
    return pars
